#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/// __________________________________________________
/// __________________________________________________
/// Class and Function declarations

void _ctrsgn(const type_info&, bool);

// Formula simplification
inline double mod1by60(double);
inline double mod1e2(double);
inline double round2(double, int);
inline double polish(double);

// Utility
template<class T, class U> 
inline vector<U> get_vec_attr(const T&, const char*);
template<class T>
inline int get_fmt_attribute(const T&);
template<class T>
inline void checkinherits(T&, const char*);
template<class T>
inline bool is_item_in_obj(const T, const int);
template<class T>
inline void prefixvecstr(vector<string>&, const vector<T>&);
inline bool prefixwithnames(vector<string>&, RObject&);
inline string str_tolower(string);
template<class T>
int nameinobj(const T, const char*);

//CoordType
enum class CoordType : char { decdeg, degmin, degminsec };

inline const CoordType get_coordtype(const int);
template<class T>
inline const CoordType get_coordtype(const T&);
inline const int coordtype_to_int(CoordType);

inline string cardpoint(bool, bool);
inline string cardi_b(bool);

//FamousFive
struct FamousFive;
struct FF_decdeg;
struct FF_degmin;
struct FF_degminsec;

//Convertor
template<CoordType type>
class Convertor;

//CoordType switch
template<class T>
vector<string> format_switch(const T&, CoordType);

//Coord
class Coord;

class Validator;

template<CoordType type>
class Format;

template<class T, CoordType type>
class FormatLL;

ostream& operator<<(ostream&, const Coord&);

// WayPoint
class WayPoint;

ostream& operator<<(ostream&, const WayPoint&);

// validation
bool check_valid(const NumericVector);
bool check_valid(const DataFrame);

template<class T>
bool validated(T, const char*, bool&);

template<class T, class U>
const T revalidate(const T);

constexpr auto revalid_Coord = &revalidate<NumericVector, Coord>;
constexpr auto revalid_WayPoint = &revalidate<DataFrame, WayPoint>;

template<class T, class U>
inline const T validate(const T);

bool valid_ll(const DataFrame);

// Conversion
template<class T, class U>
void convene(T, CoordType newtype);

template<CoordType newtype>
inline void convert(CoordType, NumericVector);
template<CoordType newtype>
inline void convert(CoordType, DataFrame);

// exported
NumericVector coords(NumericVector, const int);
NumericVector coords_replace(NumericVector, int);
NumericVector latlon(NumericVector, LogicalVector&);
NumericVector validatecoord(NumericVector);
CharacterVector formatcoord(NumericVector);
NumericVector printcoord(NumericVector);
DataFrame waypoints(DataFrame, int);
DataFrame waypoints_replace(DataFrame, int);
DataFrame validatewaypoint(DataFrame);
CharacterVector formatwaypoint(NumericVector);
DataFrame printwaypoint(DataFrame);
NumericVector as_coord(DataFrame, bool);


/// __________________________________________________
/// __________________________________________________
/// Development and Debugging functions

/// Report object construction and destruction
void _ctrsgn(const type_info& obj, bool destruct = false)
{
//	cout << (destruct ? "Destroying " : "Constructing ") << obj.name() << endl;
}


/// __________________________________________________
/// __________________________________________________
/// Formula simplification functions

/// __________________________________________________
/// Multiply integer part by sixty
inline double mod1by60(double x)
{
	return fmod(x, 1) * 60;
}


/// __________________________________________________
/// Modulus after multiplication by 100
inline double mod1e2(double x)
{
	return fmod(x, 1e2);
}


/// __________________________________________________
/// Round a floating point number to n dp
inline double round2(double x, int n = 2)
{
	int pow10n = pow(10, n);
	return round(x * pow10n) / pow10n;
}


/// __________________________________________________
/// Round a floating point number to 10 dp
inline double polish(double x)
{
	return round(x * 1e10) / 1e10;
}


/// __________________________________________________
/// __________________________________________________
/// Utility functions

/// __________________________________________________
/// Return named attribute as vector<U> or empty vector<U>
template<class T, class U> 
inline vector<U> get_vec_attr(const T& t, const char* attrname)
{
//	cout << "@get_vec_attr<T, U>(const T&, const char*) attr \"" << attrname << "\" " << boolalpha << t.hasAttribute(attrname) << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	return t.hasAttribute(attrname) ? as<vector<U>>(t.attr(attrname)) : vector<U>();
}


/// __________________________________________________
/// Return "fmt" attribute as int
template<class T>
inline int get_fmt_attribute(const T& t)
{
//	cout << "@get_fmt_attribute<T>(const T&) " << as<int>(t.attr("fmt")) << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	return as<int>(t.attr("fmt"));
}


/// __________________________________________________
/// Does object inherit given class?
template<class T>
inline void checkinherits(T& t, const char* classname)
{
//	cout << "@checkinherits<T>(T& t, const char* classname) t " << typeid(t).name() << " classname " << classname << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	if (!t.inherits(classname)) stop("Argument must be a \"%s\" object", classname);
}


/// __________________________________________________
/// Is item number present in object? (Using C++ numbering)
template<class T>
inline bool is_item_in_obj(const T t, const int item)
{
//	cout << "@is_item_in_obj(T, int)\n";
	if (NA_INTEGER == item)
		return false;
	else
		return !(item < 0) && item < t.size();
}


/// __________________________________________________
/// Prefix vector<string> elements with elements of vector<T>—default for vector<string> prefix
template<class T>
inline void prefixvecstr(vector<string>& sv, const vector<T>& prefix)
{
//	cout << "@prefixvecstr<T>(vector<string>&, const vector<T>&)\n";
	transform(sv.begin(), sv.end(), prefix.begin(), sv.begin(), [](string& lls, const string& name) { return name + "  " + lls; });	
}


/// __________________________________________________
/// Specialisation for vector<int> prefix
template<>
inline void prefixvecstr(vector<string>& sv, const vector<int>& prefix)
{
//	cout << "@prefixvecstr<>(vector<string>&, const vector<int>&)\n";
	transform(sv.begin(), sv.end(), prefix.begin(), sv.begin(), [](string& lls, const int name) { return to_string(name) + "  " + lls; });	
}


/// __________________________________________________
/// Prefix vector<string> elements with elements of RObject 
inline bool prefixwithnames(vector<string>& sv, RObject& namesobj)
{
//	cout << "@prefixwithnames(vector<string>&, RObject&)\n";
	if (is<CharacterVector>(namesobj))
		prefixvecstr(sv, as<vector<string>>(namesobj));
	else if(is<IntegerVector>(namesobj))
		prefixvecstr(sv, as<vector<int>>(namesobj));
	else
		return false;
	return true;
}


/// __________________________________________________
/// string to lower case (see cppreference.com std::tolower)
inline string str_tolower(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return tolower(c); });
    return s;
}


/// __________________________________________________
/// Find position of name within object names
template<class T>
int nameinobj(const T t, const char* name)
{
//	cout << "@nameinobj<T>(const T, const char*) name is " << name << endl;
	static_assert(std::is_same<List, T>::value || std::is_same<DataFrame, T>::value, "T must be List or DataFrame");
	vector<string> names { get_vec_attr<T, string>(t, "names") };
	if (!names.size())
		return -1;
	int i = 0;
	for (auto str : names ) {
//		cout << "@nameinobj<T>(const T, const char*) testing " << str << endl;
		if (!str_tolower(str).compare(name)) {
//			cout << "@nameinobj<T>(const T, const char*) found " << str << endl;
			break;
		}
		i++;
	}
	if (i == names.size())
		i = -1;
	return i;
}


/// __________________________________________________
/// __________________________________________________
/// CoordType enum class

/// __________________________________________________
/// Convert int to CoordType enum
inline const CoordType get_coordtype(const int i)
{
//	cout << "@get_coordtype(int) " << i << endl;
	if (i < 1 || i > 3)
		stop("\"newfmt\" must be between 1 and 3");
	return vector<CoordType>{ CoordType::decdeg, CoordType::degmin, CoordType::degminsec }[i - 1];
}


/// __________________________________________________
/// Convert "fmt" attribute to CoordType enum
template<class T>
inline const CoordType get_coordtype(const T& t)
{
//	cout << "@get_coordtype<T>(const T&) " << get_fmt_attribute(t) << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	return get_coordtype(get_fmt_attribute(t));
}


/// __________________________________________________
/// Convert CoordType enum to int
inline const int coordtype_to_int(CoordType ct)
{
//	cout << "@coordtype_to_int(CoordType ct) " << static_cast<char>(ct) + 1 << endl;
	return static_cast<char>(ct);
}


/// __________________________________________________
/// __________________________________________________
/// Cardinal points of direction
inline string cardpoint(bool negative, bool lat)
{
	return negative ? (lat ? " S" : " W") : (lat ? " N" : " E") ;
}


/// __________________________________________________
/// Cardinal points without "latlon" attribute
inline string cardi_b(bool negative)
{
	return negative ? " (S/W)" : " (N/E)";
}


/// __________________________________________________
/// __________________________________________________
/// FamousFive Class and Derived Classes

struct FamousFive {
//	FamousFive() { cout << "§FamousFive() "; _ctrsgn(typeid(*this)); }
	virtual ~FamousFive() = 0;	
	virtual int get_deg(double x) const = 0;
	virtual double get_decdeg(double x) const = 0;
	virtual int get_min(double x) const = 0;
	virtual double get_decmin(double x) const = 0;
	virtual double get_sec(double x) const = 0;
};

FamousFive::~FamousFive()
{
//	cout << "§~FamousFive(CoordType) "; _ctrsgn(typeid(*this), true); 
}	

/// __________________________________________________
/// Derived class for decimal degrees	
struct FF_decdeg : public FamousFive {
//	FF_decdeg() { cout << "§FF_decdeg() "; _ctrsgn(typeid(*this)); }	
	~FF_decdeg() = default;
//	~FF_decdeg() { cout << "§~FF_decdeg::FF_decdeg() "; _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x); }
	double get_decdeg(double x) const { return x; }
	int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
	double get_decmin(double x) const { return polish(mod1by60(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
} ff_decdeg;

/// __________________________________________________
/// Derived class for degrees and minutes
struct FF_degmin : public FamousFive {
//	FF_degmin() { cout << "§FF_degmin() "; _ctrsgn(typeid(*this)); }	
	~FF_degmin() = default;
//	~FF_degmin() { cout << "§~FF_degmin::FF_degmin() "; _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x / 1e2); }
	double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
	int get_min(double x) const { return int(x) % int(1e2); }
	double get_decmin(double x) const { return polish(mod1e2(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
} ff_degmin;

/// __________________________________________________
/// Derived class for degrees, minutes and seconds
struct FF_degminsec : public FamousFive {
//	FF_degminsec() { cout << "§FF_degminsec() "; _ctrsgn(typeid(*this)); }	
	~FF_degminsec() = default;
//	~FF_degminsec() { cout << "§~FF_degminsec::FF_degminsec() "; _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x / 1e4); }
	double get_decdeg(double x) const { return int(x / 1e4) + (double)int(fmod(x, 1e4) / 1e2) / 60 + mod1e2(x) / 3600; }
	int get_min(double x) const { return (int(x) % int(1e4)) / 1e2; }
	double get_decmin(double x) const { return int(fmod(x, 1e4) / 1e2) + mod1e2(x) / 60; }
	double get_sec(double x) const { return mod1e2(x); }
} ff_degminsec;

vector<FamousFive*> vff { &ff_decdeg, &ff_degmin, &ff_degminsec };


/// __________________________________________________
/// __________________________________________________
/// Templated coord type conversion functors

template<CoordType type>
class Convertor {
	protected:
		const FamousFive& ff; 
	public:
		Convertor(const FamousFive& _ff) : ff(_ff)
		{
//			cout << "§Convertor<CoordType>::Convertor(const FamousFive&) "; _ctrsgn(typeid(*this));
		}
		~Convertor() = default;
//		~Convertor() { cout << "§Convertor<type>::~Convertor() "; _ctrsgn(typeid(*this), true); }
		double operator()(double n);
};


/// __________________________________________________
/// Default operator(), for decimal degrees
template<CoordType type>
inline double Convertor<type>::operator()(double n)
{
//	cout << "@Convertor<CoordType>::operator() [default for CoordType::decdeg]\n";
	return ff.get_decdeg(n);
}


/// __________________________________________________
/// Specialised operator() for degrees and minutes
template<>
inline double Convertor<CoordType::degmin>::operator()(double n)
{
//	cout << "@Convertor<CoordType::degmin>::operator()\n";
	return ff.get_deg(n) * 1e2 + ff.get_decmin(n);
}


/// __________________________________________________
/// Specialised operator() for degrees, minutes and seconds
template<>
inline double Convertor<CoordType::degminsec>::operator()(double n)
{
//	cout << "@Convertor<CoordType::degminsec>::operator()\n";
	return ff.get_deg(n) * 1e4 + ff.get_min(n) * 1e2 + ff.get_sec(n);
}


/// __________________________________________________
/// __________________________________________________
/// Templated coord formatting functors

template<CoordType type>
class Format {
	protected:
		const FamousFive& ff;
		ostringstream outstrstr;
	public:
		Format(const FamousFive& _ff) : ff(_ff)
		{
//			cout << "§Format<CoordType>::Format(const FamousFive&) "; _ctrsgn(typeid(*this));
		}
		~Format() = default;
//		~Format() { cout << "§Format<CoordType>::~Format() "; _ctrsgn(typeid(*this), true); }
		string operator()(double n);
};

/// __________________________________________________
/// Default operator(), for decimal degrees
template<CoordType type>
inline string Format<type>::operator()(double n)
{
//	cout << "@Format<CoordType>::operator() [default for CoordType::decdeg]\n";
	outstrstr.str("");
	outstrstr << setw(11) << setfill(' ')  << fixed << setprecision(6) << ff.get_decdeg(n) << "\u00B0";
	return outstrstr.str();
}

/// __________________________________________________
/// Specialised operator() for degrees and minutes
template<>
inline string Format<CoordType::degmin>::operator()(double n)
{
//	cout << "@Format<CoordType::degmin>::operator()\n";
	outstrstr.str("");
	outstrstr << setw(3) << setfill(' ') << abs(ff.get_deg(n)) << "\u00B0"
					  << setw(7) << setfill('0') << fixed << setprecision(4) << abs(ff.get_decmin(n)) << "\u2032";
	return outstrstr.str();
}

/// __________________________________________________
/// Specialised operator() for degrees, minutes and seconds
template<>
inline string Format<CoordType::degminsec>::operator()(double n)
{
//	cout << "@Format<CoordType::degminsec>::operator()\n";
	outstrstr.str("");
	outstrstr << setw(3) << setfill(' ') << abs(ff.get_deg(n)) << "\u00B0"
					  << setw(2) << setfill('0') << abs(ff.get_min(n)) << "\u2032"
					  << setw(5) << fixed << setprecision(2) << abs(ff.get_sec(n)) << "\u2033";
	return outstrstr.str();
}


/// __________________________________________________
/// __________________________________________________
/// Formatting functors for latitude and longitude

/// Default functor for degrees, minutes (and seconds)
template<class T, CoordType type>
class FormatLL {
		const FamousFive& ff; 
		vector<bool>::const_iterator ll_it;
		const int ll_size;
	public:
		FormatLL(const FamousFive& _ff, const vector<bool>& ll) : ff(_ff), ll_it(ll.begin()), ll_size(ll.size())
		{
//			cout << "§FormatLL<T, CoordType>::FormatLL(const FamousFive&, vector<bool>&) "; _ctrsgn(typeid(*this));
			static_assert(std::is_same<Coord, T>::value || std::is_same<WayPoint, T>::value, "T must be Coord or WayPoint");
		}
		~FormatLL() = default;
//		~FormatLL() { cout << "§FormatLL<T, CoordType>::~FormatLL() "; _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
//			cout << "@FormatLL<T, CoordType>::operator(string, double) [default for CoordType::degmin and CoordType::degminsec]\n";
			return ostr += ll_size ? cardpoint(ff.get_decmin(n) < 0, ll_size > 1 ? *ll_it++ : *ll_it) : cardi_b(ff.get_decmin(n) < 0);
		}
};

/// __________________________________________________
/// Specialised functor for decimal degrees Coord
template<>
class FormatLL<Coord, CoordType::decdeg> {
		vector<bool>::const_iterator ll_it;
		const int ll_size;
	public:
		FormatLL(const FamousFive& _ff, const vector<bool>& ll) : ll_it(ll.begin()), ll_size(ll.size())
		{
//			cout << "§FormatLL<Coord, CoordType::decdeg>::FormatLL(const FamousFive&, vector<bool>&) "; _ctrsgn(typeid(*this));
		}
		~FormatLL() = default;
//		~FormatLL() { cout << "§FormatLL<Coord, CoordType::decdeg>::~FormatLL() "; _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
//			cout << "@FormatLL<Coord, CoordType::decdeg>::operator(string, double)\n";
			if (ll_size)
				return ostr += ((ll_size > 1 ? *ll_it++ : *ll_it) ? " lat" : " lon");
			else
				return ostr;
		}
};

/// __________________________________________________
/// Specialised functor for decimal degrees WayPoint
template<>
class FormatLL<WayPoint, CoordType::decdeg> {
		vector<bool>::const_iterator ll_it;
		const int ll_size;
	public:
		FormatLL(const FamousFive& _ff, const vector<bool>& ll) : ll_it(ll.begin()), ll_size(ll.size())
		{
//			cout << "§FormatLL<WayPoint, CoordType::decdeg>::FormatLL(const FamousFive&, vector<bool>&) "; _ctrsgn(typeid(*this));
		}
		~FormatLL() = default;
//		~FormatLL() { cout << "§FormatLL<WayPoint, CoordType::decdeg>::~FormatLL() "; _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
//			cout << "@FormatLL<WayPoint, CoordType::decdeg>::operator(string, double)\n";
			return ostr;
		}
};


/// __________________________________________________
/// __________________________________________________
/// Validate functor

class Validator {
		const FamousFive& ff;
		vector<bool>::const_iterator ll_it;
		const int ll_size;
	public:
		Validator(const FamousFive& _ff, const vector<bool>& ll) : ff(_ff), ll_it(ll.begin()), ll_size(ll.size())
		{
//			cout << "§Validator::Validator(const FamousFive&, vector<bool>&) "; _ctrsgn(typeid(*this));
		}
		~Validator() = default;
//		~Validator() { cout << "§Validator::~Validator() "; _ctrsgn(typeid(*this), true); }
		bool operator()(double n)
		{
//			cout << "@Validator() " << " validating: " << setw(9) << setfill(' ') << n << endl;
			return !((abs(ff.get_decdeg(n)) > (ll_size && (ll_size > 1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(abs(ff.get_decmin(n)) >= 60) ||
				(abs(ff.get_sec(n)) >= 60));
		}
};


/// __________________________________________________
/// Format coords or waypoints vector<string> CoordType switch 
template<class T>
vector<string> format_switch(const T& t, CoordType ct)
{
//	cout << "@format_switch<T>(const T&, CoordType) " << typeid(t).name() << endl;
	vector<string> sv; 
	switch (ct)
	{
		case CoordType::decdeg:
			sv = t.template format<CoordType::decdeg>();
			break;

		case CoordType::degmin:
			sv = t.template format<CoordType::degmin>();
			break;

		case CoordType::degminsec:
			sv = t.template format<CoordType::degminsec>();
			break;

		default:
			stop("format_switch(const T&, CoordType) my bad");
	}
	return sv;
}


/// __________________________________________________
/// __________________________________________________
/// Coord base class
class Coordbase {
	protected:
		CoordType ct;
		const FamousFive& ff;

	public:
		Coordbase(CoordType _ct);
		Coordbase(const Coordbase&) = delete;						// Disallow copying
		Coordbase& operator=(const Coordbase&) = delete;			//  ——— ditto ———
		Coordbase(Coordbase&&) = delete;							// Disallow transfer ownership
		Coordbase& operator=(Coordbase&&) = delete;					// Disallow moving
		virtual ~Coordbase() = 0;

		const FamousFive& get_ff() const;
};


Coordbase::Coordbase(CoordType _ct) :
	ct(_ct), ff(*vff[coordtype_to_int(ct)])
{
//	cout << "§Coordbase::Coordbase(CoordType) "; _ctrsgn(typeid(*this));
}


Coordbase::~Coordbase()
{
//	cout << "§Coordbase::~Coordbase() "; _ctrsgn(typeid(*this), true);
}


/// __________________________________________________
/// Get const reference to ff
inline const FamousFive& Coordbase::get_ff() const
{
//	cout << "@Coordbase::get_ff()\n";
	return ff;
}


/// __________________________________________________
/// Coordinate derived class
class Coord : public Coordbase {
	protected:
		const NumericVector nv;
		const vector<bool> valid { false };
		const vector<bool> latlon;

	public:
		Coord(CoordType, const NumericVector);
		~Coord() = default;
//		~Coord() { cout << "§Coord::~Coord() "; _ctrsgn(typeid(*this), true); }

		void validate(bool warn = true) const;
		template<CoordType type>
		vector<string> format() const;
		void print(ostream&) const;
};


Coord::Coord(CoordType ct, const NumericVector nv) :
	Coordbase(ct), nv(nv),
	latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") } //,
{
//	cout << "§Coord::Coord(CoordType, const NumericVector) "; _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Validate coords vector
void Coord::validate(bool warn) const
{
//	cout << "@Coord::validate() " << typeid(*this).name() << " latlon " << LogicalVector(wrap(latlon)) << endl;
	vector<bool>& non_const_valid { const_cast<vector<bool>&>(valid) };
	non_const_valid.assign(nv.size(), {false});
	transform(nv.begin(), nv.end(), non_const_valid.begin(), Validator(ff, latlon));
	if (all_of(valid.begin(), valid.end(), [](bool v) { return v;}))
		non_const_valid.assign({true});
	else
		if (warn)
			warning("Validation failed!");
	const_cast<NumericVector&>(nv).attr("valid") = valid;
}


/// __________________________________________________
/// Format coordinates as vector<string> of CoordType
template<CoordType type>
vector<string> Coord::format() const
{
//	cout << "@Coord::format<CoordType>() " << typeid(*this).name() << endl;
	vector<string> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), Format<type>(ff));
	transform(out.begin(), out.end(), nv.begin(), out.begin(), FormatLL<Coord, type>(ff, latlon));
	return out;
}


/// __________________________________________________
/// Print coords vector
void Coord::print(ostream& stream) const
{
//	cout << "@Coord::print() " << typeid(*this).name() << endl;
	vector<string>&& sv = format_switch(*this, ct);
	vector<string> names { get_vec_attr<NumericVector, string>(nv, "names") };
	int strwdth = 0;
	if (names.size()) {
		prefixvecstr(sv, names);
		strwdth = max_element(sv.begin(), sv.end(), [](const string& a, const string& b){ return a.size() < b.size(); })->size();
	}
	for_each(sv.begin(), sv.end(), [&stream, strwdth](const string& s) { stream << setw(strwdth) << s << "\n"; });
}


/// __________________________________________________
/// Output Coord derived object to ostream
ostream& operator<<(ostream& stream, const Coord& c)
{
//	cout << "@operator<<(ostream&, const Coord&)\n";
	c.print(stream);
	return stream;
}


/// __________________________________________________
/// __________________________________________________
/// Waypoint class

class WayPoint : public Coordbase {
	protected:
		const DataFrame df;
		const NumericVector nvlat;
		const NumericVector nvlon;
		const vector<bool> validlat { false };
		const vector<bool> validlon { false };
	public:
		explicit WayPoint(CoordType, const DataFrame);
		~WayPoint() = default;
//		~WayPoint() { cout << "§WayPoint::~WayPoint() "; _ctrsgn(typeid(*this), true); }

		void validate(bool = true) const;
		template<CoordType type>
		vector<string> format() const;
		void print(ostream& stream) const;
};


WayPoint::WayPoint(CoordType ct, const DataFrame df) :
	Coordbase(ct), df(df),
	nvlat(df[get_vec_attr<DataFrame, int>(df, "llcols")[0] - 1]), 
	nvlon(df[get_vec_attr<DataFrame, int>(df, "llcols")[1] - 1])
{
//	cout << "§WayPoint::WayPoint(CoordType ct, const DataFrame) "; _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Validate WayPoint
void WayPoint::validate(bool warn) const
{
//	cout << "@WayPoint::validate(bool)\n";

	vector<bool>& non_const_validlat { const_cast<vector<bool>&>(validlat) };
	non_const_validlat.assign(nvlat.size(), {false});
	transform(nvlat.begin(), nvlat.end(), non_const_validlat.begin(), Validator(ff, vector<bool>{ true }));

	vector<bool>& non_const_validlon { const_cast<vector<bool>&>(validlon) };
	non_const_validlon.assign(nvlon.size(), {false});
	transform(nvlon.begin(), nvlon.end(), non_const_validlon.begin(), Validator(ff, vector<bool>{ false }));

	if (all_of(validlat.begin(), validlat.end(), [](bool v) { return v;}))
		non_const_validlat.assign({true});
	else
		if (warn)
			warning("Validation of latitude failed!");
	const_cast<DataFrame&>(df).attr("validlat") = validlat;

	if (all_of(validlon.begin(), validlon.end(), [](bool v) { return v;}))
		non_const_validlon.assign({true});
	else
		if (warn)
			warning("Validation of longitude failed!");
	const_cast<DataFrame&>(df).attr("validlon") = validlon;
}


/// __________________________________________________
/// Format waypoints as vector<string> of CoordType
template<CoordType type>
vector<string> WayPoint::format() const
{
//	cout << "@WayPoint::format()\n";
	vector<string> sv_lat(nvlat.size());
	transform(nvlat.begin(), nvlat.end(), sv_lat.begin(), Format<type>(ff));
	transform(sv_lat.begin(), sv_lat.end(), nvlat.begin(), sv_lat.begin(), FormatLL<WayPoint, type>(ff, vector<bool>{ true }));
	vector<string> sv_lon(nvlon.size());
	transform(nvlon.begin(), nvlon.end(), sv_lon.begin(), Format<type>(ff));
	transform(sv_lon.begin(), sv_lon.end(), nvlon.begin(), sv_lon.begin(), FormatLL<WayPoint, type>(ff, vector<bool>{ false }));

	vector<string> out(sv_lat.size());
	transform(sv_lat.begin(), sv_lat.end(), sv_lon.begin(), out.begin(), [](string& latstr, string& lonstr) { return latstr + "  " + lonstr; });
	return out;
}


/// __________________________________________________
/// Print WayPoint
void WayPoint::print(ostream& stream) const
{
//	cout << "@WayPoint::print() " << typeid(*this).name() << endl;
	const int i { coordtype_to_int(ct) };
	int spacing[] {  5,  7,  8,
					11, 13, 14 };
	ostringstream ostrstr;
	vector<string> ttlvec;
	ostrstr << "Latitude" << string(spacing[i], ' ') << "Longitude ";
	ttlvec.push_back(ostrstr.str());

	ostrstr.str("");
	ostrstr	<< string(spacing[i + 3], '_') << string(2, ' ') << string(spacing[i + 3] + 1, '_');
	ttlvec.push_back(ostrstr.str());

	vector<string>&& sv = format_switch(*this, ct);

	vector<int> namescolvec { get_vec_attr<DataFrame, int>(df, "namescol") };
	if (1 == namescolvec.size()) {
		int namescol = namescolvec[0] - 1;
		if (is_item_in_obj(df, namescol)) {
			RObject names = df[namescol];
			if (!prefixwithnames(sv, names))
				stop("Invalid \"namescol\" attribute! (df[namescol] neither a CharacterVector nor IntegerVector)");
		} else
			stop("Invalid \"namescol\" attribute! (item not in object)");
	} else if (df.hasAttribute("row.names")) {
		RObject rownames = df.attr("row.names");
		if (!prefixwithnames(sv, rownames))
			stop("Invalid \"row.names\" attribute! (neither a CharacterVector nor IntegerVector)");
	}

	int strwdth = max_element(sv.begin(), sv.end(), [](const string& a, const string& b){ return a.size() < b.size(); })->size();
	for_each(ttlvec.begin(), ttlvec.end(), [&stream, strwdth, i](const string& s) { int fudge[] = { 2, 6, 10 }; stream << setw(strwdth - fudge[i]) << s << "\n"; });
	for_each(sv.begin(), sv.end(), [&stream, strwdth](const string& s) { stream << setw(strwdth) << s << "\n"; });
}


/// __________________________________________________
/// Output WayPoint to ostream
ostream& operator<<(ostream& stream, const WayPoint& wp)
{
//	cout << "@operator<<(ostream&, const WayPoint&)\n";
	wp.print(stream);
	return stream;
}


/// __________________________________________________
/// __________________________________________________
/// Validation functions

/// __________________________________________________
/// Check "valid" attribute of NumericVector all true

bool check_valid(const NumericVector nv)
{
//	cout << "@check_valid(const NumericVector)" << endl;
	bool unvalidated = false;
	bool valid = validated(nv, "valid", unvalidated);
	if (unvalidated)
		revalid_Coord(nv);
	return valid;
}


/// __________________________________________________
/// Check "lat_valid" and "lon_valid attributes of DataFrame are all true
bool check_valid(const DataFrame df)
{
//	cout << "@check_valid(const DataFrame)\n";
	bool unvalidated = false;

	bool latvalid = validated(df, "validlat", unvalidated);
	if (unvalidated)
		return revalid_WayPoint(df);
	bool lonvalid = validated(df, "validlon", unvalidated);
	if (unvalidated)
		return revalid_WayPoint(df);

	if (!latvalid)
		warning("Invalid latitude!");
	if (!lonvalid)
		warning("Invalid longitude!");
	return latvalid || lonvalid;
}


/// __________________________________________________
/// Check NumericVector or DataFrame has been validated and valid vector attribute all true
template<class T>
bool validated(T t, const char* attrname, bool& unvalidated)
{
//	cout << "@validated<T>(T, const char*, bool&)" << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	const vector<bool>&& validvec = get_vec_attr<T, bool>(t, attrname);
	bool valid = all_of(validvec.begin(), validvec.end(), [](bool v) { return v;});
	unvalidated = (validvec.size()) ? false : true;
	return valid;
}


/// __________________________________________________
/// Revalidate NumericVector or DataFrame
template<class T, class U>
const T revalidate(const T t)
{
//	cout << "@revalidate<T, U>(const T)\n";
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	warning("Unvalidated %s! Revalidating…", typeid(t).name());
	validate<T, U>(t);	
	return check_valid(t);
}


/// __________________________________________________
/// Validate NumericVector or DataFrame
template<class T, class U>
inline const T validate(const T t)
{
//	cout << "@validate<T, U>(const T)\n";
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	U(get_coordtype(t), t).validate();
	return t;	
}


/// __________________________________________________
/// Check df has valid "llcols" attribute
bool valid_ll(const DataFrame df)
{
//	cout << "@valid_ll(const DataFrame)\n";
	bool valid = false;
	vector<int> llcols { get_vec_attr<DataFrame, int>(df, "llcols") };
	if (2 == llcols.size()) {
		transform(llcols.begin(), llcols.end(), llcols.begin(), [](int x){ return --x; });
		if (is_item_in_obj(df, llcols[0]) && is_item_in_obj(df, llcols[1]) && llcols[0] != llcols[1])
			if (is<NumericVector>(df[llcols[0]]) && is<NumericVector>(df[llcols[1]]))
				valid = true;
	}
	return valid;
}


/// __________________________________________________
/// __________________________________________________
/// Conversion functions

template<CoordType newtype>
inline void convert(CoordType type, NumericVector nv)
{
//	cout << "@convert<CoordType>(const Coord&, NumericVector) newtype " << coordtype_to_int(newtype) + 1 << endl;
	transform(nv.begin(), nv.end(), nv.begin(), Convertor<newtype>(*vff[coordtype_to_int(type)]));
}


template<CoordType newtype>
inline void convert(CoordType type, DataFrame df)
{
//	cout << "@convert<CoordType>(const WayPoint&, DataFrame) newtype " << coordtype_to_int(newtype) + 1 << endl;
	const vector<int> llcols { get_vec_attr<DataFrame, int>(df, "llcols") };
	NumericVector nvlat(df[llcols[0] - 1]);
	NumericVector nvlon(df[llcols[1] - 1]);
	transform(nvlat.begin(), nvlat.end(), nvlat.begin(), Convertor<newtype>(*vff[coordtype_to_int(type)]));
	transform(nvlon.begin(), nvlon.end(), nvlon.begin(), Convertor<newtype>(*vff[coordtype_to_int(type)]));
}


/// __________________________________________________
/// Convene NumericVector or DataFrame
template<class T, class U>
void convene(T t, CoordType newtype)
{
//	cout << "@convene<T&, U>(T, CoordType)\n";
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	CoordType type = get_coordtype(t);
	U(type, t).validate();

	if (type != newtype) {
		switch (newtype)
		{
			case CoordType::decdeg:
				convert<CoordType::decdeg>(type, t);
				break;

			case CoordType::degmin:
				convert<CoordType::degmin>(type, t);
				break;

			case CoordType::degminsec:
				convert<CoordType::degminsec>(type, t);
				break;

			default:
				stop("convene<T&, U>(const T&, U) my bad");
		}
		t.attr("fmt") = coordtype_to_int(newtype) + 1;
	}
}


/// __________________________________________________
/// __________________________________________________
/// Exported functions

/// __________________________________________________
/// Set R vector object class to coords and return,
/// or convert format of R coords object and return
//' @title Geographic or GPS Coordinate Class
//' 
//' @name coords
//' 
//' @description
//' \code{coords()} creates a robust representation of a series of geographic or GPS
//' coordinates instantiated as an object of class \code{"coords"}.
//' 
//' \code{coords()} and replacement form \verb{coords()<-} also convert the format of existing
//' objects of class \code{"coords"} between (i) decimal degrees, (ii) degrees and minutes, and
//' (iii) degrees, minutes and seconds.
//'
//' @details
//' Individual values provided in the \code{numeric} vector argument \code{nv} should have a decimal
//' point after the number of whole degrees in the case of \emph{decimal degrees}, after the number
//' of whole minutes in the case of \emph{degrees and minutes}, and after the number of whole
//' seconds in the case of \emph{degrees, minutes and seconds}.
//'
//' The \code{fmt} argument should be \code{1L} to represent decimal degrees, \code{2L} for degrees
//' and minutes, and \code{3L} for degrees, minutes and seconds and is used to provide both the
//' format of values in \code{numeric} vector argument \code{nv} to be converted into a
//' \code{"coords"} object and the desired format if a \code{"coords"} object is to be converted to
//' a new format. Note that on conversion of a \code{"coords"} object, the original \code{numeric}
//' vector argument \code{nv} is modified such that the values are as described in the previous
//' paragraph, and may be inspected using standard R code, see examples.
//'
//' The values of a newly created \code{"coords"} object are checked to ensure they are valid
//' geographic locations as described under \code{\link[=validate]{validate}()}. Likewise, a
//' check is made to ensure that an existing \code{"coords"} object to be converted to a new format
//' has already been validated; if not, it is re-validated. 
//'
//' @family coords_waypoints
//' @seealso
//' \code{\link[base:attr]{attr()}}, \code{\link[base:attributes]{attributes}},
//'   \code{\link[=latlon]{latlon}()}, \code{\link[base:numeric]{numeric()}} and
//'   \code{\link[=validate]{validate}()}.
//'
//' @param nv \code{numeric} vector of coordinate values, optionally named.
//'
//' @param fmt,value \code{integer}, 1L, 2L or 3L, indicating the current or desired coordinate
//'   format; default 1L.
//'
//' @param cd object of class \code{"coords"} created by function \code{\link[=coords]{coords}()}.
//'
//' @return
//' An object of class \code{"coords"}, comprising the original a \code{numeric} vector argument
//' \code{nv} with values possibly converted as appropriate and additional attributes: –
//' \item{\code{"fmt"}}{the coordinate format.}
//' \item{\code{"valid"}}{a \code{logical} vector indicating whether individual coordinate values
//'   are valid geographic locations.}
//'
//' @examples
//' ## Numeric vector representing degrees and minutes
//' dm <- c(5130.4659, 4932.7726, 4806.4339, 3853.3696, 0.0000, -3706.7044, -5306.2869, -2514.4093,
//'         -007.6754, 1823.9137, -12246.7203, -7702.1145, 0.0000, -1217.3178, 7331.0370, -5731.1536)
//'
//' ## Create a "coords" object of degrees and minutes
//' coords(dm, 2)
//'
//' ## Name the "coords" object
//' names(dm) <- rep(c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
//'                    "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"), 2)
//' dm
//'
//' ## Convert to degrees, minutes and seconds
//' coords(dm) <- 3
//' dm
//'
//' ## Convert to decimal degrees
//' coords(dm) <- 1
//' dm
//'
//' ## Decimal degrees as a standard numeric vector
//' as.numeric(dm)
//'
//' ## Convert to degrees and minutes, and format as a character vector
//' coords(dm) <- 2
//' (dm_chr <- format(dm))
//'
//' ## Output using {base} cat()
//' cat(dm_chr, fill = 18, labels = paste0("{#", 1:16, "}:"))
//'
//' rm(dm, dm_chr)
//'
// [[Rcpp::export]]
NumericVector coords(NumericVector nv, const int fmt = 1)
{
//	cout << "——Rcpp::export——coords(NumericVector)\n";
	CoordType newtype = get_coordtype(fmt);
	CoordType type;
	if (nv.inherits("coords")) {
		type = get_coordtype(nv);
//		cout <<  "@coords() nv is already a \"coords\" vector of type " << coordtype_to_int(type) + 1 << endl;
		if (newtype == type) {
//			cout << "——fmt out == fmt in!——" << endl;
			if (!check_valid(nv))
				stop("Invalid coords!");
			return nv;
		}
	} else {
		type = newtype;
		nv.attr("fmt") = fmt;
	}

	convene<NumericVector, Coord>(nv, newtype);
	nv.attr("class") = "coords";
	return nv;
}


/// __________________________________________________
/// coords() as replacement function
//' @rdname coords
// [[Rcpp::export(name = "`coords<-`")]]
NumericVector coords_replace(NumericVector nv, int value)
{
//	cout << "——Rcpp::export——`coords_replace(NumericVector, int)<-`\n";
	return coords(nv, value);
}


/// __________________________________________________
/// Set latlon attribute on "coords" NumericVector and revalidate
//' @title Latitude or Longitude Attribute for Coords
//' 
//' @name latlon
//'
//' @description \code{latlon()<-} adds the attribute \code{"latlon"} to objects of class
//' \code{\link[=coords]{"coords"}}, or modifies an existing \code{"latlon"} attribute.
//'
//' @details
//' Attribute \code{"latlon"} is a \code{logical} vector of length \code{1} or \code{length(cd)}
//' for which \code{TRUE} values represent latitude and \code{FALSE} values represent longitude.
//' Setting this attribute to any other length will result in an error. A \code{logical} vector of
//' length \code{1L} signifies that values are all latitude if \code{TRUE}, or all longitude if
//' \code{FALSE}.
//'
//' This attribute is used in formatting printed output and also by
//' \code{\link[=validate]{validate}()}. Indeed, the values of \code{cd} are revalidated every time
//' attribute \code{"latlon"} is added or changed.
//'
//' @seealso
//' \code{\link[=coords]{"coords"}}.
//'
//' @param cd object of class \code{\link[=coords]{"coords"}}.
//'
//' @param value a \code{logical} vector of length \code{1} or \code{length(cd)}.
//'
//' @return
//' Argument \code{cd} is returned with \code{logical} vector attribute \code{"latlon"}
//' updated as appropriate.
//'
//' @examples
//' ## Continuing example from `coords()`, named numeric vector representing degrees and minutes
//' \dontshow{
//'    dm <-
//'        c(5130.4659, 4932.7726, 4806.4339, 3853.3696, 0.0000, -3706.7044, -5306.2869, -2514.4093,
//'		   -007.6754, 1823.9137, -12246.7203, -7702.1145, 0.0000, -1217.3178, 7331.0370, -5731.1536)
//'
//'    names(dm) <- 
//'        rep(c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
//'              "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"), 2)
//' }
//'
//' ## Create "coords" object of degrees and minutes
//' coords(dm) <- 2
//'
//' ## Set "latlon" attribute to FALSE, length 1
//' latlon(dm) <- FALSE
//' dm
//'
//' ## Set "latlon" attribute to TRUE and FALSE (each n=8)
//' latlon(dm) <- rep(c(TRUE, FALSE), each = 8)
//' dm
//'
//' ## Reversing latitude and longitude results in an
//' ## invalid longitude value and a warning
//' latlon(dm) <- rep(c(FALSE, TRUE), each = 8)
//' dm
//'
//' rm(dm)
//'
// [[Rcpp::export(name = "`latlon<-`")]]
NumericVector latlon(NumericVector cd, LogicalVector& value)
{
//	cout << "——Rcpp::export——latlon(NumericVector, LogicalVector)\n";
	checkinherits(cd, "coords");
	if (value.size() != cd.size() && value.size() != 1)
		stop("value must be either length 1 or length(cd)");
	else
		cd.attr("latlon") = value;
	validate<NumericVector, Coord>(cd);
	return cd;
}


/// __________________________________________________
/// Validate coords or waypoints vector
//' @title Validate Coords or Waypoints
//' 
//' @name validate
//' 
//' @description
//' \code{validate()} validate objects of class \code{"coords"} or \code{"waypoints"}.
//'
//' @details
//' Individual coordinate values within \code{\link[=coords]{"coords"}} or
//' \code{\link[=waypoints]{"waypoints"}} objects are validated to ensure their being plausible
//' geographic locations.
//'
//' To be valid, the absolute values of coordinates in degrees must not exceed 180, or 90 if degrees
//' of latitude and, similarly, the absolute values of the minutes and seconds components, where
//' given, must not exceed 60 degrees. Otherwise a warning will be issued and the \code{"valid"}
//' attribute in the case of a \code{"coords"} object, or \code{"validlat"} and \code{"validlon"}
//' attributes in the case of a \code{"waypoints"} object will be set to \code{FALSE} for any
//' non-compliant coordinate values.
//'
//' @seealso
//' \code{\link[=coords]{"coords"}} and \code{\link[=waypoints]{"waypoints"}}.
//'
//' @param cd object of class \code{"coords"}.
//'
//' @param df object of class \code{"waypoints"}.
//'
//' @return
//' \code{validate()} returns its argument with \code{logical} vector attribute \code{"valid"},
//' or attributes \code{"validlat"} and \code{"validlon"} updated as appropriate for
//' \code{"coords"} and' \code{"waypoints"} objects respectively.
//'
//' @examples
//' ## Continuing example from `coords()`, named numeric vector representing degrees and minutes
//' \dontshow{
//'    dm <-
//'        c(5130.4659, 4932.7726, 4806.4339, 3853.3696, 0.0000, -3706.7044, -5306.2869, -2514.4093,
//'		   -007.6754, 1823.9137, -12246.7203, -7702.1145, 0.0000, -1217.3178, 7331.0370, -5731.1536)
//'
//'    names(dm) <- 
//'        rep(c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
//'              "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"), 2)
//' }
//'
//' ## Create "coords" object of degrees and minutes
//' coords(dm) <- 2
//'
//' validate(dm)
//'
//' ## Change first value to have more than 60 minutes
//' dm[1] <- 5160.4659
//'
//' validate(dm)
//'
//' ## Examine "valid" attribute of dm
//' attr(dm, "valid")
//'
//' ###
//' ## Continuing example from `waypoints()`, data frame representing waypoint names and latitude
//' ## and longitude values in decimal degrees
//' \dontshow{
//' wp1 <- data.frame(
//'     name = c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
//'              "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"),
//'     lat = c(51.507765, 49.54621, 48.107232, 38.889494, 0, -37.11174, -53.104781, -25.240156),
//'     lon = c(-0.127924, 18.398562, -122.778671, -77.035242, 0, -12.28863, 73.517283, -57.519227)
//' )
//' }
//'
//' ## Create "waypoints" object of decimal degrees
//' waypoints(wp1) <- 1
//'
//' validate(wp1)
//'
//' ## Change penultimate latitude absolute value to have more than 90 degrees
//' wp1$lat[7] <- -93.104781
//'
//' validate(wp1)
//'
//' ## Examine "validlat" attribute of wp1
//' attr(wp1, "validlat")
//'
//' rm(dm, wp1)
//'
// [[Rcpp::export(name = "validate.coords")]]
NumericVector validatecoords(NumericVector cd)
{
//	cout << "——Rcpp::export——validatecoords(NumericVector)\n";
	checkinherits(cd, "coords");
	return validate<NumericVector, Coord>(cd);
}


/// __________________________________________________
/// Format coords vector - S3 method format.coords()
//' @rdname coords
// [[Rcpp::export(name = "format.coords")]]
CharacterVector formatcoords(NumericVector cd)
{
//	cout << "——Rcpp::export——formatcoords(NumericVector)\n";
	checkinherits(cd, "coords");
	if (!check_valid(cd))
		warning("Formatting invalid coords!");
	CoordType ct = get_coordtype(cd);
	return wrap(format_switch(Coord(ct, cd), ct));
}


/// __________________________________________________
/// Print coords vector - S3 method print.coords()
//' @rdname coords
// [[Rcpp::export(name = "print.coords", invisible = true)]]
NumericVector printcoords(NumericVector cd)
{
//	cout << "——Rcpp::export——printcoords(NumericVector) format " << get_fmt_attribute(cd) << endl;
	checkinherits(cd, "coords");
	if (!check_valid(cd))
		warning("Printing invalid coords!");
	Rcout << Coord(get_coordtype(cd), cd);
	return cd;
}


/// __________________________________________________
/// Add "waypoints" to R data.frame object class and validate,
/// or convert format of R waypoints object and return
//' @title Geographic or GPS Waypoint Class
//' 
//' @name waypoints
//' 
//' @description
//' \code{waypoints()} creates a robust representation of a series of geographic or GPS waypoints
//' instantiated as an object of class \code{"waypoints"}.
//' 
//' \code{waypoints()} and replacement form \verb{waypoints()<-} also convert the format of
//' existing objects of class \code{"waypoints"} between (i) decimal degrees, (ii) degrees and
//' minutes, and (iii) degrees, minutes and seconds.
//'
//' @details
//' By default, the names of the waypoints should be included in a "Name" column of data frame
//' argument \code{df}, and the latitude and longitude in the two columns immediately on the right
//' hand side of "Name". An alternative column for waypoint names may be specified by setting an
//' \code{integer} attribute, \code{"namescol"} indicating its position in \code{df}, while setting
//' this attribute to \code{NA} supresses printing of waypoint names. If \code{df} has neither a
//' "Name" column nor a \code{"namescol"} attribute, the \code{"row.names"} attribute is used for
//' waypoint names if present in \code{df}. Similarly, alternative columns for the latitude and
//' longitude may be specified by setting \code{"llcols"} as a length 2 \code{integer} vector
//' attribute indicating their positions in \code{df}.
//'
//' Individual values provided in the \code{numeric} vector latitude and longitude columns of data
//' frame argument \code{df} should have a decimal point after the number of whole degrees in the
//' case of \emph{decimal degrees}, after the number of whole minutes in the case of
//' \emph{degrees and minutes}, and after the number of whole seconds in the case of
//' \emph{degrees, minutes and seconds}.
//'
//' The \code{fmt} argument should be \code{1L} to represent decimal degrees, \code{2L} for degrees
//' and minutes, and \code{3L} for degrees, minutes and seconds and is used to provide both the
//' format of values in data frame argument \code{df} to be converted into a \code{"waypoints"}
//' object and the desired format if a \code{"waypoints"} object is to be converted to a new format.
//' Note that on conversion of a \code{"waypoints"} object, the original data frame argument
//' \code{df} is modified such that the latitude and longitude values are as described in the
//' previous paragraph, and may be inspected using standard R code, see examples.
//'
//' The latitude and longitude values of a newly created \code{"waypoints"} object are checked to
//' ensure they are valid geographic locations as described under
//' \code{\link[=validate]{validate}()}. Likewise, a check is made to ensure that an existing
//' \code{"waypoints"} object to be converted to a new format has already been validated; if not, it
//' is re-validated. 
//'
//' @family coords_waypoints
//' @seealso
//' \code{\link[base:attr]{attr()}}, \code{\link[base:attributes]{attributes}},
//'   \code{\link[base:data.frame]{data.frame()}}, and \code{\link[=validate]{validate}()}.
//'
//' @param df a data frame with each row representing a waypoint, comprising at least two
//'   \code{numeric} columns containing values of latitude and longitude, and optionally a
//'   \code{character} column of waypoint names (see \emph{Details}). 
//'
//' @param fmt,value an \code{integer} of value 1L, 2L or 3L, indicating the current or desired
//'   coordinate format (see \emph{Details}); default 1L.
//'
//' @param wp an object of class \code{"waypoints"} created by function
//' \code{\link[=waypoints]{waypoints}()}.
//'
//' @return
//' An object of classes \code{"waypoints"} and \code{"data.frame"}, comprising the original data
//' frame argument \code{df}, with latitude and longitude values possibly converted as appropriate
//' and additional attributes: –
//' \item{\code{"fmt"}}{the coordinate format.}
//' \item{\code{"namescol"}}{the position of waypoint names, if present within \code{df}.}
//' \item{\code{"llcols"}}{the position of latitude and longitude columns within \code{df}.}
//' \item{\code{"validlat"} and \code{"validlon"}}{\code{logical} vectors indicating whether
//'   individual latitude and longitude values are valid geographic locations.}
//'
//' @examples
//' ## Dataframe representing waypoint names, and latitude and longitude values
//' ## of degrees, minutes and seconds
//' wp1 <- data.frame(
//'     name = c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
//'              "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"),
//'     lat = c(513027.95, 493246.36, 480626.04, 385322.18, 0, -370642.26, -530617.21, -251424.56),
//'     lon = c(-00740.53, 182354.82, -1224643.22, -770206.87, 0, -121719.07, 733102.22, -573109.21)
//' )
//'
//' ## Create "waypoints" object of degrees, minutes and seconds (fmt = 3)
//' waypoints(wp1, 3)
//'
//' ## Convert to degrees and minutes (fmt = 2)
//' waypoints(wp1) <- 2
//' wp1
//'
//' ## Convert to decimal degrees (fmt = 1)
//' waypoints(wp1) <- 1
//' wp1
//'
//' ###
//' ## Dataframe representing unnamed latitude and longitude values in decimal degrees
//' wp2 <- data.frame(
//'     lat = c(51.507765, 49.54621, 48.107232, 38.889494, 0, -37.11174, -53.104781, -25.240156),
//'     lon = c(-0.127924, 18.398562, -122.778671, -77.035242, 0, -12.28863, 73.517283, -57.519227)
//' )
//'
//' ## Create "waypoints" object of decimal degrees (default fmt = 1)
//' waypoints(wp2)
//'
//' ## Add row.names
//' row.names(wp2) <-
//'     c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
//'       "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport")
//' wp2
//'
//' ## Convert to degrees and minutes (fmt = 2)
//' waypoints(wp2) <- 2
//' wp2
//'
//' ## Convert to degrees, minutes and seconds (fmt = 3)
//' waypoints(wp2) <- 3
//' wp2
//'
//' ## Degrees, minutes and seconds as a standard data frame
//' as.data.frame(wp2)
//'
//' ## Convert to decimal degrees, and format as a character vector
//' waypoints(wp2) <- 1
//' (wp2_chr <- format(wp2))
//'
//' ## Output using {base} cat()
//' cat(wp2_chr, fill = 26, labels = paste0("{#", 1:8, "}:"))
//'
//' rm(wp1, wp2)
//'
// [[Rcpp::export]]
DataFrame waypoints(DataFrame df, int fmt = 1)
{
//	cout << "——Rcpp::export——waypoints(DataFrame, int)\n";
	CoordType newtype = get_coordtype(fmt);
	CoordType type;
	if (df.inherits("waypoints")) {
		type = get_coordtype(df);
//		cout << "argument df is already a \"waypoints\" vector of type " << coordtype_to_int(type) + 1 << endl;
		if (newtype == type) {
//			cout << "——fmt out == fmt in!——" << endl;
			if (!check_valid(df))
				stop("Invalid waypoints!");
			return df;
		}
	} else {
		type = newtype;
		df.attr("fmt") = fmt;
		int namescol = 0;
		if (!df.hasAttribute("namescol")) {
			namescol = nameinobj(df, "name");
			if (++namescol)
				df.attr("namescol") = namescol;
		}
		if (!df.hasAttribute("llcols")) {
			const vector<int> llcols { namescol + 1, namescol + 2 };
			df.attr("llcols") = llcols;
		}
	}

	if(!valid_ll(df))
		stop("Invalid llcols attribute!");
	convene<DataFrame, WayPoint>(df, newtype);
	df.attr("class") = CharacterVector{"waypoints", "data.frame"};
	return df;
}


/// __________________________________________________
/// waypoints() as replacement function
//' @rdname waypoints
// [[Rcpp::export(name = "`waypoints<-`")]]
DataFrame waypoints_replace(DataFrame df, int value)
{
//	cout << "——Rcpp::export——`waypoints_replace(DataFrame, int)<-`\n";
	return waypoints(df, value);
}


/// __________________________________________________
/// Validate waypoints vector
//' @rdname validate
// [[Rcpp::export(name = "validate.waypoints")]]
DataFrame validatewaypoints(DataFrame df)
{
//	cout << "——Rcpp::export——validatewaypoints(DataFrame) format " << get_fmt_attribute(df) << endl;
	checkinherits(df, "waypoints");
	if(!valid_ll(df))
		stop("Invalid llcols attribute!");
	return validate<DataFrame, WayPoint>(df);
}


/// __________________________________________________
/// Format waypoints vector - S3 method format.waypoints()
//' @rdname waypoints
// [[Rcpp::export(name = "format.waypoints")]]
CharacterVector formatwaypoints(DataFrame wp)
{
//	cout << "——Rcpp::export——formatwaypoints(NumericVector)\n";
	checkinherits(wp, "waypoints");
	if(!valid_ll(wp))
		stop("Invalid llcols attribute!");
	if (!check_valid(wp))
		warning("Formatting invalid waypoints!");
//	return wrap(WayPoint(get_coordtype(wp), wp).format_switch());
	CoordType ct = get_coordtype(wp);
	return wrap(format_switch(WayPoint(ct, wp), ct));
}


/// __________________________________________________
/// Print waypoints vector - S3 method print.waypoints()
//' @rdname waypoints
// [[Rcpp::export(name = "print.waypoints", invisible = true)]]
DataFrame printwaypoints(DataFrame wp)
{
//	cout << "——Rcpp::export——printwaypoints(DataFrame) format " << get_fmt_attribute(wp) << endl;
	checkinherits(wp, "waypoints");
	if(!valid_ll(wp))
		stop("Invalid llcols attribute!");
	if (!check_valid(wp))
		warning("Printing Invalid waypoints!");
	Rcout << WayPoint(get_coordtype(wp), wp);	
	return wp;
}


/// __________________________________________________
/// Clone coords object from waypoints vector
// [[Rcpp::export]]
NumericVector as_coord(DataFrame df, bool latlon)
{
//	cout << "——Rcpp::export——as_coord(DataFrame)\n;
	checkinherits(df, "waypoints");
	NumericVector nv = df[get_vec_attr<DataFrame, int>(df, "llcols")[latlon ? 0 : 1]];
	nv = clone(nv);
	nv.attr("class") = "coords";
	nv.attr("fmt") = df.attr("fmt");
	nv.attr("latlon") = latlon ? vector<bool>{ TRUE } : vector<bool>{ FALSE };
	nv.attr("valid") = df.attr(latlon ? "validlat" : "validlon");
	return nv;
}


/// __________________________________________________
/// __________________________________________________