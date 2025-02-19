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
template<class T, class V>
inline vector<int> getllcolsattr(const DataFrame);

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
NumericVector printcoord(NumericVector);
NumericVector validatecoord(NumericVector);
CharacterVector formatcoord(NumericVector);
DataFrame waypoints(DataFrame, int);
DataFrame waypoints_replace(DataFrame, int);
DataFrame printwaypoint(DataFrame);
DataFrame validatewaypoint(DataFrame);
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
//	cout << "checkinherits<T>(T& t, const char* classname) t " << typeid(t).name() << " classname " << classname << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	if (!t.inherits(classname)) stop("Argument must be a \"%s\" object", classname);
}


/// __________________________________________________
/// get "llcols" attribute for DataFrame object
inline vector<int> getllcolsattr(const DataFrame df)
{
//	return as<vector<int>>(df.attr("llcols"));
	return get_vec_attr<DataFrame, int>(df, "llcols");
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
//	cout << "@get_coordtype<T>(const T&) " << as<int>(t.attr("fmt")) << endl;
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
					  << setw(7) << setfill('0') << fixed << setprecision(4) << abs(ff.get_decmin(n)) << "'";
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
					  << setw(2) << setfill('0') << abs(ff.get_min(n)) << "'"
					  << setw(5) << fixed << setprecision(2) << abs(ff.get_sec(n)) << "\"";
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
/// Formatted coordinate strings for printing
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
	vector<string> sv; 
	switch (ct)
	{
		case CoordType::decdeg:
			sv = format<CoordType::decdeg>();
			break;

		case CoordType::degmin:
			sv = format<CoordType::degmin>();
			break;

		case CoordType::degminsec:
			sv = format<CoordType::degminsec>();
			break;

		default:
			stop("Coord::print(ostream&) my bad");
	}

	vector<string> names { get_vec_attr<NumericVector, string>(nv, "names") };
	if (names.size())
		transform(sv.begin(), sv.end(), names.begin(), sv.begin(), [](string& lls, const string& name) { return lls + "  " + name; });
	for_each(sv.begin(), sv.end(), [&stream](const string& s) { stream << s << "\n"; });
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
		const char latcol;
		const char loncol;
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
	latcol(getllcolsattr(df)[0]),
	loncol(getllcolsattr(df)[1]),
//	latcol(1), loncol(2),											/////// !!!!!!! Temporary Solution !!!!!!! ///////
	nvlat(df[latcol]), nvlon(df[loncol])
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
/// Formatted character strings for printing
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
	vector<int> spacing { 5, 7, 8, 11, 13, 14, 2, 2, 2 };
	stream << " Latitude" << string(spacing[i], ' ') << "Longitude\n"
		   << string(1, ' ') << string(spacing[i + 3], '_')
		   << string(spacing[i + 6], ' ') << string(spacing[i + 3] + 1, '_') << endl;

	vector<string> sv; 
	switch (ct)
	{
		case CoordType::decdeg:
			sv = format<CoordType::decdeg>();
			break;

		case CoordType::degmin:
			sv = format<CoordType::degmin>();
			break;

		case CoordType::degminsec:
			sv = format<CoordType::degminsec>();
			break;

		default:
			stop("WayPoint::print(ostream&) my bad");
	}

	vector<string> names { as<vector<string>>(df[as<int>(df.attr("namescol"))]) };	/////// !!! revise this !!! ///////
	transform(sv.begin(), sv.end(), names.begin(), sv.begin(), [](string& lls, const string& name) { return lls + "  " + name; });
	for_each(sv.begin(), sv.end(), [&stream](const string& s) { stream << s << "\n"; });
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
	const vector<int> llcols = getllcolsattr(df);
	NumericVector nvlat(df[llcols[0]]);
	NumericVector nvlon(df[llcols[1]]);
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
//' \code{coords()} and replacement form function \verb{coords()<-} also convert the format of
//' existing objects of class \code{"coords"} between (i) decimal degrees, (ii) degrees and minutes,
//' and (iii) degrees, minutes and seconds.
//'
//' @details
//' Individual values provided in the numeric vector argument \code{nv} should have a decimal
//' point after the number of whole degrees in the case of decimal degrees,  after the number of
//' whole minutes in the case of degrees and minutes, and after the number of whole seconds in the
//' case of degrees, minutes and seconds.
//'
//' The \code{fmt} argument is used to provide the format of values in a numeric vector to be
//' converted into a \code{"coords"} object, and the desired  format if a \code{"coords"} object is
//' to be converted to a new format.  \code{fmt} should be 1 for decimal degrees, 2 for degrees and
//' minutes, and 3 for degrees, minutes and seconds.
//'
//' The values of a newly created \code{"coords"} object are validated to ensure their being
//' plausible geographic locations as described under \code{\link[=validate]{validate()}}. Likewise,
//' a check is made to ensure that that an existing \code{"coords"} object to be converted to a new
//' format has already been validated; if not, it is re-validated. 
//'
//' @family coords_waypoints
//'
//' @param nv \code{numeric vector} of coordinate values, optionally named.
//' @param fmt,value \code{integer}, 1L, 2L or 3L, indicating the current or desired coordinate format;
//'   default 1L.
//' @param cd object of class \code{"coords"}, as created by function \code{\link[=coords]{coords()}}.
//'
//' @return
//' An object of class \code{"coords"} comprising a \code{numeric vector} with a
//' \code{boolean vector} attribute \code{"valid"}, indicating whether the individual coordinate
//' values are indeed valid, as described  under \emph{Details}.
//'
//' @examples
//' ## Named numeric vector representing degrees and minutes
//' (num_dm <- c(
//'         5130.4659, 4932.7726, 4806.4339, 3853.3696, 0.0000, -3706.7044, -5306.2869, -5731.1536,
//'         -007.6754, 1823.9137, -12246.7203, -7702.1145,0.0000, -1217.3178, 7331.0370, -2514.4093
//'     ) |>
//'     setNames(rep(
//'         c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
//'           "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"), 2
//'     ))
//' )
//'
//' ## Create "coords" object of degrees and minutes
//' coords(num_dm, 2)
//'
//' ## Convert "coords" object to degrees, minutes and seconds
//' coords(num_dm) <- 3
//' num_dm
//'
//' ## Convert to decimal degrees
//' coords(num_dm) <- 1
//' num_dm
//'
//' ## Convert back to degrees and minutes
//' coords(num_dm) <- 2
//' num_dm
//'
//' rm(num_dm)
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
//' @name latlon-set
//'
//' @description \code{latlon()<-} adds the attribute \code{"latlon"} to objects of class
//' \code{\link[=coords]{"coords"}}, or modifies an existing \code{"latlon"} attribute.
//'
//' @details
//' Attribute \code{"latlon"} is a \code{boolean vector} of length \code{1L} or \code{length(cd)}
//' for which \code{TRUE} values represent latitude and \code{FALSE} values represent longitude.
//' Setting this attribute to any other length will result in an error. A \code{boolean vector} of
//' length \code{1L} signifies that values are all latitude if \code{TRUE}, or longitude if
//' \code{FALSE}.
//'
//' This attribute is used in formatting printed output and also by
//' \code{\link[=validate]{validate()}}. Indeed, the values of \code{cd} are revalidated every time
//' attribute \code{"latlon"} is added or changed.
//'
//' @param cd object of class \code{\link[=coords]{"coords"}}.
//' @param value a \code{boolean vector} of length \code{1L} or of \code{length(cd)}.
//'
//' @return
//' Argument \code{cd} is returned with \code{boolean vector} attribute \code{"latlon"}
//' updated as appropriate.
//'
//' @examples
//' ## Continuing example from `coords()`, named numeric vector representing degrees and minutes
//' \dontshow{
//'     (num_dm <- c(5130.4659, 4932.7726, 4806.4339, 3853.3696, 0.0000, -3706.7044, -5306.2869, -5731.1536,
//'    	             -007.6754, 1823.9137, -12246.7203, -7702.1145, 0.0000, -1217.3178, 7331.0370, -2514.4093) |>
//'         setNames(
//'             rep(
//'                 c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
//'                   "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"), 2
//'             )
//'         )
//'     )
//' }
//'
//' ## Create "coords" object of degrees and minutes
//' coords(num_dm) <- 2
//'
//' ## Set "latlon" attribute to FALSE, length 1
//' latlon(num_dm) <- FALSE
//' num_dm
//'
//' ## Set "latlon" attribute to TRUE and FALSE (each n=8)
//' latlon(num_dm) <- rep(c(TRUE, FALSE), each = 8)
//' num_dm
//'
//' ## Reversing latitude and longitude results in an implausible 
//' ## longitude value and a warning message 
//' latlon(num_dm) <- rep(c(FALSE, TRUE), each = 8)
//' num_dm
//'
//' rm(num_dm)
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
/// Print coords vector - S3 method print.coords()	  /////// "invisible" not working ///////
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
//' The absolute values of coordinates in degrees must not exceed 180, or 90 if degrees of
//' latitude. Likewise the absolute values of the minutes and seconds components, where given,
//' must not exceed 60 degrees, otherwise a warning will be given and the \code{"valid"} attribute
//' in the case of a \code{"coords"} object, or \code{"validlat"} and \code{"validlon"} attributes
//' in the case of a \code{"waypoints"} object set to \code{FALSE} for any non-compliant coordinate
//' values.
//'
//' @param cd object of class \code{"coords"}.
//' @param df object of class \code{"waypoints"}.
//'
//' @return
//' \code{validate()} returns its argument with \code{boolean vector} attribute \code{"valid"},
//' or attributes \code{"validlat"} and \code{"validlon"} updated as appropriate for
//' \code{"coords"} and' \code{"waypoints"} objects respectively.
//'
//' @examples
//' ## Continuing example from `coords()`, named numeric vector representing degrees and minutes
//' \dontshow{
//' num_dm <- c(
//'         5130.4659, 4932.7726, 4806.4339, 3853.3696, 0.0000, -3706.7044, -5306.2869, -5731.1536,
//'         -007.6754, 1823.9137, -12246.7203, -7702.1145,0.0000, -1217.3178, 7331.0370, -2514.4093
//'     ) |>
//'     setNames(rep(
//'         c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
//'           "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"), 2
//'     ))
//' }
//'
//' ## Create "coords" object of degrees and minutes
//' coords(num_dm) <- 2
//'
//' validate(num_dm)
//'
//' ## Change first value to have more than 60 minutes
//' num_dm[1] <- 5160.4659
//'
//' validate(num_dm)
//'
//' ## Continuing example from `waypoints()`, dataframe representing waypoint names and latitude
//' ## and longitude values in decimal degrees
//' \dontshow{
//' wp1 <- data.frame(
//'     name = c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument",
//'              "Null Island", "Tristan da Cunha", "Mawson Peak",
//'              "Silvio Pettirossi International Airport"),
//'     lat = c(51.507765, 49.54621, 48.107232, 38.889494,
//'             0, -37.11174, -53.104781, -57.519227),
//'     lon = c(-0.127924, 18.398562, -122.778671, -77.035242,
//'             0, -12.28863, 73.517283, -25.240156)
//' )
//' }
//'
//' # Create "waypoints" object of decimal degrees
//' waypoints(wp1) <- 1
//'
//' validate(wp1)
//'
//' ## Change third longitude absolute value to have more than 180 degrees
//' wp1$lon[3] <- -182.778671
//'
//' validate(wp1)
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
// [[Rcpp::export(name = "format.coords")]]
CharacterVector formatcoords(NumericVector nv)
{
//	cout << "——Rcpp::export——formatcoords(NumericVector)\n";
	checkinherits(nv, "coords");
	if (!check_valid(nv))
		warning("Formatting invalid coords!");
	Coord c(get_coordtype(nv), nv);

	switch (get_coordtype(nv))
	{
		case CoordType::decdeg:
			return wrap(c.format<CoordType::decdeg>());

		case CoordType::degmin:
			return wrap(c.format<CoordType::degmin>());

		case CoordType::degminsec:
			return wrap(c.format<CoordType::degminsec>());

		default:
			stop("formatcoord(NumericVector) my bad");
	}
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
//' \code{waypoints()} and replacement form function \verb{waypoints()<-} also convert the format
//' of existing objects of class \code{"waypoints"} between (i) decimal degrees, (ii) degrees and
//' minutes, and (iii) degrees, minutes and seconds.
//'
//' @details
//' Individual values provided in the numeric vector latitude and longitude columns of argument
//' \code{df} should have a decimal point after the number of whole degrees in the case of decimal degrees,
//' after the number of whole minutes in the case of degrees and minutes, and after the number of whole
//' seconds in the case of degrees, minutes and seconds.
//'
//' The \code{fmt} argument is used to provide the format of values in a dataframe to be converted into a
//' \code{"waypoints"} object, and the desired  format if a \code{"waypoints"} object is to be converted
//' to a new format.  \code{fmt} should be 1 for decimal degrees, 2 for degrees and minutes, and 3 for degrees, 
//' minutes and seconds.
//'
//' The latitude and longitude values of a newly created \code{"waypoints"} object are validated to ensure
//' their being plausible geographic locations as described under \code{\link[=validate]{validate()}}.
//' Likewise, a check is made to ensure that an existing \code{"coords"} object to be converted to a new
//' format has already been validated; if not, it is re-validated. 
//'
//' @family coords_waypoints
//'
//' @param df a \code{dataframe} containing at least two numeric columns containing corresponding values
//' of latitude and longitude in each row representing a waypoint.
//' @param fmt,value integer, 1, 2 or 3, indicating the current or desired coordinate format; default 1.
//'
//' @return
//' An object of class \code{"waypoints"} comprising a \code{data.frame} with two \code{boolean vector}
//' attributes \code{"validlat"} and \code{"validlon"} indicating whether the individual coordinate values
//' are indeed valid, as described under \emph{Details}.
//'
//' @examples
//' # Dataframe representing waypoint names and latitude and longitude values in decimal degrees
//' (wp1 <- data.frame(
//'     name = c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
//'              "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"),
//'     lat = c(51.507765, 49.54621, 48.107232, 38.889494, 0, -37.11174, -53.104781, -57.519227),
//'     lon = c(-0.127924, 18.398562, -122.778671, -77.035242, 0, -12.28863, 73.517283, -25.240156)
//' ))
//'
//' # Create "waypoints" object of decimal degrees
//' waypoints(wp1)
//'
//' ## Convert "waypoints" object to degrees and minutes
//' waypoints(wp1) <- 2
//' wp1
//'
//' ## Convert "waypoints" object to degrees, minutes and seconds
//' waypoints(wp1) <- 3
//' wp1
//'
//' ## Convert back to decimal degrees
//' waypoints(wp1) <- 1
//' wp1
//'
//' rm(wp1)
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
		if (!df.hasAttribute("llcols")) {
			const vector<int> llcols { 1, 2 };									// !!!!!!!! Temporary Solution !!!!!!
			df.attr("llcols") = llcols;
		}
		if (!df.hasAttribute("namescol")) {
			constexpr int namescol = 0;											// !!!!!!!! Temporary Solution !!!!!!
			df.attr("namescol") = namescol;
		}
	}

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
/// Print waypoints vector - S3 method print.waypoints()	  /////// "invisible" not working ///////
// [[Rcpp::export(name = "print.waypoints", invisible = true)]]
DataFrame printwaypoints(DataFrame df)
{
//	cout << "——Rcpp::export——printwaypoints(DataFrame) format " << get_fmt_attribute(df) << endl;
	checkinherits(df, "waypoints");
	if (!check_valid(df))
		warning("Invalid waypoints!");
	Rcout << WayPoint(get_coordtype(df), df);	
	return df;
}


/// __________________________________________________
/// Validate waypoints vector
//' @rdname validate
// [[Rcpp::export(name = "validate.waypoints")]]
DataFrame validatewaypoints(DataFrame df)
{
//	cout << "——Rcpp::export——validatewaypoints(DataFrame) format " << get_fmt_attribute(df) << endl;
	checkinherits(df, "waypoints");
	return validate<DataFrame, WayPoint>(df);
}


/// __________________________________________________
/// Clone coords object from waypoints vector
// [[Rcpp::export]]
NumericVector as_coord(DataFrame df, bool latlon)
{
//	cout << "——Rcpp::export——as_coord(DataFrame)\n;
	checkinherits(df, "waypoints");
	NumericVector nv = df[getllcolsattr(df)[latlon ? 0 : 1]];
	nv = clone(nv);
	nv.attr("class") = "coords";
	nv.attr("fmt") = df.attr("fmt");
	nv.attr("latlon") = latlon ? vector<bool>{ TRUE } : vector<bool>{ FALSE };
	nv.attr("valid") = df.attr(latlon ? "validlat" : "validlon");
	return nv;
}


/// __________________________________________________
/// __________________________________________________