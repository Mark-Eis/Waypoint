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
void setcolattr(const T&, int, const char*, V&&);			/////// deprecate ///////
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

// waypoint
class WayPoint;

ostream& operator<<(ostream&, const WayPoint&);

// validation
inline bool check_valid(const NumericVector);
// bool check_valid(const DataFrame&);

template<class T>
bool check_valid(T, const char*);

const vector<bool> validate(const NumericVector);

/*
template<class T>
vector<bool> validate(const T&);

template<CoordType type>
vector<bool> validatelet(const DataFrame&);

template<CoordType type>
void setvalidattr(const DataFrame&, WayPoint<type>&);
*/

// conversion
template<CoordType newtype>
inline void convertlet(const Coord&, NumericVector);

/*
template<CoordType type>
void waypointlet(DataFrame&, CoordType newtype);

template<CoordType type, CoordType newtype>
inline void wpconvertlet(DataFrame&, WayPoint<type>&);
*/

// exported
NumericVector coords(NumericVector, const int);
NumericVector coords_replace(NumericVector, int);
NumericVector latlon(NumericVector, LogicalVector&);
NumericVector printcoord(NumericVector);
const vector<bool> validatecoord(NumericVector);
vector<string> formatcoord(NumericVector);
/*
vector<int> get_deg(NumericVector);
vector<double> get_decdeg(NumericVector);
vector<int> get_min(NumericVector);
vector<double> get_decmin(NumericVector);
vector<double> get_sec(NumericVector);
DataFrame waypoints(DataFrame, int);
DataFrame waypoints_replace(DataFrame df, int value);
DataFrame printwaypoint(DataFrame);
const DataFrame validatewaypoint(DataFrame);
*/

/// __________________________________________________
/// __________________________________________________
/// Development and Debugging functions

/// Report object construction and destruction
void _ctrsgn(const type_info& obj, bool destruct = false)
{
	cout << (destruct ? "Destroying " : "Constructing ") << obj.name() << endl;
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
/// set attributes for vector column within object
template<class T, class V>
void setcolattr(const T& t, int col, const char* attrib, V&& val)
{
//	cout << "@setcolattr<T, V>(const T&, int, const char*, V&&) attrib " << attrib << ", col " << col << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	as<NumericVector>(t[col]).attr(attrib) = std::forward<V>(val);
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
	FamousFive() { cout << "§FamousFive() "; _ctrsgn(typeid(*this)); }
	virtual ~FamousFive() = 0;	
	virtual int get_deg(double x) const = 0;
	virtual double get_decdeg(double x) const = 0;
	virtual int get_min(double x) const = 0;
	virtual double get_decmin(double x) const = 0;
	virtual double get_sec(double x) const = 0;
};

FamousFive::~FamousFive() { cout << "§~FamousFive(CoordType) "; _ctrsgn(typeid(*this), true); }	

/// __________________________________________________
/// Derived class for decimal degrees	
struct FF_decdeg : public FamousFive {
	FF_decdeg() { cout << "§FF_decdeg() "; _ctrsgn(typeid(*this)); }	
//	~FF_decdeg() = default;
	~FF_decdeg() { cout << "§~FF_decdeg::FF_decdeg() "; _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x); }
	double get_decdeg(double x) const { return x; }
	int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
	double get_decmin(double x) const { return polish(mod1by60(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
} ff_decdeg;

/// __________________________________________________
/// Derived class for degrees and minutes
struct FF_degmin : public FamousFive {
	FF_degmin() { cout << "§FF_degmin() "; _ctrsgn(typeid(*this)); }	
//	~FF_degmin() = default;
	~FF_degmin() { cout << "§~FF_degmin::FF_degmin() "; _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x / 1e2); }
	double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
	int get_min(double x) const { return int(x) % int(1e2); }
	double get_decmin(double x) const { return polish(mod1e2(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
} ff_degmin;

/// __________________________________________________
/// Derived class for degrees, minutes and seconds
struct FF_degminsec : public FamousFive {
	FF_degminsec() { cout << "§FF_degminsec() "; _ctrsgn(typeid(*this)); }	
//	~FF_degminsec() = default;
	~FF_degminsec() { cout << "§~FF_degminsec::FF_degminsec() "; _ctrsgn(typeid(*this), true); }
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
		const T& t; 
		vector<bool>::const_iterator ll_it;
	public:
		FormatLL(const T& _t) : t(_t), ll_it(t.latlon.begin())
		{
//			cout << "§FormatLL<T, CoordType>::FormatLL(const T&) "; _ctrsgn(typeid(*this));
			static_assert(std::is_same<Coord, T>::value || std::is_same<WayPoint, T>::value, "T must be Coord or WayPoint");
		}
		~FormatLL() = default;
//		~FormatLL() { cout << "§FormatLL<T, CoordType>::~FormatLL() "; _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
//			cout << "@FormatLL<T, CoordType>::operator(string, double) [default for CoordType::degmin and CoordType::degminsec]\n";
			return ostr += t.latlon.size() ? cardpoint(t.ff.get_decmin(n) < 0, t.llgt1 ? *ll_it++ : *ll_it) : cardi_b(t.ff.get_decmin(n) < 0);
		}
};

/// __________________________________________________
/// Specialised functor for decimal degrees
template<class T>
class FormatLL<T, CoordType::decdeg> {
		const T& t; 
		vector<bool>::const_iterator ll_it;
	public:
		FormatLL(const T& _t) : t(_t), ll_it(t.latlon.begin())
		{
//			cout << "§FormatLL<T, CoordType::decdeg>::FormatLL(const T&) "; _ctrsgn(typeid(*this));
			static_assert(std::is_same<Coord, T>::value || std::is_same<WayPoint, T>::value, "T must be Coord or WayPoint");
		}
		~FormatLL() = default;
//		~FormatLL() { cout << "§FormatLL<T, CoordType::decdeg>::~FormatLL() "; _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
//			cout << "@FormatLL<T, CoordType::decdeg>::operator(string, double) t.waypoint " << boolalpha << t.waypoint << endl;
			if (t.latlon.size() && !t.waypoint)
				return ostr += ((t.llgt1 ? *ll_it++ : *ll_it) ? " lat" : " lon");
			else
				return ostr;
		}
};


/// __________________________________________________
/// __________________________________________________
/// Validate Coord functor

class Validator {
		const FamousFive& ff;
		vector<bool>::const_iterator ll_it;
		const int ll_size;
	public:
		Validator(const FamousFive& _ff, const vector<bool>& ll) : ff(_ff), ll_it(ll.begin()), ll_size(ll.size())
		{
			cout << "§Validator::Validator(const FamousFive&, vector<bool>) "; _ctrsgn(typeid(*this));
		}
//		~Validator() = default;
		~Validator() { cout << "§Validator::~Validator() "; _ctrsgn(typeid(*this), true); }
		bool operator()(double n)
		{
			cout << "@Validator() " << " validating: " << setw(9) << setfill(' ') << n << endl;
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
		const vector<bool> latlon;

	public:
		Coordbase(CoordType _ct);
		Coordbase(const Coordbase&) = delete;						// Disallow copying
		Coordbase& operator=(const Coordbase&) = delete;				//  ——— ditto ———
		Coordbase(Coordbase&&) = delete;								// Disallow transfer ownership
		Coordbase& operator=(Coordbase&&) = delete;					// Disallow moving
		virtual ~Coordbase() = 0;

		const FamousFive& get_ff() const;
		virtual void validate(bool warn = true) const = 0;
		virtual const NumericVector get_nv(bool) const = 0;
		virtual const vector<bool>& get_valid(bool) const = 0;
//		virtual const vector<string>& get_names() const = 0;
//		virtual bool all_valid() const = 0;
//		virtual void set_waypoint() const = 0;
//		template<CoordType type>
//		virtual vector<string> format() const = 0;
		virtual void print(ostream&) const = 0;

		template<class T, CoordType type>
		friend class FormatLL;
//		template<class T>
		friend class Validator;
};


Coordbase::Coordbase(CoordType _ct) :
	ct(_ct), ff(*vff[coordtype_to_int(ct)])
{
	cout << "§Coordbase::Coordbase(CoordType) "; _ctrsgn(typeid(*this));
}


Coordbase::~Coordbase()
{
	cout << "§Coordbase::~Coordbase() "; _ctrsgn(typeid(*this), true);
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
		const vector<string> names;
		const bool llgt1 = false;
		bool waypoint = false;

	public:
		Coord(CoordType, const NumericVector);
		Coord(const Coord&) = delete;					// Disallow copying
		Coord& operator=(const Coord&) = delete;			//  ——— ditto ———
		Coord(Coord&&) = delete;							// Disallow transfer ownership
		Coord& operator=(Coord&&) = delete;				// Disallow moving
//		~Coord() = default;
		~Coord() { cout << "§Coord::~Coord() "; _ctrsgn(typeid(*this), true); }

		void validate(bool warn = true) const;
		const NumericVector get_nv(bool) const;
		const vector<bool>& get_valid(bool) const;
		const vector<string>& get_names() const;
//		bool all_valid() const;
		void set_waypoint() const;
		template<CoordType type>
		vector<string> format() const;
		void print(ostream&) const;

		template<class T, CoordType type>
		friend class FormatLL;
//		template<class T>
		friend class Validator;
};


Coord::Coord(CoordType ct, const NumericVector nv) :
	Coordbase(ct), nv(nv),
	latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") },
	names{ get_vec_attr<NumericVector, string>(nv, "names") },
	llgt1(latlon.size() > 1)
{
	cout << "§Coord::Coord(CoordType, const NumericVector) "; _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Validate coords vector
void Coord::validate(bool warn) const
{
//	cout << "@Coord::validate() " << typeid(*this).name() << " latlon " << LogicalVector(wrap(latlon)) << endl;
	vector<bool>& non_const_valid { const_cast<vector<bool>&>(valid) };
	non_const_valid.assign(nv.size(), {false});
	transform(nv.begin(), nv.end(), non_const_valid.begin(), Validator(ff, latlon));
//	if (all_valid())
	if (all_of(valid.begin(), valid.end(), [](bool v) { return v;}))
		non_const_valid.assign({true});
	else
		if (warn)
			warning("Validation failed!");
	const_cast<NumericVector&>(nv).attr("valid") = valid;
}

/* DEPRECATED!
/// __________________________________________________
/// All valid are true
inline bool Coord::all_valid() const
{
//	cout << "@Coord::all_valid()\n";
	return all_of(valid.begin(), valid.end(), [](bool v) { return v;});
} */


/// __________________________________________________
/// Get const reference to nv
inline const NumericVector Coord::get_nv(bool ll = true) const
{
//	cout << "@Coord::get_nv()\n";
	return nv;
}


/// __________________________________________________
/// Get const reference to valid
inline const vector<bool>& Coord::get_valid(bool ll = true) const
{
//	cout << "@Coord::get_valid()\n";
	return valid;
}


/// __________________________________________________
/// Get const reference to names
inline const vector<string>& Coord::get_names() const
{
	return names;
}


/// __________________________________________________
/// Set waypoint flag
inline void Coord::set_waypoint() const
{
//	cout << "@Coord::set_waypoint()\n";
	bool& wpt = const_cast<bool&>(waypoint);
	wpt = true;
}


/// __________________________________________________
/// Formatted coordinate strings for printing
template<CoordType type>
vector<string> Coord::format() const
{
//	cout << "@Coord::format<CoordType>() " << typeid(*this).name() << endl;
	vector<string> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), Format<type>(ff));
	transform(out.begin(), out.end(), nv.begin(), out.begin(), FormatLL<Coord, type>(*this));
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

	if (names.size()) {
		vector<string>::const_iterator nm_it(names.begin());
		for_each(sv.begin(), sv.end(),
			[&stream,& nm_it](const string& s) { stream << s << " " << *nm_it++ << "\n"; });
	} else
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
		const int latcol;
		const int loncol;
		const NumericVector nvlat;
		const NumericVector nvlon;
		const vector<bool> validlat { false };
		const vector<bool> validlon { false };
	public:
		explicit WayPoint(CoordType, const DataFrame);
		WayPoint(const WayPoint&) = delete;					// Disallow copying
		WayPoint& operator=(const WayPoint&) = delete;		//  ——— ditto ———
		WayPoint(WayPoint&&) = delete;						// Disallow transfer ownership
		WayPoint& operator=(WayPoint&&) = delete;			// Disallow moving
//		~WayPoint() = default;
		~WayPoint() { cout << "§WayPoint::~WayPoint() "; _ctrsgn(typeid(*this), true); }

		void validate(bool = true) const;
		const NumericVector get_nv(bool) const;
		const vector<bool>& get_valid(bool) const;
//		bool all_valid() const;
		template<CoordType type>
		vector<string> format() const;
		void print(ostream& stream) const;
};


WayPoint::WayPoint(CoordType ct, const DataFrame df) :
	Coordbase(ct), df(df),
//	latcol(getllcolsattr(df)[0]),
//	loncol(getllcolsattr(df)[1])
	latcol(1), loncol(2),											/////// !!!!!!! Temporary Solution !!!!!!! ///////
	nvlat(df[latcol]), nvlon(df[loncol])
{
	cout << "§WayPoint::WayPoint(CoordType ct, const DataFrame) "; _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Validate WayPoint
void WayPoint::validate(bool warn) const
{
	cout << "@WayPoint::validate(bool)\n";

	vector<bool>& non_const_validlat { const_cast<vector<bool>&>(validlat) };
	non_const_validlat.assign(nvlat.size(), {false});
	transform(nvlat.begin(), nvlat.end(), non_const_validlat.begin(), Validator(ff, vector<bool>{ true }));

	vector<bool>& non_const_validlon { const_cast<vector<bool>&>(validlon) };
	non_const_validlon.assign(nvlat.size(), {false});
	transform(nvlat.begin(), nvlat.end(), non_const_validlon.begin(), Validator(ff, vector<bool>{ false }));

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
			warning("Validation of latitude failed!");
	const_cast<DataFrame&>(df).attr("validlon") = validlon;
}


/// __________________________________________________
/// Get const reference to nv
inline const NumericVector WayPoint::get_nv(bool ll = true) const
{
//	cout << "@WayPoint::get_nv(bool)\n";
	return df[ll ? latcol : loncol];
}


/// __________________________________________________
/// WayPoint validity
const vector<bool>& WayPoint::get_valid(bool latlon) const
{
//	cout << "@WayPoint::get_valid(bool) latlon " << boolalpha << latlon << endl;
	return latlon ? validlat : validlon;
}


/*
/// __________________________________________________
/// WayPoint validity warning
void WayPoint::warn_invalid() const
{
//	cout << "@WayPoint::warn_invalid()\n";
	if (any_of(validlat.begin(), validlat.end(), [](bool v) { return !v;})) {
		warning("Invalid latitude");
	}
	if (any_of(validlon.begin(), validlon.end(), [](bool v) { return !v;})) {
		warning("Invalid longitude");
	}
}
*/

/// __________________________________________________
/// Formatted character strings for printing
template<CoordType type>
vector<string> WayPoint::format() const
{
	cout << "@WayPoint::format()\n";
	vector<string> sv_lat(df.nrows());
	transform(nvlat.begin(), nvlat.end(), sv_lat.begin(), Format<type>(ff));
	vector<string> sv_lon(df.nrows());
	transform(nvlon.begin(), nvlon.end(), sv_lon.begin(), Format<type>(ff));

	// FormatLL goes here…

	vector<string> out(sv_lat.size());
	transform(
		sv_lat.begin(), sv_lat.end(), sv_lon.begin(), out.begin(),
		[](string& latstr, string& lonstr) { return latstr + "  " + lonstr; }
	);
/*	vector<string> names { as<vector<string>>(df[as<int>(df.attr("namescol"))]) };	/////// !!! revise this !!! ///////
	transform(
		out.begin(), out.end(), names.begin(), out.begin(),
		[](string& lls, const string& name) { return lls + "  " + name; }
	); */
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

/*	if (names.size()) {
		vector<string>::const_iterator nm_it(names.begin());
		for_each(sv.begin(), sv.end(),
			[&stream,& nm_it](const string& s) { stream << s << " " << *nm_it++ << "\n"; });
	} else */
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
inline bool check_valid(const NumericVector nv)
{
//	cout << "@check_valid(const NumericVector)" << endl;
	return check_valid(nv, "valid");
}


/// __________________________________________________
/// Check "valid" attribute of object of class T all true
template<class T>
bool check_valid(T t, const char* attrname)
{
//	cout << "@check_valid<T>(T, const char*)" << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	const vector<bool>&& valid = get_vec_attr<T, bool>(t, attrname);
	if (valid.size())
		return all_of(valid.begin(), valid.end(), [](bool v) { return v;});
	else {
		warning("Unvalidated %s! Revalidating…", typeid(t).name());
		validate(t);	
		return check_valid(t, attrname);
	}
}


/// __________________________________________________
/// Validate NumericVector
const vector<bool> validate(const NumericVector nv)
{
//	cout << "@validate(const NumericVector)\n";
	Coord c(get_coordtype(nv), nv);
	c.validate();
	return c.get_valid();	
}


/*
/// __________________________________________________
/// Check "lat_valid" and "lon_valid attributes of DataFrame are all true
bool check_valid(const DataFrame& df)
{
//	cout << "@check_valid(const DataFrame&)\n";

	const vector<bool>&& latvalid = get_vecbool_attr(df, "lat_valid");
	bool boolat = all_of(latvalid.begin(), latvalid.end(), [](bool v) { return v;});
	if (!boolat)
		warning("Invalid latitude!");

	const vector<bool>&& lonvalid = get_vecbool_attr(df, "lon_valid");
	bool boolon = all_of(lonvalid.begin(), lonvalid.end(), [](bool v) { return v;});
	if (!boolon)
		warning("Invalid longitude!");

	return boolat || boolon;
}


/// __________________________________________________
/// Validate DataFrame or NumericVector
template<class T>
vector<bool> validate(const T& t)
{
//	cout << "@validate<>(const T&)\n";
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");

	switch (get_coordtype(t))
	{
		case CoordType::decdeg:
			return validatelet<CoordType::decdeg>(t);

		case CoordType::degmin:
			return validatelet<CoordType::degmin>(t);

		case CoordType::degminsec:
			return validatelet<CoordType::degminsec>(t);

		default:
			stop("validate<>(const T&) my bad");
	}
}


template<CoordType type>
vector<bool> validatelet(const DataFrame& df)
{
//	cout << "@validatelet(const DataFrame&)\n";
	WayPoint<type> wp(df);
	wp.validate(true);
	wp.warn_invalid();
	setvalidattr(df, wp);
	return vector<bool>();							/////// temporary solution ///////
}


template<CoordType type>
void setvalidattr(const DataFrame& df, WayPoint<type>& wp)
{
//	cout << "@setvalidattr(DataFrame&, WayPoint<type>&)\n";
	DataFrame& non_const_df { const_cast<DataFrame&>(df) };
	for (const auto x : { 0, 1 } )
		non_const_df.attr(vector<string>{ "lat_valid", "lon_valid" }[x]) = wp.get_valid(1 - x);
}
*/


/// __________________________________________________
/// __________________________________________________
/// Conversion functions

template<CoordType newtype>
inline void convertlet(const Coord& c, NumericVector nv)
{
//	cout << "@convertlet<CoordType>(const Coord&, NumericVector) newtype " << coordtype_to_int(newtype) + 1 << endl;
	transform(c.get_nv().begin(), c.get_nv().end(), nv.begin(), Convertor<newtype>(c.get_ff()));
}


/*
template<CoordType type>
void waypointlet(DataFrame& df, CoordType newtype)
{
//	cout << "@waypointlet<type>(DataFrame&, const vector<int>&, CoordType) type " << coordtype_to_int(type) + 1
//		 << " newtype " << coordtype_to_int(newtype) + 1 << endl;

	WayPoint<type> wp(df);
	wp.validate();
	wp.warn_invalid();

	if (type != newtype) {
		switch (newtype)
		{
			case CoordType::decdeg:
				wpconvertlet<type, CoordType::decdeg>(df, wp);
				break;

			case CoordType::degmin:
				wpconvertlet<type, CoordType::degmin>(df, wp);
				break;

			case CoordType::degminsec:
				wpconvertlet<type, CoordType::degminsec>(df, wp);
				break;

			default:
				stop("@waypointlet<type>(DataFrame&, CoordType) my bad");
		}
	}

	setvalidattr(df, wp);
	df.attr("class") = CharacterVector{"waypoints", "data.frame"};
}


template<CoordType type, CoordType newtype>
inline void wpconvertlet(DataFrame& df, WayPoint<type>& wp)
{
//	cout << "@wpconvertlet(DataFrame&, WayPoint<type>&)\n";
	const vector<int> llcols = getllcolsattr(df);
	for (const auto x : llcols) {
		convertlet<type, newtype>(wp.get_c(llcols[1] - x), as<NumericVector>(df[x]));
	}
}

*/
/// __________________________________________________
/// __________________________________________________
/// Exported functions

/// __________________________________________________
/// Set R vector object class to coords and return,
/// or convert format of R coords object and return
// [[Rcpp::export]]
NumericVector coords(NumericVector nv, const int fmt = 1)
{
//	cout << "——Rcpp::export——coords()\n";
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
	} else
		type = newtype;

	Coord c(type, nv);
	c.validate();

	if (type != newtype) {
		switch (newtype)
		{
			case CoordType::decdeg:
				convertlet<CoordType::decdeg>(c, nv);
				break;

			case CoordType::degmin:
				convertlet<CoordType::degmin>(c, nv);
				break;

			case CoordType::degminsec:
				convertlet<CoordType::degminsec>(c, nv);
				break;

			default:
				stop("coords(NumericVector nv, const int) my bad");
		}
	}

	nv.attr("class") = "coords";
	nv.attr("fmt") = fmt;
	return nv;
}


/// __________________________________________________
/// coords() as replacement function
// [[Rcpp::export(name = "`coords<-`")]]
NumericVector coords_replace(NumericVector nv, int value)
{
//	cout << "——Rcpp::export——`coords_replace()<-`\n";
	return coords(nv, value);
}


/// __________________________________________________
/// Set latlon attribute on "coords" NumericVector and revalidate
// [[Rcpp::export(name = "`latlon<-`")]]
NumericVector latlon(NumericVector nv, LogicalVector& value)
{
//	cout << "——Rcpp::export——set_latlon()\n";
	checkinherits(nv, "coords");
	if (value.size() != nv.size() && value.size() != 1)
		stop("value must be either length 1 or length(nv)");
	else
		nv.attr("latlon") = value;
	validate(nv);
	return nv;
}


/// __________________________________________________
/// Print coords vector - S3 method print.coords()	  /////// "invisible" not working ///////
// [[Rcpp::export(name = "print.coords", invisible = true)]]
NumericVector printcoord(NumericVector nv)
{
//	cout << "——Rcpp::export——printcoord() format " << get_fmt_attribute(nv) << endl;
	checkinherits(nv, "coords");
	if (!check_valid(nv))
		warning("Printing invalid coords!");
	Rcout << Coord(get_coordtype(nv), nv);
	return nv;
}


/// __________________________________________________
/// Validate coords vector
// [[Rcpp::export(name = "validate.coords")]]
const vector<bool> validatecoord(NumericVector nv)
{
//	cout << "——Rcpp::export——validatecoord()\n";
	checkinherits(nv, "coords");
	return validate(nv);
}


/// __________________________________________________
/// Format coords vector - S3 method format.coords()
// [[Rcpp::export(name = "format.coords")]]
vector<string> formatcoord(NumericVector nv)
{
//	cout << "——Rcpp::export——format()\n";
	checkinherits(nv, "coords");
	if (!check_valid(nv))
		warning("Formatting invalid coords!");
	Coord c(get_coordtype(nv), nv);

	switch (get_coordtype(nv))
	{
		case CoordType::decdeg:
			return c.format<CoordType::decdeg>();

		case CoordType::degmin:
			return c.format<CoordType::degmin>();

		case CoordType::degminsec:
			return c.format<CoordType::degminsec>();

		default:
			stop("formatcoord(NumericVector) my bad");
	}
}


/// __________________________________________________
/// dummy() exported function
// [[Rcpp::export]]
DataFrame dummy(DataFrame df, int fmt)
{
	cout << "——Rcpp::export——`dummy(DataFrame, int)`\n";
	CoordType newtype = get_coordtype(fmt);
	WayPoint wp(newtype, df);
	wp.validate();
//	wp.warn_invalid();
	cout << wp;
	df.attr("fmt") = fmt;
	return df;
}


/*
/// __________________________________________________
/// Add "waypoints" to R data.frame object class and validate,
/// or convert format of R waypoints object and return
// [[Rcpp::export]]
DataFrame waypoints(DataFrame df, int fmt = 1)
{
//	cout << "——Rcpp::export——waypoints()\n";
	CoordType newtype = get_coordtype(fmt);
	CoordType type;
	const bool inheritswaypoints { df.inherits("waypoints") };
	if (inheritswaypoints) {
		type = get_coordtype(df);
//		cout << "argument df is already a \"waypoints\" vector of type " << coordtype_to_int(type) + 1 << endl;
		if (!check_valid(df))
			stop("Invalid waypoints!");
		if (newtype == type) {
		//	cout << "——fmt out == fmt in!——" << endl;
			return df;
		}
	} else {
		type = newtype;
		const vector<int> llcols { 1, 2 };									// !!!!!!!! Temporary Solution !!!!!!
		df.attr("llcols") = llcols;
		for (const auto x : llcols)
			setcolattr(df, x, "latlon", vector<bool>(1, llcols[1] - x));
		constexpr int namescol = 0;											// !!!!!!!! Temporary Solution !!!!!!
		df.attr("namescol") = namescol;
	}

	switch (type)
	{
    		case CoordType::decdeg:
			waypointlet<CoordType::decdeg>(df, newtype);
			break;

		case CoordType::degmin:
			waypointlet<CoordType::degmin>(df, newtype);
			break;

		case CoordType::degminsec:
			waypointlet<CoordType::degminsec>(df, newtype);
			break;

		default:
			stop("waypoints(DataFrame, int) my bad");
	}

	df.attr("fmt") = fmt;
	return df;
}


/// __________________________________________________
/// waypoints() as replacement function
// [[Rcpp::export(name = "`waypoints<-`")]]
DataFrame waypoints_replace(DataFrame df, int value)
{
//	cout << "——Rcpp::export——`waypoints_replace()<-`\n";
	return waypoints(df, value);
} */


/// __________________________________________________
/// Print waypoints vector - S3 method print.waypoints()	  /////// "invisible" not working ///////
// [[Rcpp::export(name = "print.waypoints", invisible = true)]]
DataFrame printwaypoint(DataFrame df)
{
//	cout << "——Rcpp::export——printwaypoint() format " << get_fmt_attribute(df) << endl;
	checkinherits(df, "waypoints");
//	if (!check_valid(df))
//		warning("Invalid waypoints!");

	Rcout << WayPoint(get_coordtype(df), df);	

	return df;
}

/*
/// __________________________________________________
/// Validate waypoints vector
// [[Rcpp::export(name = "validate.waypoints")]]
const DataFrame validatewaypoint(DataFrame df)
{
//	cout << "——Rcpp::export——validatewaypoint(DataFrame) format " << get_fmt_attribute(df) << endl;
	checkinherits(df, "waypoints");
	return validate(df);
}

**************************/


/// __________________________________________________
/// __________________________________________________