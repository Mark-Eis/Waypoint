#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/// __________________________________________________
/// __________________________________________________
/// Class and Function declarations

void _ctrsgn(const type_info&, bool);

inline double mod1by60(double);
inline double mod1e2(double);
inline double round2(double, int);
inline double polish(double);

enum class CoordType : char { decdeg, degmin, degminsec };

inline const CoordType get_coordtype(int);
template<class T>
inline const CoordType get_coordtype(const T&);
inline const int coordtype_to_int(CoordType);

inline string cardpoint(bool, bool);
inline string cardi_b(bool);

//Coord
template<CoordType>
struct FamousFive;
template<>
struct FamousFive<CoordType::decdeg>;
template<>
struct FamousFive<CoordType::degmin>;
template<>
struct FamousFive<CoordType::degminsec>;

template<CoordType type, CoordType newtype>
class Convertor;
template<CoordType type>
class Convertor<type, CoordType::decdeg>;
template<CoordType type>
class Convertor<type, CoordType::degmin>;
template<CoordType type>
class Convertor<type, CoordType::degminsec>;

template<CoordType type>
class Coord;

template<CoordType type>
class Validator;

template<CoordType type>
class Format;

template<CoordType type>
class FormatLL;

template<CoordType type>
ostream& operator<<(ostream&, const Coord<type>&);

// waypoint
template<CoordType type>
class WayPoint;

template<CoordType type>
ostream& operator<<(ostream&, const WayPoint<type>&);

// convenience
template<class T> 
inline vector<bool> get_attr(const T&, const char*);
template<class T>
inline int get_fmt_attribute(const T&);
template<class T>
inline void checkinherits(T&, const char*);
template<class T, class V>
void setcolattr(const T&, int, const char*, V&&);			/////// deprecate ///////
inline vector<int> getllcolsattr(const DataFrame&);

// validation
inline bool validcoord(NumericVector&);

inline bool check_valid(const NumericVector&);
bool check_valid(const DataFrame&);

template<class T>
bool check_valid(const T&, const char*);

template<class T>
vector<bool> validate(const T&);

template<CoordType type>
vector<bool> validatelet(const NumericVector&);

template<CoordType type>
vector<bool> validatelet(const DataFrame&);

template<CoordType type>
void setvalidattr(const NumericVector&, Coord<type>&);

template<CoordType type>
void setvalidattr(const DataFrame&, WayPoint<type>&);

// conversion
template<CoordType type>
void coordlet(NumericVector&, CoordType);

template<CoordType type, CoordType newtype>
inline void convertlet(const Coord<type>&, NumericVector&);

template<CoordType type, CoordType newtype>
inline void convertlet(const Coord<type>&, NumericVector&&);

template<CoordType type>
void waypointlet(DataFrame&, CoordType newtype);

template<CoordType type, CoordType newtype>
inline void wpconvertlet(DataFrame&, WayPoint<type>&);

// exported
NumericVector coords(NumericVector&, int);
NumericVector coords_replace(NumericVector&, int);
NumericVector latlon(NumericVector&, LogicalVector&);
NumericVector printcoord(NumericVector&);
vector<bool> Rvalidatecoord(NumericVector&);
vector<string> formatcoord(NumericVector&);
vector<int> get_deg(NumericVector&);
vector<double> get_decdeg(NumericVector&);
vector<int> get_min(NumericVector&);
vector<double> get_decmin(NumericVector&);
vector<double> get_sec(NumericVector&);
DataFrame waypoints(DataFrame&, int);
DataFrame waypoints_replace(DataFrame& df, int value);
DataFrame printwaypoint(DataFrame&);
const DataFrame validatewaypoint(DataFrame&);

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
/// CoordType enum class

/// __________________________________________________
/// Convert int to CoordType enum
inline const CoordType get_coordtype(int i)
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
/// Famous Five functions
template<CoordType>
struct FamousFive {
};

template<>
struct FamousFive<CoordType::decdeg> {
//	FamousFive() { cout << "§FamousFive<CoordType::decdeg>() "; _ctrsgn(typeid(*this)); 	}	
	int get_deg(double x) const { return int(x); }
	double get_decdeg(double x) const { return x; }
	int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
	double get_decmin(double x) const { return polish(mod1by60(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

template<>
struct FamousFive<CoordType::degmin> {
//	FamousFive() { cout << "§FamousFive<CoordType::degmin>() "; _ctrsgn(typeid(*this)); 	}	
	int get_deg(double x) const { return int(x / 1e2); }
	double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
	int get_min(double x) const { return int(x) % int(1e2); }
	double get_decmin(double x) const { return polish(mod1e2(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

template<>
struct FamousFive<CoordType::degminsec> {
//	FamousFive() { cout << "§FamousFive<CoordType::degminsec>() "; _ctrsgn(typeid(*this)); }	
	int get_deg(double x) const { return int(x / 1e4); }
	double get_decdeg(double x) const { return int(x / 1e4) + (double)int(fmod(x, 1e4) / 1e2) / 60 + mod1e2(x) / 3600; }
	int get_min(double x) const { return (int(x) % int(1e4)) / 1e2; }
	double get_decmin(double x) const { return int(fmod(x, 1e4) / 1e2) + mod1e2(x) / 60; }
	double get_sec(double x) const { return mod1e2(x); }
};


/// __________________________________________________
/// __________________________________________________
/// Templated Coord conversion functors
template<CoordType type, CoordType newtype>
class Convertor {
};

/// __________________________________________________
/// Specialised Coord conversion functor for decimal degrees
template<CoordType type>
class Convertor<type, CoordType::decdeg> {
	protected:
		const Coord<type>& c; 
	public:
		Convertor(const Coord<type>& _c) : c(_c)
		{
//			cout << "§Convertor<type, CoordType::decdeg>::Convertor(const Coord<type>&) "; _ctrsgn(typeid(*this));
		}
		~Convertor() = default;
//		~Convertor() { cout << "§Convertor<type, CoordType::decdeg>::~Convertor() "; _ctrsgn(typeid(*this), true); }
		double operator()(double n) { return c.ff.get_decdeg(n); }
};

/// __________________________________________________
/// Specialised Coord conversion functor for degrees and minutes
template<CoordType type>
class Convertor<type, CoordType::degmin> {
	protected:
		const Coord<type>& c; 
	public:
		Convertor(const Coord<type>& _c) : c(_c)
		{
//			cout << "§Convertor<type, CoordType::degmin>::Convertor(const Coord<type>&) "; _ctrsgn(typeid(*this));
		}
		~Convertor() = default;
//		~Convertor() { cout << "§Convertor<type, CoordType::degmin>::~Convertor() "; _ctrsgn(typeid(*this), true); }
		double operator()(double n) { return c.ff.get_deg(n) * 1e2 + c.ff.get_decmin(n); }
};

/// __________________________________________________
/// Specialised Coord conversion functor for degrees, minutes and seconds
template<CoordType type>
class Convertor<type, CoordType::degminsec> {
	protected:
		const Coord<type>& c; 
	public:
		Convertor(const Coord<type>& _c) : c(_c)
		{
//			cout << "§Convertor<type, CoordType::degminsec>::Convertor(const Coord<type>&) "; _ctrsgn(typeid(*this));
		}
		~Convertor() = default;
//		~Convertor() { cout << "§Convertor<type, CoordType::degminsec>::~Convertor() "; _ctrsgn(typeid(*this), true); }
		double operator()(double n) { return c.ff.get_deg(n) * 1e4 + c.ff.get_min(n) * 1e2 + c.ff.get_sec(n); }
};


/// __________________________________________________
/// __________________________________________________
/// Coordinate class
template<CoordType type>
class Coord {
	protected:
		vector<double> nv;
		FamousFive<type> ff;
		const vector<bool> valid { false };
		const vector<bool> latlon;
		const vector<string> names;
		const bool llgt1 = false;
		bool all_valid() const;
		bool waypoint = false;

	public:
		Coord(const NumericVector&);
		Coord(const Coord&) = delete;					// Disallow copying
		Coord& operator=(const Coord&) = delete;			//  ——— ditto ———
		Coord(Coord&&) = delete;							// Disallow transfer ownership
		Coord& operator=(Coord&&) = delete;				// Disallow moving
		~Coord() = default;
//		~Coord<type>() { cout << "§Coord<type>::~Coord() "; _ctrsgn(typeid(*this), true); }

		void validate(bool = true) const;
		const vector<double>& get_nv() const;
		const vector<bool>& get_valid() const;
		const vector<string>& get_names() const;
		void warn_invalid() const;
		void set_waypoint() const;
		vector<string> format() const;
		void print(ostream&) const;

		friend class Coord<CoordType::decdeg>;
		friend class Coord<CoordType::degmin>;
		friend class Coord<CoordType::degminsec>;
		friend class Convertor<type, CoordType::decdeg>;
		friend class Convertor<type, CoordType::degmin>;
		friend class Convertor<type, CoordType::degminsec>;
		friend class Format<type>;
		friend class FormatLL<type>;
		friend class Validator<type>;
};


template<CoordType type>
Coord<type>::Coord(const NumericVector& nv) :
	nv(as<vector<double>>(nv)),
	latlon{ nv.hasAttribute("latlon") ? as<vector<bool>>(nv.attr("latlon")) : vector<bool>() },
	names{ nv.hasAttribute("names") ? as<vector<string>>(nv.attr("names")) : vector<string>() },
	llgt1(latlon.size() > 1)
{
//	cout << "§Coord<type>::Coord(const NumericVector&) "; _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// __________________________________________________
/// Validate Coord functor

template<CoordType type>
class Validator {
		const Coord<type>& c; 
		vector<bool>::const_iterator ll_it;
	public:
		Validator(const Coord<type>& _c) : c(_c), ll_it(c.latlon.begin())
		{
		//	cout << "§Validator<type>::Validator(const Coord<type>&) "; _ctrsgn(typeid(*this));
		}
		~Validator() = default;
//		~Validator() { cout << "§Validator<type>::~Validator(const Coord<type>&) "; _ctrsgn(typeid(*this), true); }
		bool operator()(double n)
		{
		//	cout << "@Validator() " << " validating: " << setw(9) << setfill(' ') << n << endl;
			return !((abs(c.ff.get_decdeg(n)) > (c.latlon.size() && (c.llgt1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(abs(c.ff.get_decmin(n)) >= 60) ||
				(abs(c.ff.get_sec(n)) >= 60));
		}
};


/// __________________________________________________
/// Validate coords vector
template<CoordType type>
void Coord<type>::validate(bool warn) const
{
//	cout << "@Coord<type>::validate() " << typeid(*this).name() << " latlon " << LogicalVector(wrap(latlon)) << endl;
	vector<bool>& non_const_valid { const_cast<vector<bool>&>(valid) };
	non_const_valid.assign(nv.size(), {false});
	transform(nv.begin(), nv.end(), non_const_valid.begin(), Validator(*this));
	if (all_valid())
		non_const_valid.assign({true});
	else
		if (warn)
			warning("Validation failed!");
}


/// __________________________________________________
/// All valid are true
template<CoordType type>
bool Coord<type>::all_valid() const
{
//	cout << "@Coord<type>::all_valid()\n";
	return all_of(valid.begin(), valid.end(), [](bool v) { return v;});
}


/// __________________________________________________
/// Get const reference to nv
template<CoordType type>
inline const vector<double>& Coord<type>::get_nv() const
{
//	cout << "@Coord<type>::get_nv()\n";
	return nv;
}


/// __________________________________________________
/// Get const reference to valid
template<CoordType type>
inline const vector<bool>& Coord<type>::get_valid() const
{
	return valid;
}


/// __________________________________________________
/// Get const reference to names
template<CoordType type>
inline const vector<string>& Coord<type>::get_names() const
{
	return names;
}


/// __________________________________________________
/// Warn if any valid are false
template<CoordType type>
void Coord<type>::warn_invalid() const
{
//	cout << "@Coord<type>::warn_invalid()\n";
	if (!all_valid())
		warning("Validation failed!");
}


/// __________________________________________________
/// Set waypoint flag
template<CoordType type>
inline void Coord<type>::set_waypoint() const
{
//	cout << "@Coord<type>::set_waypoint()\n";
	bool& wpt = const_cast<bool&>(waypoint);
	wpt = true;
}


/// __________________________________________________
/// __________________________________________________
/// Templated formatting functors
template<CoordType type>
class Format {
	protected:
		const Coord<type>& c;
		ostringstream outstrstr;
	public:
		Format(const Coord<type>& _c) : c(_c)
		{
		//	cout << "§Format<type>::Format(const Coord<type>&) "; _ctrsgn(typeid(*this));
		}
		~Format() = default;
//		~Format() { cout << "§Format<type>::~Format() "; _ctrsgn(typeid(*this), true); }
		string operator()(double n);
};

/// __________________________________________________
/// Default operator(), for decimal degrees
template<CoordType type>
inline string Format<type>::operator()(double n)
{
//	cout << "@Format<type>::operator() [default for CoordType::decdeg]\n";
	outstrstr.str("");
	outstrstr << setw(11) << setfill(' ')  << fixed << setprecision(6) << c.ff.get_decdeg(n) << "\u00B0";
	return outstrstr.str();
}

/// __________________________________________________
/// Specialised operator() for degrees and minutes
template<>
inline string Format<CoordType::degmin>::operator()(double n)
{
//	cout << "@Format<CoordType::degmin>::operator()\n";
	outstrstr.str("");
	outstrstr << setw(3) << setfill(' ') << abs(c.ff.get_deg(n)) << "\u00B0"
					  << setw(7) << setfill('0') << fixed << setprecision(4) << abs(c.ff.get_decmin(n)) << "'";
	return outstrstr.str();
}

/// __________________________________________________
/// Specialised operator() for degrees, minutes and seconds
template<>
inline string Format<CoordType::degminsec>::operator()(double n)
{
//	cout << "@Format<CoordType::degminsec>::operator()\n";
	outstrstr.str("");
	outstrstr << setw(3) << setfill(' ') << abs(c.ff.get_deg(n)) << "\u00B0"
					  << setw(2) << setfill('0') << abs(c.ff.get_min(n)) << "'"
					  << setw(5) << fixed << setprecision(2) << abs(c.ff.get_sec(n)) << "\"";
	return outstrstr.str();
}


/// __________________________________________________
/// __________________________________________________
/// Formatting functors for Latitude and Longitude
template<CoordType type>
class FormatLL {
		const Coord<type>& c; 
		vector<bool>::const_iterator ll_it;
	public:
		FormatLL(const Coord<type>& _c) : c(_c), ll_it(c.latlon.begin())
		{
		//	cout << "§FormatLL<type>::FormatLL(const Coord<type>&) "; _ctrsgn(typeid(*this));
		}
		~FormatLL() = default;
//		~FormatLL() { cout << "§FormatLL<type>::~FormatLL() "; _ctrsgn(typeid(*this), true); }
		string operator()(string, double);
};

/// __________________________________________________
/// Default operator(), for degrees, minutes (and seconds)
template<CoordType type>
inline string FormatLL<type>::operator()(string ostr, double n)
{
//	cout << "@FormatLL<type>::operator(string, double) [default for CoordType::degmin and CoordType::degminsec]\n";
	return ostr += c.latlon.size() ? cardpoint(c.ff.get_decmin(n) < 0, c.llgt1 ? *ll_it++ : *ll_it) : cardi_b(c.ff.get_decmin(n) < 0);
}

/// __________________________________________________
/// Specialised operator(), for decimal degrees
template<>
inline string FormatLL<CoordType::decdeg>::operator()(string ostr, double n)
{
//	cout << "@FormatLL<CoordType::decdeg>::operator(string, double) c.waypoint " << boolalpha << c.waypoint << endl;
	if (c.latlon.size() && !c.waypoint)
		return ostr += ((c.llgt1 ? *ll_it++ : *ll_it) ? " lat" : " lon");
	else
		return ostr;
}


/// __________________________________________________
/// Formatted coordinate strings for printing
template<CoordType type>
vector<string> Coord<type>::format() const
{
//	cout << "@Coord<type>::format() " << typeid(*this).name() << endl;
	vector<string> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), Format<type>(*this));
	transform(out.begin(), out.end(), nv.begin(), out.begin(), FormatLL<type>(*this));
	return out;
}


/// __________________________________________________
/// Print coords vector
template<CoordType type>
void Coord<type>::print(ostream& stream) const
{
//	cout << "@Coord<type>::print() " << typeid(*this).name() << endl;
	vector<string> sv(format()); 
	if (names.size()) {
		vector<string>::const_iterator nm_it(names.begin());
		for_each(sv.begin(), sv.end(),
			[&stream,& nm_it](const string& s) { stream << s << " " << *nm_it++ << "\n"; });
	} else
		for_each(sv.begin(), sv.end(), [&stream](const string& s) { stream << s << "\n"; });
}


/// __________________________________________________
/// Output Coord derived object to ostream
template<CoordType type>
ostream& operator<<(ostream& stream, const Coord<type>& c)
{
//	cout << "@operator<<(ostream&, const Coord<type>&)\n";
	c.print(stream);
	return stream;
}


/// __________________________________________________
/// __________________________________________________
/// Waypoint class

template<CoordType type>
class WayPoint {
	protected:
		const DataFrame& df;
		const Coord<type> c_lat;
		const Coord<type> c_lon;
		const vector<bool>& validlat;
		const vector<bool>& validlon;
	public:
		explicit WayPoint(const DataFrame&);
		~WayPoint() = default;
//		~WayPoint() { cout << "§WayPoint::~WayPoint() "; _ctrsgn(typeid(*this), true); }

		const Coord<type>& get_c(bool) const;
		void validate(bool = true) const;
		const vector<bool>& get_valid(bool) const;
		void warn_invalid() const;
		void print(ostream& stream) const;
		vector<string> format() const;
};


template<CoordType type>
WayPoint<type>::WayPoint(const DataFrame& _df) :
	df(_df), c_lat(as<NumericVector>(df[getllcolsattr(df)[0]])), c_lon(as<NumericVector>(df[getllcolsattr(df)[1]])),
	validlat(c_lat.get_valid()), validlon(c_lon.get_valid())
{
//	cout << "§WayPoint<type>::WayPoint(const DataFrame) "; _ctrsgn(typeid(*this));
	c_lat.set_waypoint();
	c_lon.set_waypoint();
}


/// __________________________________________________
/// Get const reference to c_lat or c_lon
template<CoordType type>
inline const Coord<type>& WayPoint<type>::get_c(bool latlon) const
{
//	cout << "@WayPoint<type>::get_c(bool) const\n";
	return latlon ? c_lat : c_lon;
}


/// __________________________________________________
/// Validate WayPoint
template<CoordType type>
void WayPoint<type>::validate(bool warn) const
{
//	cout << "@WayPoint<type>::validate(bool)\n";
	c_lat.validate(warn);
	c_lon.validate(warn);
}


/// __________________________________________________
/// WayPoint validity
template<CoordType type>
const vector<bool>& WayPoint<type>::get_valid(bool latlon) const
{
//	cout << "@WayPoint<type>::get_valid(bool) latlon " << boolalpha << latlon << endl;
	return latlon ? validlat : validlon;
}


/// __________________________________________________
/// WayPoint validity warning
template<CoordType type>
void WayPoint<type>::warn_invalid() const
{
//	cout << "@WayPoint<type>::warn_invalid()\n";
	if (any_of(validlat.begin(), validlat.end(), [](bool v) { return !v;})) {
		warning("Invalid latitude");
	}
	if (any_of(validlon.begin(), validlon.end(), [](bool v) { return !v;})) {
		warning("Invalid longitude");
	}
}


/// __________________________________________________
/// Formatted character strings for printing
template<CoordType type>
vector<string> WayPoint<type>::format() const
{
//	cout << "@WayPoint<type>::format()\n";
	vector<string> sv_lat{ c_lat.format() };
	vector<string> sv_lon{ c_lon.format() };
	vector<string> out(sv_lat.size());
	transform(
		sv_lat.begin(), sv_lat.end(), sv_lon.begin(), out.begin(),
		[](string& latstr, string& lonstr) { return latstr + "  " + lonstr; }
	);
	vector<string> names { as<vector<string>>(df[as<int>(df.attr("namescol"))]) };	/////// !!! revise this !!! ///////
	transform(
		out.begin(), out.end(), names.begin(), out.begin(),
		[](string& lls, const string& name) { return lls + "  " + name; }
	);
	return out;
}


/// __________________________________________________
/// Print WayPoint
template<CoordType type>
void WayPoint<type>::print(ostream& stream) const
{
//	cout << "@WayPoint<type>::print() " << typeid(*this).name() << endl;
	const int i { coordtype_to_int(type) };
	vector<int> spacing { 5, 7, 8, 11, 13, 14, 2, 2, 2 };
	stream << " Latitude" << string(spacing[i], ' ') << "Longitude\n"
		   << string(1, ' ') << string(spacing[i + 3], '_')
		   << string(spacing[i + 6], ' ') << string(spacing[i + 3] + 1, '_') << endl;
	vector<string> sv(format());
	for_each(sv.begin(), sv.end(), [&stream](const string& s) { stream << s << "\n"; });
}


/// __________________________________________________
/// Output WayPoint to ostream
template<CoordType type>
ostream& operator<<(ostream& stream, const WayPoint<type>& wp)
{
//	cout << "@operator<<(ostream&, const WayPoint<type>&)\n";
	wp.print(stream);
	return stream;
}


/// __________________________________________________
/// __________________________________________________
/// Convenience functions


/// __________________________________________________
/// Return named attribute as vector<bool> or empty vector<bool>
template<class T> 
inline vector<bool> get_attr(const T& t, const char* attrname)
{
//	cout << "@get_attr<>(const T&, const char*) attr \"" << attrname << "\" " << boolalpha << t.hasAttribute(attrname) << endl;
	return (t.hasAttribute(attrname) ? as<vector<bool>>(t.attr(attrname)) : vector<bool>());
}


/// __________________________________________________
/// Return "fmt" attribute as int
template<class T>
inline int get_fmt_attribute(const T& t)
{
//	cout << "@get_fmt_attribute<T>(const T&) " << as<int>(t.attr("fmt")) << endl;
	return as<int>(t.attr("fmt"));
}


/// __________________________________________________
/// Does object inherit given class?
template<class T>
inline void checkinherits(T& t, const char* classname)
{
//	cout << "checkinherits(T& t, const char* classname) t " << typeid(t).name() << " classname " << classname << endl;
	if (!t.inherits(classname)) stop("Argument must be a \"%s\" object", classname);
}


/// __________________________________________________
/// set attributes for vector column within object
template<class T, class V>
void setcolattr(const T& t, int col, const char* attrib, V&& val)
{
//	cout << "@setcolattr(const T&, int, const char*, V&&) attrib " << attrib << ", col " << col << endl;
	as<NumericVector>(t[col]).attr(attrib) = std::forward<V>(val);
}


/// __________________________________________________
/// get "llcols" attribute for DataFrame object
inline vector<int> getllcolsattr(const DataFrame& df)
{
	return as<vector<int>>(df.attr("llcols"));
}


/// __________________________________________________
/// __________________________________________________
/// Validation functions

/// __________________________________________________
/// Has R coords object been validated?
inline bool validcoord(NumericVector& nv)
{
//	cout << "@validcoord(NumericVector&)\n";
	LogicalVector lv { as<LogicalVector>(nv.attr("valid")) };
	return 1 == lv.size() && lv[0];
}


/// __________________________________________________
/// Check "valid" attribute of NumericVector all true
inline bool check_valid(const NumericVector& nv)
{
//	cout << "@check_valid(const NumericVector&)" << endl;
	return check_valid(nv, "valid");
}


/// __________________________________________________
/// Check "valid" attribute of object of class T all true
template<class T>
bool check_valid(const T& t, const char* attrname)
{
//	cout << "@check_valid(const T&, const char*)" << endl;
	const vector<bool>&& valid = get_attr(t, attrname);
	if (valid.size())
		return all_of(valid.begin(), valid.end(), [](bool v) { return v;});
	else {
		warning("Unvalidated %s! Revalidating…", typeid(t).name());
		validate(t);	
		return check_valid(t, attrname);
	}
}


/// __________________________________________________
/// Check "lat_valid" and "lon_valid attributes of DataFrame are all true
bool check_valid(const DataFrame& df)
{
//	cout << "@check_valid(const DataFrame&)\n";

	const vector<bool>&& latvalid = get_attr(df, "lat_valid");
	bool boolat = all_of(latvalid.begin(), latvalid.end(), [](bool v) { return v;});
	if (!boolat)
		warning("Invalid latitude!");

	const vector<bool>&& lonvalid = get_attr(df, "lon_valid");
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

	// Code to constrain template to DataFrame or NumericVector

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
vector<bool> validatelet(const NumericVector& nv)
{
//	cout << "@validatelet(const NumericVector&)\n";
	Coord<type> c(nv);
	c.validate();
	setvalidattr(nv, c);
	return c.get_valid();	
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
void setvalidattr(const NumericVector& nv, Coord<type>& c)
{
//	cout << "@setvalidattr(const NumericVector&, Coord<type>&)\n";
	const_cast<NumericVector&>(nv).attr("valid") = c.get_valid();
}


template<CoordType type>
void setvalidattr(const DataFrame& df, WayPoint<type>& wp)
{
//	cout << "@setvalidattr(DataFrame&, WayPoint<type>&)\n";
	DataFrame& non_const_df { const_cast<DataFrame&>(df) };
	for (const auto x : { 0, 1 } )
		non_const_df.attr(vector<string>{ "lat_valid", "lon_valid" }[x]) = wp.get_valid(1 - x);
}


/// __________________________________________________
/// __________________________________________________
/// Conversion functions

template<CoordType type>
void coordlet(NumericVector& nv, CoordType newtype)
{
//	cout << "@coordlet(NumericVector&, CoordType) type " << coordtype_to_int(type) + 1 << " newtype " << coordtype_to_int(newtype) + 1 << endl;
	Coord<type> c(nv);
	c.validate();
	if (type != newtype) {
		switch (newtype)
		{
			case CoordType::decdeg:
				convertlet<type, CoordType::decdeg>(c, nv);
				break;

			case CoordType::degmin:
				convertlet<type, CoordType::degmin>(c, nv);
				break;

			case CoordType::degminsec:
				convertlet<type, CoordType::degminsec>(c, nv);
				break;

			default:
				stop("coordlet(NumericVector&, CoordType) my bad");
		}
	}

	setvalidattr(nv, c);
	nv.attr("class") = "coords";
}


template<CoordType type, CoordType newtype>
inline void convertlet(const Coord<type>& c, NumericVector& nv)
{
//	cout << "@convertlet(const Coord<type>&, NumericVector&) type " << coordtype_to_int(type) + 1 << " newtype " << coordtype_to_int(newtype) + 1 << endl;
	transform(c.get_nv().begin(), c.get_nv().end(), nv.begin(), Convertor<type, newtype>(c));
}


template<CoordType type, CoordType newtype>
inline void convertlet(const Coord<type>& c, NumericVector&& nv)
{
//	cout << "@convertlet(const Coord<type>&, NumericVector&&) type " << coordtype_to_int(type) + 1 << " newtype " << coordtype_to_int(newtype) + 1 << endl;
	transform(c.get_nv().begin(), c.get_nv().end(), nv.begin(), Convertor<type, newtype>(c));
}


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


/// __________________________________________________
/// __________________________________________________
/// Exported functions


/// __________________________________________________
/// Set R vector object class to coords and return,
/// or convert format of R coords object and return
// [[Rcpp::export]]
NumericVector coords(NumericVector& nv, int fmt = 1)
{
//	cout << "——Rcpp::export——coords()\n";
	CoordType newtype = get_coordtype(fmt);
	const bool inheritscoords { nv.inherits("coords") };
	CoordType type;
	if (inheritscoords) {
		type = get_coordtype(nv);
//		cout <<  "coords() argument nv is already a \"coords\" vector of type "
//			 << coordtype_to_int(type) + 1 << endl;
		if (!check_valid(nv))
			stop("Invalid coords!");
		if (newtype == type) {
//			cout << "——fmt out == fmt in!——" << endl;
			return nv;
		}
	} else
		type = newtype;

	switch (type)
	{
		case CoordType::decdeg:
			coordlet<CoordType::decdeg>(nv, newtype);
			break;

		case CoordType::degmin:
			coordlet<CoordType::degmin>(nv, newtype);
			break;

		case CoordType::degminsec:
			coordlet<CoordType::degminsec>(nv, newtype);
			break;

		default:
			stop("coords(NumericVector& nv, int) my bad");
	}
	nv.attr("fmt") = fmt;
	return nv;
}


/// __________________________________________________
/// coords() as replacement function
// [[Rcpp::export(name = "`coords<-`")]]
NumericVector coords_replace(NumericVector& nv, int value)
{
//	cout << "——Rcpp::export——`coords_replace()<-`\n";
	return coords(nv, value);
}


/// __________________________________________________
/// Set latlon attribute on "coords" NumericVector and revalidate
// [[Rcpp::export(name = "`latlon<-`")]]
NumericVector latlon(NumericVector& nv, LogicalVector& value)
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
NumericVector printcoord(NumericVector& nv)
{
//	cout << "——Rcpp::export——printcoord() format " << get_fmt_attribute(nv) << endl;
	checkinherits(nv, "coords");
	if (!check_valid(nv))
		warning("Printing invalid coords!");

	switch (get_coordtype(nv))
	{
		case CoordType::decdeg:
			Rcout << Coord<CoordType::decdeg>(nv);
			break;

		case CoordType::degmin:
			Rcout << Coord<CoordType::degmin>(nv);
			break;

		case CoordType::degminsec:
			Rcout << Coord<CoordType::degminsec>(nv);
			break;

		default:
			stop("printcoord(NumericVector&) my bad");
	}
	
	return nv;
}


/// __________________________________________________
/// Validate coords vector
// [[Rcpp::export(name = "validate.coords")]]
vector<bool> validatecoord(NumericVector& nv)
{
//	cout << "——Rcpp::export——validatecoord()\n";
	checkinherits(nv, "coords");
	return validate(nv);
}


/// __________________________________________________
/// Format coords vector - S3 method format.coords()
// [[Rcpp::export(name = "format.coords")]]
vector<string> formatcoord(NumericVector& nv)
{
//	cout << "——Rcpp::export——format()\n";
	checkinherits(nv, "coords");
	if (!check_valid(nv))
		warning("Formatting invalid coords!");

	switch (get_coordtype(nv))
	{
		case CoordType::decdeg:
			return Coord<CoordType::decdeg>(nv).format();

		case CoordType::degmin:
			return Coord<CoordType::degmin>(nv).format();

		case CoordType::degminsec:
			return Coord<CoordType::degminsec>(nv).format();

		default:
			stop("formatcoord(NumericVector&) my bad");
	}
}

/*
/// __________________________________________________
/// Return degrees as integer
// [[Rcpp::export]]
vector<int> get_deg(NumericVector& nv)
{
//	cout << "——Rcpp::export——get_deg()\n";
	checkinherits(nv, "coords");
	unique_ptr<const CoordBase> c{newconstCoordBase(nv, get_coordtype(nv))};
	vector<int> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), [&c](double n) { return c->get_deg(n); });
	return out;
}


/// __________________________________________________
/// Return decimal degrees as double
// [[Rcpp::export]]
vector<double> get_decdeg(NumericVector& nv)
{
//	cout << "——Rcpp::export——get_decdeg()\n";
	checkinherits(nv, "coords");
	unique_ptr<const CoordBase> c{newconstCoordBase(nv, get_coordtype(nv))};
	vector<double> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), [&c](double n) { return c->get_decdeg(n); });
	return out;
}


/// __________________________________________________
/// Return minutes as integer
// [[Rcpp::export]]
vector<int> get_min(NumericVector& nv)
{
//	cout << "——Rcpp::export——get_min()\n";
	checkinherits(nv, "coords");
	unique_ptr<const CoordBase> c{newconstCoordBase(nv, get_coordtype(nv))};
	vector<int> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), [&c](double n) { return c->get_min(n); });
	return out;
}


/// __________________________________________________
/// Return decimal minutes as double
// [[Rcpp::export]]
vector<double> get_decmin(NumericVector& nv)
{
//	cout << "——Rcpp::export——get_decmin()\n";
	checkinherits(nv, "coords");
	unique_ptr<const CoordBase> c{newconstCoordBase(nv, get_coordtype(nv))};
	vector<double> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), [&c](double n) { return c->get_decmin(n); });
	return out;
}


/// __________________________________________________
/// Return decimal seconds as double
// [[Rcpp::export]]
vector<double> get_sec(NumericVector& nv)
{
//	cout << "——Rcpp::export——get_sec()\n";
	checkinherits(nv, "coords");
	unique_ptr<const CoordBase> c{newconstCoordBase(nv, get_coordtype(nv))};
	vector<double> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), [&c](double n) { return c->get_sec(n); });
	return out;
}
*/

/// __________________________________________________
/// Add "waypoints" to R data.frame object class and validate,
/// or convert format of R waypoints object and return
// [[Rcpp::export]]
DataFrame waypoints(DataFrame& df, int fmt = 1)
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
			stop("waypoints(DataFrame&, int) my bad");
	}

	df.attr("fmt") = fmt;
	return df;
}


/// __________________________________________________
/// waypoints() as replacement function
// [[Rcpp::export(name = "`waypoints<-`")]]
DataFrame waypoints_replace(DataFrame& df, int value)
{
//	cout << "——Rcpp::export——`waypoints_replace()<-`\n";
	return waypoints(df, value);
}


/// __________________________________________________
/// Print waypoints vector - S3 method print.waypoints()	  /////// "invisible" not working ///////
// [[Rcpp::export(name = "print.waypoints", invisible = true)]]
DataFrame printwaypoint(DataFrame& df)
{
//	cout << "——Rcpp::export——printwaypoint() format " << get_fmt_attribute(df) << endl;
	checkinherits(df, "waypoints");
	if (!check_valid(df))
		warning("Invalid waypoints!");

	switch (get_coordtype(df))
	{
   		case CoordType::decdeg:
			Rcout << WayPoint<CoordType::decdeg>(df);
			break;

		case CoordType::degmin:
			Rcout << WayPoint<CoordType::degmin>(df);
			break;

		case CoordType::degminsec:
			Rcout << WayPoint<CoordType::degminsec>(df);
			break;

		default:
			stop("printwaypoint(DataFrame&) my bad");
	}
	return df;
}


/// __________________________________________________
/// Validate waypoints vector
// [[Rcpp::export(name = "validate.waypoints")]]
const DataFrame validatewaypoint(DataFrame& df)
{
//	cout << "——Rcpp::export——validatewaypoint(DataFrame&) format " << get_fmt_attribute(df) << endl;
	checkinherits(df, "waypoints");
	return validate(df);
}

/// __________________________________________________
/// __________________________________________________