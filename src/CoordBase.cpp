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

template<class T, class A1>
inline unique_ptr<T> factory(A1&&);
template<class T, class A1, class A2>
inline unique_ptr<T> factory(A1&&, A2&&);


enum class CoordType : char;
inline const CoordType get_coordtype(int);
template<class T>
inline const CoordType get_coordtype(const T&);
inline const int coordtype_to_int(CoordType);

inline string cardpoint(bool, bool);
inline string cardi_b(bool);

enum class CoordType : char { decdeg, degmin, degminsec };

template<CoordType>
struct FamousFive;
template<>
struct FamousFive<CoordType::decdeg>;
template<>
struct FamousFive<CoordType::degmin>;
template<>
struct FamousFive<CoordType::degminsec>;

template<class FF>
class Convert;

template<class FF>
class ConvertDD;
using convert_dd_dd = ConvertDD<FamousFive<CoordType::decdeg>>;
using convert_dm_dd = ConvertDD<FamousFive<CoordType::degmin>>;
using convert_dms_dd = ConvertDD<FamousFive<CoordType::degminsec>>;

template<class FF>
class ConvertDM;
using convert_dd_dm = ConvertDM<FamousFive<CoordType::decdeg>>;
using convert_dm_dm = ConvertDM<FamousFive<CoordType::degmin>>;
using convert_dms_dm = ConvertDM<FamousFive<CoordType::degminsec>>;

template<class FF>
class ConvertDMS;
using convert_dd_dms = ConvertDMS<FamousFive<CoordType::decdeg>>;
using convert_dm_dms = ConvertDMS<FamousFive<CoordType::degmin>>;
using convert_dms_dms = ConvertDMS<FamousFive<CoordType::degminsec>>;

template<CoordType type>
class Coord;
template<CoordType type>
class newValidator;

template<CoordType type>
class Format;

template<CoordType type>
class FormatLL;

class WayPoint;
template<class T>
unique_ptr<const WayPoint> newconstWayPoint(const T&);
ostream& operator<<(ostream&, const WayPoint&);

template<class T>
inline int get_fmt_attribute(const T&);
template<class T>
inline void checkinherits(T&, const char*);
template<class T, class V>
void colattrset(const T&, int, const char*, V&&);

inline bool validcoord(NumericVector&);
inline LogicalVector get_valid(const NumericVector&);
template<class T>
bool check_valid(const T&);
template<>
bool check_valid<NumericVector>(const NumericVector&);
template<CoordType type>
vector<bool> validatecoord(const NumericVector&);

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
/// Generic factory functions returning unique_ptr
/// see: Perfect Forwarding (http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2006/n2027.html)

template<class T, class A1>
inline unique_ptr<T>
factory(A1&& a1)   // one argument version
{
///§
	cout << "@factory(A1& a1)\n";
	return unique_ptr<T>(new T(std::forward<A1>(a1)));
}

template<class T, class A1, class A2>
inline unique_ptr<T>
factory(A1&& a1, A2&& a2)   // two argument version
{
///§
	cout << "@factory(A1& a1, A2& a2)\n";
	return unique_ptr<T>(new T(std::forward<A1>(a1), std::forward<A2>(a2)));
}


/// __________________________________________________
/// __________________________________________________
/// CoordType enum class

// enum class CoordType : char { decdeg, degmin, degminsec };

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
	FamousFive()
	{
		cout << "§FamousFive<CoordType::decdeg>() "; _ctrsgn(typeid(*this));
	}	
	int get_deg(double x) const { return int(x); }
	double get_decdeg(double x) const { return x; }
	int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
	double get_decmin(double x) const { return polish(mod1by60(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

template<>
struct FamousFive<CoordType::degmin> {
	FamousFive()
	{
		cout << "§FamousFive<CoordType::degmin>() "; _ctrsgn(typeid(*this));
	}	
	int get_deg(double x) const { return int(x / 1e2); }
	double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
	int get_min(double x) const { return int(x) % int(1e2); }
	double get_decmin(double x) const { return polish(mod1e2(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

template<>
struct FamousFive<CoordType::degminsec> {
	FamousFive()
	{
		cout << "§FamousFive<CoordType::degminsec>() "; _ctrsgn(typeid(*this));
	}	
	int get_deg(double x) const { return int(x / 1e4); }
	double get_decdeg(double x) const { return int(x / 1e4) + (double)int(fmod(x, 1e4) / 1e2) / 60 + mod1e2(x) / 3600; }
	int get_min(double x) const { return (int(x) % int(1e4)) / 1e2; }
	double get_decmin(double x) const { return int(fmod(x, 1e4) / 1e2) + mod1e2(x) / 60; }
	double get_sec(double x) const { return mod1e2(x); }
};



/// __________________________________________________
/// __________________________________________________
/// Convert functor base class
template<class FF>
class Convert {
	protected:
//		const FamousFive& ff; 
		const FF ff; 
	public:
/*		Convert()
		{
			cout << "§Convert<class FF>() "; _ctrsgn(typeid(*this));
		} */
		Convert() = default;
		Convert(const Convert&) = delete;				// Disallow copying
		Convert& operator=(const Convert&) = delete;	//  ——— ditto ———
		Convert(Convert&&) = delete;					// Disallow transfer ownership
		Convert& operator=(Convert&&) = delete;			// Disallow moving
		virtual ~Convert() = 0;
};

template<class FF>
inline Convert<FF>::~Convert() { /* cout << "§Convert::~Convert() "; _ctrsgn(typeid(*this), true); */ }

/// __________________________________________________
/// Convert functor for decimal degrees
template<class FF>
class ConvertDD : public Convert<FF> {
	public:
		ConvertDD()
		{
			cout << "§ConvertDD<class FF>() "; _ctrsgn(typeid(*this));
		}
	    using Convert<FF>::ff;
		double operator()(double n) { cout << "@ConvertDD::operator()(double)\n"; return ff.get_decdeg(n); }
};


/// __________________________________________________
/// Convert functor for degrees and minutes
template<class FF>
class ConvertDM : public Convert<FF> { 
	public:
		ConvertDM()
		{
			cout << "§ConvertDM<class FF>() "; _ctrsgn(typeid(*this));
		}
	    using Convert<FF>::ff;
		double operator()(double n) { cout << "@ConvertDM::operator()(double)\n"; return ff.get_deg(n) * 1e2 + ff.get_decmin(n); }
};


/// __________________________________________________
/// Convert functor for degrees, minutes and seconds
template<class FF>
class ConvertDMS : public Convert<FF> { 
	public:
		ConvertDMS()
		{
			cout << "§ConvertDMS<class FF>() "; _ctrsgn(typeid(*this));
		}
	    using Convert<FF>::ff;
		double operator()(double n) { cout << "@ConvertDMS::operator()(double)\n"; return ff.get_deg(n) * 1e4 + ff.get_min(n) * 1e2 + ff.get_sec(n); }
};


/// __________________________________________________
/// __________________________________________________
/// Coordinate class
template<CoordType type>
class Coord {
	protected:
		vector<double> nv;
//		unique_ptr<FamousFive> ff;
//		FamousFive<CoordType::degmin> fly5;
		FamousFive<type> ff;
		const vector<bool> valid { false };
		const vector<bool> latlon;
		const vector<string> names;
		const bool llgt1 = false;
		bool all_valid() const;
		bool waypoint = false;

	public:
		Coord(const vector<double>, const vector<bool>&, const vector<string>&);
		Coord(const NumericVector&);
//		template<class CV>
//		explicit Coord(const Coord&, in_place_type_t<CV>);
		Coord<type>& operator=(const Coord<type>&) = delete;
		~Coord<type>() {
			cout << "§Coord<type>::~Coord() "; _ctrsgn(typeid(*this), true);		
		}

		void validate(bool = true) const;
		const vector<double>& get_nv() const;
		const vector<bool>& get_valid() const;
		const vector<string>& get_names() const;

		void warn_invalid() const;
		void set_waypoint() const;
		vector<string> format() const;
		void print(ostream&) const;

		friend class Format<type>;
		friend class FormatLL<type>;
		friend class newValidator<type>;

};


template<CoordType type>
Coord<type>::Coord(const vector<double> n, const vector<bool>& ll, const vector<string>& _names) :
	nv(std::move(n)), latlon{ ll }, names{ std::move(_names) }, llgt1(latlon.size() > 1)
{
	cout << "§Coord<type>::Coord(const vector<double>, const vector<bool>&, const vector<string>&) latlon.size() " << latlon.size() << " "; _ctrsgn(typeid(*this));
}


template<CoordType type>
Coord<type>::Coord(const NumericVector& nv) :
	Coord(
		as<vector<double>>(nv),
		nv.hasAttribute("latlon") ? as<vector<bool>>(nv.attr("latlon")) : vector<bool>(),
		nv.hasAttribute("names") ? as<vector<string>>(nv.attr("names")) : vector<string>()
	)
{
	cout << "§Coord<type>::Coord(const NumericVector&) "; _ctrsgn(typeid(*this));
}

/*
template<class CV>
Coord::Coord(const Coord& c, in_place_type_t<CV>) :
	Coord(vector<double>(c.nv.size()), vector<bool>{ c.latlon }, vector<string>{ c.names })
{
	cout << "§Coord::Coord<Convert_type>(const Coord&, in_place_type_t<Convert_type>) "; _ctrsgn(typeid(*this));
//	using convert_ff = CV<c::FF>;
	transform(c.nv.begin(), c.nv.end(), nv.begin(), CV());
} */


/// __________________________________________________
/// __________________________________________________
/// Validate Coord functor

template<CoordType type>
class newValidator {
		const Coord<type>& c; 
		vector<bool>::const_iterator ll_it;
	public:
		newValidator(const Coord<type>& _c) : c(_c), ll_it(c.latlon.begin())
		{
			cout << "§newValidator(const Coord<type>&) "; _ctrsgn(typeid(*this));
		}
		bool operator()(double n)
		{
			cout << "@newValidator() " << " n: " << setw(9) << setfill(' ') << n << endl;
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
	cout << "@Coord<type>::validate() " << typeid(*this).name() << " latlon " << LogicalVector(wrap(latlon)) << endl;
	vector<bool>& non_const_valid { const_cast<vector<bool>&>(valid) };
	non_const_valid.assign(nv.size(), {false});
	transform(nv.begin(), nv.end(), non_const_valid.begin(), newValidator(*this));
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
		Format(const Coord<type>&);
		Format(const Format&) = delete;				// Disallow copying
		Format& operator=(const Format&) = delete;	//  ——— ditto ———
		Format(Format&&) = delete;					// Disallow transfer ownership
		Format& operator=(Format&&) = delete;		// Disallow moving
		~Format();
		string operator()(double n);
};

template<CoordType type>
Format<type>::Format(const Coord<type>& _c) : c(_c)
{
	cout << "§Format(const Coord<type>&) "; _ctrsgn(typeid(*this));
}

template<CoordType type>
Format<type>::~Format()
{
	cout << "§Format::~Format() "; _ctrsgn(typeid(*this), true);		
}

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
		FormatLL(const Coord<type>& _c);
		FormatLL(const FormatLL&) = delete;				// Disallow copying
		FormatLL& operator=(const FormatLL&) = delete;	//  ——— ditto ———
		FormatLL(FormatLL&&) = delete;					// Disallow transfer ownership
		FormatLL& operator=(FormatLL&&) = delete;		// Disallow moving
		~FormatLL();
		string operator()(string, double);
};

template<CoordType type>
FormatLL<type>::FormatLL(const Coord<type>& _c) : c(_c), ll_it(c.latlon.begin())
{
	cout << "§FormatLL(const Coord<type>&) "; _ctrsgn(typeid(*this));
}

template<CoordType type>
FormatLL<type>::~FormatLL()
{
	cout << "§FormatLL::~FormatLL() "; _ctrsgn(typeid(*this), true);		
}

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
	cout << "@Coord<type>::format()\n";
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
	cout << "@Coord<type>::print() type " << typeid(*this).name() << endl;
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
	cout << "@operator<<(ostream&, const Coord<type>&)\n";
	c.print(stream);
	return stream;
}


/// __________________________________________________
/// __________________________________________________
/// Convenience functions

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
void colattrset(const T& t, int col, const char* attrib, V&& val)
{
//	cout << "@colattrset(const T&, int, const char*, V&&) attrib " << attrib << ", col " << col << endl;
	as<NumericVector>(t[col]).attr(attrib) = std::forward<V>(val);
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
/// Return "valid" attribute or empty LogicalVector    !!!!!!! Generalise with template and specialisation !!!!!!!
inline LogicalVector get_valid(const NumericVector& nv)
{
//	cout << "@get_valid(const NumericVector&) has attr \"valid\" " << boolalpha << nv.hasAttribute("valid") << endl;
	return (nv.hasAttribute("valid") ? LogicalVector(nv.attr("valid")) : LogicalVector());
}


/// __________________________________________________
/// Check first two columns have attribute "valid" all true
template<class T>
bool check_valid(const T& t)
{
//	cout << "@check_valid(const T&) fmt " << get_fmt_attribute(t) << endl;
	bool boolat = check_valid(as<NumericVector>(t[1]));
	if (!boolat)
		warning("Invalid latitude!");
	bool boolon = check_valid(as<NumericVector>(t[2]));
	if (!boolat)
		warning("Invalid longitude!");
	return boolat || boolon;
}


/// __________________________________________________
/// Specialisation for <NumericVector>
template<>
bool check_valid<NumericVector>(const NumericVector& nv)
{
//	cout << "@check_valid<NumericVector>(const NumericVector&)" << endl;
	LogicalVector valid = std::move(get_valid(nv));
	if (valid.size())
		return all_of(valid.begin(), valid.end(), [](bool v) { return v;});
	else {
		warning("Unvalidated coords! Revalidating…");
//		validatecoord(nv);
		return check_valid(nv);
	}
}


/// __________________________________________________
/// Validate coords vector
template<CoordType type>
vector<bool> validatecoord(const NumericVector& nv)
{
	cout << "@validatecoord()\n";
	Coord<type> c(nv);
	c.validate();
	const_cast<NumericVector&>(nv).attr("valid") = c.get_valid();
	return c.get_valid();
}


/// __________________________________________________
/// __________________________________________________
/// Exported functions

/*
/// __________________________________________________
/// Set R vector object class to coords and return,
/// or convert format of R coords object and return  !!!!!!! Template and Specialise [with waypoints()] !!!!!!!
// [[Rcpp::export]]
NumericVector coords(NumericVector& nv, int fmt = 1)
{
//	cout << "——Rcpp::export——coords()\n";
	CoordType newtype = get_coordtype(fmt);
	CoordType oldtype;
	const bool inheritscoords { nv.inherits("coords") };
	if (inheritscoords) {
		oldtype = get_coordtype(nv);
//		cout <<  "coords() argument nv is already a \"coords\" vector of type "
//			 << coordtype_to_int(oldtype) + 1 << endl;
		if (newtype == oldtype) {
//			cout << "——fmt out == fmt in!——" << endl;
			if (!check_valid(nv))
				warning("Invalid coords!");
			return nv;
		}
	}
//	unique_ptr<const CoordBase> cb1{ newconstCoordBase(nv, inheritscoords ? oldtype : newtype) };
//	unique_ptr<const Coord<CoordType::degmin>> cb1{ newconstCoord<CoordType::degmin>(nv, inheritscoords ? oldtype : newtype) };
	unique_ptr<const Coord<CoordType::degmin>> cb1{ newconstCoord<CoordType::degmin>(nv, CoordType::degmin) };
	if (inheritscoords) {
//		unique_ptr<const CoordBase> cb2{ cb1->convert(newtype) };
		unique_ptr<const Coord<CoordType::degmin>> cb2{ newconstCoord<CoordType::degmin>(nv, CoordType::degmin) };  // Placeholder for compiler…
		cb1.swap(cb2);
		copy((cb1->get_nv()).begin(), (cb1->get_nv()).end(), nv.begin());
	} else {
		nv.attr("class") = "coords";
	}
	nv.attr("fmt") = fmt;
//	cb1->validate_tmpl();
	cb1->validate();
	nv.attr("valid") = cb1->get_valid();
	return nv;
}


/// __________________________________________________
/// coords() as replacement function
// [[Rcpp::export(name = "`coords<-`")]]
NumericVector coords_replace(NumericVector& nv, int value)
{
//	cout << "——Rcpp::export——`coords_replace()<-`\n";
	return coords(nv, value);
} */


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
//	validatecoord(nv);
	return nv;
}


/// __________________________________________________
/// Print coords vector - S3 method print.coords()      /////// "invisible" not working ///////
// [[Rcpp::export(name = "print.coords", invisible = true)]]
NumericVector printcoord(NumericVector& nv)
{
	cout << "——Rcpp::export——printcoord() format " << get_fmt_attribute(nv) << endl;
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

/*
/// __________________________________________________
/// Validate coords vector
// [[Rcpp::export(name = "validate.coords")]]
vector<bool> Rvalidatecoord(NumericVector& nv)
{
//	cout << "——Rcpp::export——Rvalidatecoord()\n";
	checkinherits(nv, "coords");
	return validatecoord(nv);
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
	return newconstCoordBase(nv, get_coordtype(nv))->fmt_fctr_tmpl();
}


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


/// __________________________________________________
/// Add "waypoints" to R data.frame object class and validate. !!!!!!! Template and Specialise !!!!!!!
// [[Rcpp::export]]
DataFrame waypoints(DataFrame& df, int fmt = 1)
{
//	cout << "——Rcpp::export——waypoints()\n";
	bool boolio;
	vector<int> llcols { 1, 2 };
	CoordType newtype = get_coordtype(fmt);
	CoordType oldtype;
	const bool inheritswaypoints { df.inherits("waypoints") };
	if (inheritswaypoints) {
		oldtype = get_coordtype(df);
//		cout << "argument df is already a \"waypoints\" vector of type " << coordtype_to_int(oldtype) + 1 << endl;
		if (!check_valid(df))
			stop("Invalid waypoints!");
		if (newtype == oldtype) {
//			cout << "——fmt out == fmt in!——" << endl;
			return df;
		}
	} else {
		df.attr("fmt") = fmt;
		for (const auto x : llcols)
			colattrset(df, x, "latlon", vector<bool>(1, llcols[1] - x));
	}
	colattrset(df, llcols[0], "names", df[0]);				// !!!!!!!! Temporary Solution !!!!!!
	unique_ptr<const WayPoint> wp1{ newconstWaypoint(df) };
	if (inheritswaypoints) {
		unique_ptr<const WayPoint> wp2{ wp1->convert(newtype) };
		wp1.swap(wp2);
		for (const auto x : llcols) {
			boolio = llcols[2] - x;
			copy((wp1->get_cbp(boolio).get_nv()).begin(), (wp1->get_cbp(boolio).get_nv()).end(), as<NumericVector>(df[x]).begin());
		}
	} else {
		df.attr("class") = CharacterVector{"waypoints", "data.frame"};
	}
	wp1->validate(true);
	wp1->warn_invalid();
	for (const auto x : llcols) {
		boolio = llcols[2] - x;
		colattrset(df, x, "valid", LogicalVector(wrap(wp1->get_valid(boolio))));
		colattrset(df, x, "fmt", fmt);
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
/// Print waypoints vector - S3 method print.waypoints()      /////// "invisible" not working ///////
// [[Rcpp::export(name = "print.waypoints", invisible = true)]]
DataFrame printwaypoint(DataFrame& df)
{
//	cout << "——Rcpp::export——printwaypoint() format " << get_fmt_attribute(df) << endl;
	checkinherits(df, "waypoints");
	if (!check_valid(df))
		warning("Invalid waypoints!");
	Rcout << *newconstWaypoint(df) << endl;
	return df;
}


/// __________________________________________________
/// Validate waypoints vector
// [[Rcpp::export(name = "validate.waypoints")]]
const DataFrame validatewaypoint(DataFrame& df)
{
//	cout << "——Rcpp::export——validatewaypoint()\n";
	checkinherits(df, "waypoints");
	unique_ptr<const WayPoint> wp{ newconstWaypoint(df) };
	wp->validate(true);
	wp->warn_invalid();
	vector<int> llcols { 1, 2 };
	for (const auto x : llcols)
		colattrset(df, x, "valid", wp->get_valid(llcols[2] - x));
	return df;
}
*/
/// __________________________________________________
/// __________________________________________________