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

class FamousFive;
class FamousFiveDD;
class FamousFiveDM;
class FamousFiveDMS;

template<class FF>
class Convert;

template<class FF>
class ConvertDD;
using convert_dd_dd = ConvertDD<FamousFiveDD>;
using convert_dm_dd = ConvertDD<FamousFiveDM>;
using convert_dms_dd = ConvertDD<FamousFiveDMS>;

template<class FF>
class ConvertDM;
using convert_dd_dm = ConvertDM<FamousFiveDD>;
using convert_dm_dm = ConvertDM<FamousFiveDM>;
using convert_dms_dm = ConvertDM<FamousFiveDMS>;

template<class FF>
class ConvertDMS;
using convert_dd_dms = ConvertDMS<FamousFiveDD>;
using convert_dm_dms = ConvertDMS<FamousFiveDM>;
using convert_dms_dms = ConvertDMS<FamousFiveDMS>;

template<class FF>
class Format;

template<class FF>
class FormatDD;
using format_decdeg = FormatDD<FamousFiveDD>;

template<class FF>
class FormatDM;
using format_degmin = FormatDM<FamousFiveDM>;

template<class FF>
class FormatDMS;
using format_degminsec = FormatDMS<FamousFiveDMS>;

class Coord;
class newValidator;

class CoordBase;
ostream& operator<<(ostream&, const CoordBase&);

template<class FF>
class Validator;

class DecDeg;
class DegMin;
class DegMinSec;

class FormatLL;
class FormatLL_DD;
template<class FF>
class FormatLL_DM_S;
using formatll_dm = FormatLL_DM_S<FamousFiveDM>;
using formatll_dms = FormatLL_DM_S<FamousFiveDMS>;

template<class T>
unique_ptr<const Coord> newconstCoord(const T&, const CoordType);

template<class T>
unique_ptr<const CoordBase> newconstCoordBase(const T&, const CoordType);

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

enum class CoordType : char { decdeg, degmin, degminsec };

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
/// FamousFive abstract base class
class FamousFive {
	public:
		virtual int get_deg(double) const = 0;
		virtual double get_decdeg(double) const = 0;
		virtual int get_min(double) const = 0;
		virtual double get_decmin(double) const = 0;
		virtual double get_sec(double) const = 0;
};


/// __________________________________________________
/// Decimal degrees derived class
class FamousFiveDD : public FamousFive {
	public:
		int get_deg(double x) const { return int(x); }
		double get_decdeg(double x) const { return x; }
		int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
		double get_decmin(double x) const { return polish(mod1by60(x)); }
		double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};


/// __________________________________________________
/// Degrees and minutes derived class
class FamousFiveDM : public FamousFive {
	public:
		int get_deg(double x) const { return int(x / 1e2); }
		double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
		int get_min(double x) const { return int(x) % int(1e2); }
		double get_decmin(double x) const { return polish(mod1e2(x)); }
		double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};


/// __________________________________________________
/// Degrees minutes and seconds derived class
class FamousFiveDMS : public FamousFive {
	public:
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
/// Formatting functor base class
template<class FF>
class Format {
	protected:
		const FF ff;
		ostringstream outstrstr;
	public:
/*		Format()
		{
			cout << "§Format<class FF>() "; _ctrsgn(typeid(*this));
		} */
		Format() = default;
		Format(const Format&) = delete;				// Disallow copying
		Format& operator=(const Format&) = delete;	//  ——— ditto ———
		Format(Format&&) = delete;					// Disallow transfer ownership
		Format& operator=(Format&&) = delete;		// Disallow moving
		virtual ~Format() = 0;
};

template<class FF>
inline Format<FF>::~Format() { /* cout << "§Format::~Format() "; _ctrsgn(typeid(*this), true); */ }


/// __________________________________________________
/// Formatting functor for decimal degrees [rescue version]
template<class FF>
class FormatDD : public Format<FF> {
	public:
/*		FormatDD()
		{
			cout << "§FormatDD<class FF>() "; _ctrsgn(typeid(*this));
		} */
		FormatDD() = default;
	    using Format<FF>::outstrstr;
	    using Format<FF>::ff;
		string operator()(double n)
		{
		//	cout << "@FormatDD::operator()\n";
			outstrstr.str("");
			outstrstr << setw(11) << setfill(' ')  << fixed << setprecision(6) << ff.get_decdeg(n) << "\u00B0";
			return outstrstr.str();
		}
};


/// __________________________________________________
/// Formatting functor for degrees and minutes
template<class FF>
class FormatDM : public Format<FF> {
	public:
/*		FormatDM()
		{
			cout << "§FormatDM<class FF>() "; _ctrsgn(typeid(*this));
		} */
		FormatDM() = default;
	    using Format<FF>::outstrstr;
	    using Format<FF>::ff;
		string operator()(double n)
		{
		//	cout << "@FormatDM::operator()\n";
			outstrstr.str("");
			outstrstr << setw(3) << setfill(' ') << abs(ff.get_deg(n)) << "\u00B0"
					  << setw(7) << setfill('0') << fixed << setprecision(4) << abs(ff.get_decmin(n)) << "'";
			return outstrstr.str();
		}
};


/// __________________________________________________
/// Formatting functor for degrees, minutes and seconds
template<class FF>
class FormatDMS : public Format<FF> {
	public:
/*		FormatDMS()
		{
			cout << "§FormatDMS<class FF>() "; _ctrsgn(typeid(*this));
		} */
		FormatDMS() = default;
	    using Format<FF>::outstrstr;
	    using Format<FF>::ff;
		string operator()(double n)
		{
		//	cout << "@FormatDMS::operator()\n";
			outstrstr.str("");
			outstrstr << setw(3) << setfill(' ') << abs(ff.get_deg(n)) << "\u00B0"
					  << setw(2) << setfill('0') << abs(ff.get_min(n)) << "'"
					  << setw(5) << fixed << setprecision(2) << abs(ff.get_sec(n)) << "\"";
			return outstrstr.str();
		}
};


/// __________________________________________________
/// __________________________________________________
/// Coordinate class
class Coord {
	protected:
		vector<double> nv;
		unique_ptr <FamousFive> ff;
		const vector<bool> valid { false };
		const vector<bool> latlon;
		const vector<string> names;
		const bool llgt1 = false;
		bool all_valid() const;
		bool waypoint = false;

	public:
		template<class FF>
		Coord(const vector<double>, in_place_type_t<FF>, const vector<bool>&, const vector<string>&);
		template<class FF>
		Coord(const NumericVector&, in_place_type_t<FF>);
		template<class CV>
//		explicit Coord(const Coord&, in_place_type_t<CV>);
		Coord& operator=(const Coord&) = delete;

		void validate(bool) const;
		const vector<double>& get_nv() const;
		const vector<bool>& get_valid() const;
		const vector<string>& get_names() const;

		void warn_invalid() const;
		void set_waypoint() const;
		vector<string> format() const;
		void print(ostream&) const;

		friend class FormatLL;
		friend class FormatLL_DD;
		template<class FF>
		friend class FormatLL_DM_S;
		friend class newValidator;

};


template<class FF>
Coord::Coord(const vector<double> n, in_place_type_t<FF>, const vector<bool>& ll, const vector<string>& _names) :
	nv(std::move(n)), ff(std::move(unique_ptr<FF>(new FF))), latlon{ ll }, names{ std::move(_names) }, llgt1(latlon.size() > 1)
{
	cout << "§Coord::Coord(const vector<double>, const LogicalVector&, const vector<string>&) "; _ctrsgn(typeid(*this));
}


template<class FF>
Coord::Coord(const NumericVector& nv, in_place_type_t<FF>) :
	Coord(
		as<vector<double>>(nv), in_place_type<FF>,
		nv.hasAttribute("latlon") ? as<vector<bool>>(nv.attr("latlon")) : vector<bool>(),
		nv.hasAttribute("names") ? as<vector<string>>(nv.attr("names")) : vector<string>()
	)
{
	cout << "§Coord::Coord(const NumericVector&, in_place_type_t<FF>) "; _ctrsgn(typeid(*this));
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

class newValidator {
		const Coord& c; 
		vector<bool>::const_iterator ll_it;
	public:
		newValidator(const Coord& _c) : c(_c), ll_it(c.latlon.begin())
		{
			cout << "§newValidator(const Coord&) "; _ctrsgn(typeid(*this));
		}
		bool operator()(double n)
		{
			cout << "@newValidator() " << " n: " << setw(9) << setfill(' ') << n << endl;
			return !((abs(c.ff->get_decdeg(n)) > (c.latlon.size() && (c.llgt1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(abs(c.ff->get_decmin(n)) >= 60) ||
				(abs(c.ff->get_sec(n)) >= 60));
		}
};


/// __________________________________________________
/// Validate coords vector
void Coord::validate(bool warn = true) const
{
	cout << "@Coord::validate() " << typeid(*this).name() << " latlon " << LogicalVector(wrap(latlon)) << endl;
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
bool Coord::all_valid() const
{
//	cout << "@Coord::all_valid()\n";
	return all_of(valid.begin(), valid.end(), [](bool v) { return v;});
}


/// __________________________________________________
/// Get const reference to nv
inline const vector<double>& Coord::get_nv() const
{
//	cout << "@Coord::get_nv()\n";
	return nv;
}


/// __________________________________________________
/// Get const reference to valid
inline const vector<bool>& Coord::get_valid() const
{
	return valid;
}


/// __________________________________________________
/// Get const reference to names
inline const vector<string>& Coord::get_names() const
{
	return names;
}


/// __________________________________________________
/// Warn if any valid are false
void Coord::warn_invalid() const
{
//	cout << "@Coord::warn_invalid()\n";
	if (!all_valid())
		warning("Validation failed!");
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
/// __________________________________________________
/// Formatting functors for Latitude and Longitude

/// __________________________________________________
/// FormatLL functor base class
class FormatLL {
	protected:
		const Coord& c; 
		vector<bool>::const_iterator ll_it;
		ostringstream outstrstr;
	public:
		FormatLL(const Coord& _c) : c(_c), ll_it(c.latlon.begin())
		{
			cout << "§FormatLL(const Coord&) "; _ctrsgn(typeid(*this));
		}
		FormatLL(const FormatLL&) = delete;				// Disallow copying
		FormatLL& operator=(const FormatLL&) = delete;	//  ——— ditto ———
		FormatLL(FormatLL&&) = delete;					// Disallow transfer ownership
		FormatLL& operator=(FormatLL&&) = delete;	    // Disallow moving
		virtual ~FormatLL() = 0;
};

inline FormatLL::~FormatLL() { /* cout << "§FormatLL::~FormatLL() "; _ctrsgn(typeid(*this), true); */ }

/// __________________________________________________
/// FormatLL functor for latitude and longitude strings for decimal degrees
class FormatLL_DD : public FormatLL {
	public:
		FormatLL_DD(const Coord& c) : FormatLL(c)
		{
			cout << "§FormatLL_DD<FF>(const Coord&) "; _ctrsgn(typeid(*this));
		}
		string operator()(string ostr, double n)
		{
			cout << "@FormatLL_DD::operator() c.waypoint " << boolalpha << c.waypoint << endl;
			if (c.latlon.size() && !c.waypoint)
				return ostr += ((c.llgt1 ? *ll_it++ : *ll_it) ? " lat" : " lon");
			else
				return ostr;
		}
};


/// __________________________________________________
/// FormatLL functor for latitude and longitude strings for degrees, minutes (and seconds)
template<class FF>
class FormatLL_DM_S : public FormatLL {
	public:
		FormatLL_DM_S(const Coord& c) : FormatLL(c)
		{
			cout << "§FormatLL_DM_S<FF>(const Coord&) "; _ctrsgn(typeid(*this));
		}
		string operator()(string ostr, double n)
		{
			cout << "@FormatLL_DM_S::operator()\n";
			return ostr += c.latlon.size() ? cardpoint(c.ff->get_decmin(n) < 0, c.llgt1 ? *ll_it++ : *ll_it) : cardi_b(c.ff->get_decmin(n) < 0);
		}
};


/// __________________________________________________
/// Formatted coordinate strings for printing
vector<string> Coord::format() const
{
//	cout << "@Coord::format<FT, FL>()\n";
	cout << "@Coord::format()\n";
	vector<string> out(nv.size());
//	transform(nv.begin(), nv.end(), out.begin(), FT());
//	transform(out.begin(), out.end(), nv.begin(), out.begin(), FL(*this));
	transform(nv.begin(), nv.end(), out.begin(), format_degmin());
	transform(out.begin(), out.end(), nv.begin(), out.begin(), formatll_dm(*this));
	return out;
}


/// __________________________________________________
/// Print coords vector
void Coord::print(ostream& stream) const
{
//	cout << "@Coord::print() type " << typeid(*this).name() << endl;
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
ostream& operator<<(ostream& stream, const Coord& c)
{
//	cout << "@operator<<(ostream&, const Coord&)\n";
	c.print(stream);
	return stream;
}


/// __________________________________________________
/// __________________________________________________
/// Create unique_ptr<Coord> to new Coord object
template<class T>
unique_ptr<const Coord> newconstCoord(const T& t, const CoordType type)
{
	cout << "@newconstCoord<T>(const T&, const CoordType) of type " << coordtype_to_int(type) + 1 << endl;

	switch (type)
	{
		case CoordType::decdeg:
					return factory<const Coord>(t, in_place_type<FamousFiveDD>);

		case CoordType::degmin:
					return factory<const Coord>(t, in_place_type<FamousFiveDM>);

		case CoordType::degminsec:
					return factory<const Coord>(t, in_place_type<FamousFiveDMS>);
		default:
					stop("newconstCoord(const T&, const CoordType) my bad");
	}
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
		validatecoord(nv);
		return check_valid(nv);
	}
}


/// __________________________________________________
/// Validate coords vector
vector<bool> validatecoord(const NumericVector& nv)
{
//	cout << "@validatecoord()\n";
	unique_ptr<const Coord> c{ newconstCoord(nv, get_coordtype(nv)) };
	c->validate();
	const_cast<NumericVector&>(nv).attr("valid") = c->get_valid();
	return c->get_valid();
}


/// __________________________________________________
/// __________________________________________________
/// Exported functions


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
	unique_ptr<const Coord> cb1{ newconstCoord(nv, inheritscoords ? oldtype : newtype) };
	if (inheritscoords) {
//		unique_ptr<const CoordBase> cb2{ cb1->convert(newtype) };
		unique_ptr<const Coord> cb2{ newconstCoord(nv, inheritscoords ? oldtype : newtype) };  // Placeholder for compiler…
		cb1.swap(cb2);
		copy((cb1->get_nv()).begin(), (cb1->get_nv()).end(), nv.begin());
	} else {
		nv.attr("class") = "coords";
	}
	nv.attr("fmt") = fmt;
//	cb1->validate_tmpl();
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
	validatecoord(nv);
	return nv;
}


/// __________________________________________________
/// Print coords vector - S3 method print.coords()      /////// "invisible" not working ///////
// [[Rcpp::export(name = "print.coords", invisible = true)]]
NumericVector printcoord(NumericVector& nv)
{
//	cout << "——Rcpp::export——printcoord() format " << get_fmt_attribute(nv) << endl;
	checkinherits(nv, "coords");
	if (!check_valid(nv))
		warning("Printing invalid coords!");
//	Rcout << *newconstCoordBase(nv, get_coordtype(nv)) << endl;
	Rcout << *newconstCoord(nv, get_coordtype(nv)) << endl;
	return nv;
}


/// __________________________________________________
/// Validate coords vector
// [[Rcpp::export(name = "validate.coords")]]
vector<bool> Rvalidatecoord(NumericVector& nv)
{
//	cout << "——Rcpp::export——Rvalidatecoord()\n";
	checkinherits(nv, "coords");
	return validatecoord(nv);
}

/*
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