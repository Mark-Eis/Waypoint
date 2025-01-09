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

template <class T, class A1>
inline unique_ptr<T> factory(A1&&);
template <class T, class A1, class A2>
inline unique_ptr<T> factory(A1&&, A2&&);


enum class CoordType : char;
inline const CoordType get_coordtype(int);
template<class T>
inline const CoordType get_coordtype(const T&);
inline const int coordtype_to_int(CoordType);

inline string cardpoint(bool, bool);
inline string cardi_b(bool);

struct FamousFiveDD;
struct FamousFiveDM;
struct FamousFiveDMS;

template <class FamousFive_type>
class ConvertDD;
template <class FamousFive_type>
class ConvertDM;
template <class FamousFive_type>
class ConvertDMS;

template <class FamousFive_type>
class Format;

template <class FamousFive_type>
class FormatDD;
using format_decdeg = FormatDD<FamousFiveDD>;

template <class FamousFive_type>
class FormatDM;
using format_degmin = FormatDM<FamousFiveDM>;

template <class FamousFive_type>
class FormatDMS;
using format_degminsec = FormatDMS<FamousFiveDMS>;

class CoordBase;
ostream& operator<<(ostream&, const CoordBase&);

template<class FamousFive_type>
class Validator;

class DecDeg;
class DegMin;
class DegMinSec;

class FormatLL;
class FormatLL_DD;
template <class FamousFive_type>
class FormatLL_DM_S;

template<class T>
unique_ptr<const CoordBase> newconstCoordBase(const T&, const CoordType);

class WayPoint;
template <class T>
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

template <class T, class A1>
inline unique_ptr<T>
factory(A1&& a1)   // one argument version
{
///§
	cout << "@factory(A1& a1)\n";
	return unique_ptr<T>(new T(std::forward<A1>(a1)));
}

template <class T, class A1, class A2>
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
/// Famous five functions class for decimal degrees
struct FamousFiveDD {
		int get_deg(double x) const { return int(x); }
		double get_decdeg(double x) const { return x; }
		int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
		double get_decmin(double x) const { return polish(mod1by60(x)); }
		double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};


/// __________________________________________________
/// Famous five functions class for degrees and minutes
struct FamousFiveDM {
		int get_deg(double x) const { return int(x / 1e2); }
		double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
		int get_min(double x) const { return int(x) % int(1e2); }
		double get_decmin(double x) const { return polish(mod1e2(x)); }
		double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};


/// __________________________________________________
/// Famous five functions class for degrees minutes and seconds
struct FamousFiveDMS {
		int get_deg(double x) const { return int(x / 1e4); }
		double get_decdeg(double x) const { return int(x / 1e4) + (double)int(fmod(x, 1e4) / 1e2) / 60 + mod1e2(x) / 3600; }
		int get_min(double x) const { return (int(x) % int(1e4)) / 1e2; }
		double get_decmin(double x) const { return int(fmod(x, 1e4) / 1e2) + mod1e2(x) / 60; }
		double get_sec(double x) const { return mod1e2(x); }
};

/// __________________________________________________
/// __________________________________________________
/// Convert functor for decimal degrees
template <class FamousFive_type>
class ConvertDD {
		const FamousFive_type ff; 
	public:
		ConvertDD()
		{
			cout << "§ConvertDD<class FamousFive_type>() "; _ctrsgn(typeid(*this));
		}
		double operator()(double n) { cout << "@ConvertDD::operator()(double)\n"; return ff.get_decdeg(n); }
};


/// __________________________________________________
/// Convert functor for degrees and minutes
template <class FamousFive_type>
class ConvertDM { 
		const FamousFive_type ff; 
	public:
		ConvertDM()
		{
			cout << "§ConvertDM<class FamousFive_type>() "; _ctrsgn(typeid(*this));
		}
		double operator()(double n) { cout << "@ConvertDM::operator()(double)\n"; return ff.get_deg(n) * 1e2 + ff.get_decmin(n); }
};


/// __________________________________________________
/// Convert functor for degrees, minutes and seconds
template <class FamousFive_type>
class ConvertDMS { 
		const FamousFive_type ff; 
	public:
		ConvertDMS()
		{
			cout << "§ConvertDMS<class FamousFive_type>() "; _ctrsgn(typeid(*this));
		}
		double operator()(double n) { cout << "@ConvertDMS::operator()(double)\n"; return ff.get_deg(n) * 1e4 + ff.get_min(n) * 1e2 + ff.get_sec(n); }
};


/// __________________________________________________
/// __________________________________________________
/// Formatting functor base class
template <class FamousFive_type>
class Format {
	protected:
		FamousFive_type ff;
		ostringstream outstrstr;
	public:
/*		Format()
		{
			cout << "§Format<class FamousFive_type>() "; _ctrsgn(typeid(*this));
		} */
		Format() = default;
		Format(const Format&) = delete;				// Disallow copying
		Format& operator=(const Format&) = delete;	//  ——— ditto ———
		Format(Format&&) = delete;					// Disallow transfer ownership
		Format& operator=(Format&&) = delete;	    // Disallow moving
		virtual ~Format() = 0;
};

template <class FamousFive_type>
inline Format<FamousFive_type>::~Format() { /* cout << "§Format::~Format() "; _ctrsgn(typeid(*this), true); */ }


/// __________________________________________________
/// Formatting functor for decimal degrees [rescue version]
template <class FamousFive_type>
class FormatDD : public Format<FamousFive_type> {
	public:
/*		FormatDD()
		{
			cout << "§FormatDD<class FamousFive_type>() "; _ctrsgn(typeid(*this));
		} */
		FormatDD() = default;
	    using Format<FamousFive_type>::outstrstr;
	    using Format<FamousFive_type>::ff;
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
template <class FamousFive_type>
class FormatDM : public Format<FamousFive_type> {
	public:
/*		FormatDM()
		{
			cout << "§FormatDM<class FamousFive_type>() "; _ctrsgn(typeid(*this));
		} */
		FormatDM() = default;
	    using Format<FamousFive_type>::outstrstr;
	    using Format<FamousFive_type>::ff;
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
template <class FamousFive_type>
class FormatDMS : public Format<FamousFive_type> {
	public:
/*		FormatDMS()
		{
			cout << "§FormatDMS<class FamousFive_type>() "; _ctrsgn(typeid(*this));
		} */
		FormatDMS() = default;
	    using Format<FamousFive_type>::outstrstr;
	    using Format<FamousFive_type>::ff;
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
/// Coordinate abstract base class

class CoordBase {
	protected:
		vector<double> nv;
		const vector<bool> valid { false };
		const vector<bool> latlon;
		const vector<string> names;
		const bool llgt1 = false;
		template<typename validate_type>
		void validate(bool) const;
		bool all_valid() const;
		bool waypoint = false;

	public:
		CoordBase(const NumericVector&);
		template <class FamousFive_type>
		explicit CoordBase(const CoordBase&, FamousFive_type&&);
		CoordBase(const vector<double>, const vector<bool>&, const vector<string>&);
		CoordBase& operator=(const CoordBase&) = delete;

		virtual ~CoordBase() = 0;
		virtual const CoordType getfmt() const = 0;
		virtual void validate_tmpl(bool = true) const = 0;
		virtual vector<string> fmt_fctr_tmpl() const = 0;

		const vector<double>& get_nv() const;
		unique_ptr<const CoordBase> convert(const CoordType) const;
		void newconvert(const CoordType) const;
		const vector<bool>& get_valid() const;
		const vector<string>& get_names() const;
		void warn_invalid() const;
		void set_waypoint() const;
		template <class Format_type, class FormatLL_type>
		vector<string> format() const;
		void print(ostream&) const;

		template<class FamousFive_type>
		friend class Validator;
		friend class DecDeg;
		friend class DegMin;
		friend class DegMinSec;
		friend class FormatLL;
		friend class FormatLL_DD;
		template<class FamousFive_type>
		friend class FormatLL_DM_S;

		friend ostream& operator<<(ostream&, const CoordBase&);
};

CoordBase::CoordBase(const NumericVector& nv) :
	CoordBase(
		as<vector<double>>(nv),
		nv.hasAttribute("latlon") ? as<vector<bool>>(nv.attr("latlon")) : vector<bool>(),
		nv.hasAttribute("names") ? as<vector<string>>(nv.attr("names")) : vector<string>()
	)
{
//
	cout << "§CoordBase::CoordBase(const NumericVector&) "; _ctrsgn(typeid(*this));
}


template <class FamousFive_type>
CoordBase::CoordBase(const CoordBase& cb, FamousFive_type&& ff) :
	CoordBase(vector<double>(cb.nv.size()), vector<bool>{ cb.latlon }, vector<string>{ cb.names })
{
//
	cout << "§CoordBase::CoordBase(const CoordBase&, FunctObj) "; _ctrsgn(typeid(*this));
	transform(cb.nv.begin(), cb.nv.end(), nv.begin(), std::move(ff));
}


CoordBase::CoordBase(const vector<double> n, const vector<bool>& ll, const vector<string>& _names) :
	nv(std::move(n)), latlon{ ll }, names{ std::move(_names) }, llgt1(latlon.size() > 1)
{
//
	cout << "§CoordBase::CoordBase(const vector<double>, const LogicalVector&, const vector<string>&) "; _ctrsgn(typeid(*this));
}


CoordBase::~CoordBase()
{
//
	cout << "§CoordBase::~CoordBase() "; _ctrsgn(typeid(*this), true);
}


/// __________________________________________________
/// Get const reference to nv
inline const vector<double>& CoordBase::get_nv() const
{
//	cout << "@CoordBase::get_nv()\n";
	return nv;
}


/// __________________________________________________
/// Get const reference to valid
inline const vector<bool>& CoordBase::get_valid() const
{
	return valid;
}


/// __________________________________________________
/// Get const reference to names
inline const vector<string>& CoordBase::get_names() const
{
	return names;
}


/// __________________________________________________
/// Convert to degrees, minutes and seconds, to degrees and minutes or to decimal degrees
inline unique_ptr<const CoordBase> CoordBase::convert(const CoordType type) const
{
//	cout << "@CoordBase::convert(const CoordType type) " << typeid(*this).name()
//       << " to type "<< coordtype_to_int(type) + 1 << endl;
	return newconstCoordBase(*this, type);
}


/// __________________________________________________
/// Convert to degrees, minutes and seconds, to degrees and minutes or to decimal degrees
inline void CoordBase::newconvert(const CoordType type) const
{
//	cout << "@CoordBase::newconvert(const CoordType type) " << typeid(*this).name()
//       << " to type "<< coordtype_to_int(type) + 1 << endl;
	

}


/// __________________________________________________
/// __________________________________________________
/// Validate coord value functor
template<class FamousFive_type>
class Validator {
		const CoordBase& cb; 
		const FamousFive_type ff;
		vector<bool>::const_iterator ll_it;
	public:
		Validator(const CoordBase& _cb) : cb(_cb), ll_it(cb.latlon.begin())
		{
			cout << "§Validator(const CoordBase&) "; _ctrsgn(typeid(*this));
		}
		bool operator()(double n)
		{
			cout << "@Validator() " << " n: " << setw(9) << setfill(' ') << n << endl;
			return !((abs(ff.get_decdeg(n)) > (cb.latlon.size() && (cb.llgt1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(abs(ff.get_decmin(n)) >= 60) ||
				(abs(ff.get_sec(n)) >= 60));
		}
};


/// __________________________________________________
/// Validate coords vector
template<typename validate_type>
void CoordBase::validate(bool warn) const
{
	cout << "@CoordBase::<typename validate_type>validate() " << typeid(*this).name() << " latlon " << LogicalVector(wrap(latlon)) << endl;
	vector<bool>& non_const_valid { const_cast<vector<bool>&>(valid) };
	non_const_valid.assign(nv.size(), {false});
	transform(nv.begin(), nv.end(), non_const_valid.begin(), Validator<validate_type>(*this));
	if (all_valid())
		non_const_valid.assign({true});
	else
		if (warn)
			warning("Validation failed!");
}


/// __________________________________________________
/// All valid are true
bool CoordBase::all_valid() const
{
//	cout << "@CoordBase::all_valid()\n";
	return all_of(valid.begin(), valid.end(), [](bool v) { return v;});
}


/// __________________________________________________
/// Warn if any valid are false
void CoordBase::warn_invalid() const
{
//	cout << "@CoordBase::warn_invalid()\n";
	if (!all_valid())
		warning("Validation failed!");
}


/// __________________________________________________
/// Set waypoint flag
inline void CoordBase::set_waypoint() const
{
//	cout << "@CoordBase::set_waypoint()\n";
	bool& wpt = const_cast<bool&>(waypoint);
	wpt = true;
}


/// __________________________________________________
/// Formatted coordinate strings for printing
template <class Format_type, class FormatLL_type>
vector<string> CoordBase::format() const
{
//	cout << "@CoordBase::format<class Format_type, class FormatLL_type>()\n";
	vector<string> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), Format_type());
	transform(out.begin(), out.end(), nv.begin(), out.begin(), FormatLL_type(*this));
	return out;
}

/*
/// __________________________________________________
/// __________________________________________________
/// Print coords vector functor
class Printor {
		ostream& stream; 
		vector<string>::const_iterator nm_it;
	public:
		Printor(ostream& ostr) : stream(ostr), nm_it(cb.latlon.begin())
		{
	//		cout << "@Printor(ostream& stream) ";
			if (names.size()) {
				vector<string>::const_iterator nm_it(names.begin());		}
		}
		string operator()(double n)
		{
	//		cout << "@Printor()\n";
			stream << s << " " << *nm_it++ << "\n"; ;
		}
};


/// __________________________________________________
/// Print coords vector
void CoordBase::print(ostream& stream) const
{
//	cout << "@CoordBase::print() type " << typeid(*this).name() << endl;
	vector<string> fmtstr(fmt_fctr_tmpl()); 
	
		for_each(fmtstr.begin(), fmtstr.end(),
			[&stream,& nm_it](const string& s) { stream << s << " " << *nm_it++ << "\n"; });
	} else
		for_each(fmtstr.begin(), fmtstr.end(), [&stream](const string& s) { stream << s << "\n"; });
} */


/// __________________________________________________
/// Print coords vector
void CoordBase::print(ostream& stream) const
{
//	cout << "@CoordBase::print() type " << typeid(*this).name() << endl;
	vector<string> sv(fmt_fctr_tmpl()); 
	if (names.size()) {
		vector<string>::const_iterator nm_it(names.begin());
		for_each(sv.begin(), sv.end(),
			[&stream,& nm_it](const string& s) { stream << s << " " << *nm_it++ << "\n"; });
	} else
		for_each(sv.begin(), sv.end(), [&stream](const string& s) { stream << s << "\n"; });
}


/// __________________________________________________
/// Output CoordBase derived object to ostream
ostream& operator<<(ostream& stream, const CoordBase& c)
{
//	cout << "@operator<<(ostream&, const CoordBase&)\n";
	c.print(stream);
	return stream;
}


/// __________________________________________________
/// Decimal degrees derived class
class DecDeg : public CoordBase {
	public:
		DecDeg(const NumericVector&);
		DecDeg(const CoordBase&);
		~DecDeg();

		const CoordType getfmt() const { return CoordType::decdeg; }
		void validate_tmpl(bool) const;
		vector<string> fmt_fctr_tmpl() const;
};


DecDeg::DecDeg(const NumericVector& nv) : CoordBase(nv)
{
//
	cout << "§DecDeg::DecDeg(NumericVector&) "; _ctrsgn(typeid(*this));
}


DecDeg::DecDeg(const CoordBase& c) : CoordBase(c, ConvertDD<FamousFiveDMS>())
{
//
	cout << "§DecDeg::DecDeg(const CoordBase&) "; _ctrsgn(typeid(*this));
}


DecDeg::~DecDeg()
{
//
	cout << "§DecDeg::~DecDeg() "; _ctrsgn(typeid(*this), true);
}


/// __________________________________________________
/// Instantiate functor template for validating decimal degrees
inline void DecDeg::validate_tmpl(bool warn) const
{
	cout << "@DecDeg::validate_tmpl(bool)\n";
	return validate<FamousFiveDD>(warn);
}


/// __________________________________________________
/// Instantiate functor template for formatting decimal degrees
inline vector<string> DecDeg::fmt_fctr_tmpl() const
{
//	cout << "@DecDeg::fmt_fctr_tmpl()\n";
	return format<format_decdeg, FormatLL_DD>();
}


/// __________________________________________________
/// Degrees and minutes derived class
class DegMin : public CoordBase {
	public:
		DegMin(const NumericVector&);
		DegMin(const CoordBase&);
		~DegMin();

		const CoordType getfmt() const { return CoordType::degmin; }
		void validate_tmpl(bool) const;
		vector<string> fmt_fctr_tmpl() const;
};


DegMin::DegMin(const NumericVector& nv) : CoordBase(nv)
{
//
	cout << "§DegMin::DegMin(NumericVector&) "; _ctrsgn(typeid(*this));
}


DegMin::DegMin(const CoordBase& c) : CoordBase(c, ConvertDM<FamousFiveDD>())
{
//
	cout << "§DegMin::DegMin(const CoordBase&) "; _ctrsgn(typeid(*this));
}

DegMin::~DegMin()
{
//
	cout << "§DegMin::~DegMin() "; _ctrsgn(typeid(*this), true);
}


/// __________________________________________________
/// Instantiate functor template for validating degrees and minutes
inline void DegMin::validate_tmpl(bool warn) const
{
	cout << "@DegMin::validate_tmpl(bool)\n";
	return validate<FamousFiveDM>(warn);
}


/// __________________________________________________
/// Instantiate functor template for formatting degrees and minutes
inline vector<string> DegMin::fmt_fctr_tmpl() const
{
//	cout << "@DegMin::fmt_fctr_tmpl()\n";
	return format<format_degmin, FormatLL_DM_S<FamousFiveDM>>();
}


/// __________________________________________________
/// Degrees minutes and seconds derived class
class DegMinSec : public CoordBase {
	public:
		DegMinSec(const NumericVector&);
		DegMinSec(const CoordBase&);
		~DegMinSec();

		const CoordType getfmt() const { return CoordType::degminsec; }
		void validate_tmpl(bool) const;
		vector<string> fmt_fctr_tmpl() const;
};


DegMinSec::DegMinSec(const NumericVector& nv) : CoordBase(nv)
{
//
	cout << "§DegMinSec::DegMinSec(NumericVector&) "; _ctrsgn(typeid(*this));
}


DegMinSec::DegMinSec(const CoordBase& c) : CoordBase(c, ConvertDMS<FamousFiveDM>())
{
//
	cout << "§DegMinSec::DegMinSec(const CoordBase&) "; _ctrsgn(typeid(*this));
}

DegMinSec::~DegMinSec()
{
//
	cout << "§DegMinSec::~DegMinSec() "; _ctrsgn(typeid(*this), true);
}


/// __________________________________________________
/// Instantiate functor template for validating degrees, minutes and seconds
inline void DegMinSec::validate_tmpl(bool warn) const
{
	cout << "@DegMinSec::validate_tmpl(bool)\n";
	return validate<FamousFiveDMS>(warn);
}


/// __________________________________________________
/// Instantiate functor template for formatting degrees, minutes and seconds
inline vector<string> DegMinSec::fmt_fctr_tmpl() const
{
//	cout << "@DegMinSec::fmt_fctr_tmpl()\n";
	return format<format_degminsec, FormatLL_DM_S<FamousFiveDMS>>();
}


/// __________________________________________________
/// __________________________________________________
/// Formatting functors for Latitude and Longitude

/// __________________________________________________
/// FormatLL functor base class
class FormatLL {
	protected:
		const CoordBase& cb; 
		vector<bool>::const_iterator ll_it;
		ostringstream outstrstr;
	public:
		FormatLL(const CoordBase& _cb) : cb(_cb), ll_it(cb.latlon.begin())
		{
		//	cout << "§FormatLL(const CoordBase&) "; _ctrsgn(typeid(*this));
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
		FormatLL_DD(const CoordBase& cb) : FormatLL(cb)
		{
		//	cout << "§FormatLL_DD<FamousFive_type>(const CoordBase&) "; _ctrsgn(typeid(*this));
		}
		string operator()(string ostr, double n)
		{
		//	cout << "@FormatLL_DD::operator() cb.waypoint " << boolalpha << cb.waypoint << endl;
			if (cb.latlon.size() && !cb.waypoint)
				return ostr += ((cb.llgt1 ? *ll_it++ : *ll_it) ? " lat" : " lon");
			else
				return ostr;
		}
};


/// __________________________________________________
/// FormatLL functor for latitude and longitude strings for degrees, minutes (and seconds)
template<class FamousFive_type>
class FormatLL_DM_S : public FormatLL {
		const FamousFive_type ff;
	public:
		FormatLL_DM_S(const CoordBase& cb) : FormatLL(cb)
		{
		//	cout << "§FormatLL_DM_S<FamousFive_type>(const CoordBase&) "; _ctrsgn(typeid(*this));
		}
		string operator()(string ostr, double n)
		{
		//	cout << "@FormatLL_DM_S::operator()\n";
			return ostr += cb.latlon.size() ? cardpoint(ff.get_decmin(n) < 0, cb.llgt1 ? *ll_it++ : *ll_it) : cardi_b(ff.get_decmin(n) < 0);
		}
};


/// __________________________________________________
/// __________________________________________________
/// Create unique_ptr<CoordBase> to new DecDeg, DegMin or DegMinSec object
template<class T>
unique_ptr<const CoordBase> newconstCoordBase(const T& t, const CoordType type)
{
	cout << "@newconstCoordBase<T>(const T&, const CoordType) of type " << coordtype_to_int(type) + 1 << endl;

	switch (type)
	{
		case CoordType::decdeg:
					return factory<const DecDeg>(t);

		case CoordType::degmin:
					return factory<const DegMin>(t);

		case CoordType::degminsec:
					return factory<const DegMinSec>(t);
		default:
					stop("newconstCoordBase<t>(const T&, const CoordType) my bad");
	}
}


/// __________________________________________________
/// __________________________________________________
/// Waypoint class

class WayPoint {
	protected:
		unique_ptr<const CoordBase> cbp_lat;
		unique_ptr<const CoordBase> cbp_lon;
		const vector<bool>& validlat;
		const vector<bool>& validlon;
	public:
		explicit WayPoint(unique_ptr<const CoordBase>, unique_ptr<const CoordBase>);
		explicit WayPoint(const WayPoint&, CoordType);
		~WayPoint();

		const CoordBase& get_cbp(bool) const;
		unique_ptr<const WayPoint> convert(const CoordType) const;
		void validate(bool) const;
		const vector<bool>& get_valid(bool) const;
		void warn_invalid() const;
		void print(ostream& stream) const;
		vector<string> format() const;
		friend ostream& operator<<(ostream&, const WayPoint&);
};


WayPoint::WayPoint(unique_ptr<const CoordBase> _cbp_lat, unique_ptr<const CoordBase> _cbp_lon) :
	cbp_lat{std::move(_cbp_lat)}, cbp_lon{std::move(_cbp_lon)},
	validlat(cbp_lat->get_valid()), validlon(cbp_lon->get_valid())
{
//
	cout << "§WayPoint(unique_ptr<const CoordBase>, unique_ptr<const CoordBase>) "; _ctrsgn(typeid(*this));
	cbp_lat->set_waypoint();
	cbp_lon->set_waypoint();
}

WayPoint::WayPoint(const WayPoint& wp, CoordType type) :
	WayPoint{ wp.get_cbp(true).convert(type), wp.get_cbp(false).convert(type) }
{
//
	cout << "§WayPoint(const WayPoint&) "; _ctrsgn(typeid(*this));
}


WayPoint::~WayPoint()
{
//
	cout << "§WayPoint::~WayPoint() "; _ctrsgn(typeid(*this), true);
}


/// __________________________________________________
/// Convert waypoint format
unique_ptr<const WayPoint> WayPoint::convert(const CoordType type) const
{
//	cout << "@WayPoint::convert(const CoordType type) " << typeid(*this).name()
//       << " to type "<< coordtype_to_int(type) + 1 << endl;
	return factory<const WayPoint>(*this, type);
}


/// __________________________________________________
/// Get const reference to cbp_lat or cbp_lon
inline const CoordBase& WayPoint::get_cbp(bool latlon) const
{
//	cout << "@CoordBase::get_cbp(bool)\n";
	return *(latlon ? cbp_lat : cbp_lon).get() ;
}


/// __________________________________________________
/// Formatted strings for printing
vector<string> WayPoint::format() const
{
//	cout << "@WayPoint::format()\n";
	vector<string> sv_lat{ cbp_lat->fmt_fctr_tmpl() };
	vector<string> sv_lon{ cbp_lon->fmt_fctr_tmpl() };
	vector<string> out(sv_lat.size());
	transform(
		sv_lat.begin(), sv_lat.end(), sv_lon.begin(), out.begin(),
		[](string& latstr, string& lonstr) { return latstr + "  " + lonstr; }
	);
	if (cbp_lat->get_names().size())
		transform(
			out.begin(), out.end(), cbp_lat->get_names().begin(), out.begin(),
			[](string& lls, const string& name) { return lls + "  " + name; }
		);
	return out;
}


/// __________________________________________________
/// Print WayPoint
void WayPoint::print(ostream& stream) const
{
//	cout << "@WayPoint::print() " << typeid(*this).name() << endl;
	const int i { coordtype_to_int(cbp_lat->getfmt()) };
	vector<int> spacing { 5, 7, 8, 11, 13, 14, 2, 2, 2 };
	stream << " Latitude" << string(spacing[i], ' ') << "Longitude\n"
		   << string(1, ' ') << string(spacing[i + 3], '_')
		   << string(spacing[i + 6], ' ') << string(spacing[i + 3] + 1, '_') << endl;
	vector<string> sv(format());
	for_each(sv.begin(), sv.end(), [&stream](const string& s) { stream << s << "\n"; });
}


/// __________________________________________________
/// Validate WayPoint
void WayPoint::validate(bool warn = true) const
{
//	cout << "@WayPoint::validate(bool)\n";
	cbp_lat->validate_tmpl(warn);
	cbp_lon->validate_tmpl(warn);
}


/// __________________________________________________
/// WayPoint validity
const vector<bool>& WayPoint::get_valid(bool latlon) const
{
//	cout << "@WayPoint::get_valid(bool)\n";
	return latlon ? validlat : validlon;
}


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


/// __________________________________________________
/// Create unique_ptr<const WayPoint> to new WayPoint object
template<class T>
unique_ptr<const WayPoint> newconstWaypoint(const T& t)
{
//	cout << "@newconstWaypoint(const T&) fmt " << get_fmt_attribute(t) << endl;
	return factory<const WayPoint>(
		newconstCoordBase(as<NumericVector>(t[1]), get_coordtype(t)),
		newconstCoordBase(as<NumericVector>(t[2]), get_coordtype(t))
	);
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
	unique_ptr<const CoordBase> cb{ newconstCoordBase(nv, get_coordtype(nv)) };
	cb->validate_tmpl();
	const_cast<NumericVector&>(nv).attr("valid") = cb->get_valid();
	return cb->get_valid();
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
	unique_ptr<const CoordBase> cb1{ newconstCoordBase(nv, inheritscoords ? oldtype : newtype) };
	if (inheritscoords) {
		unique_ptr<const CoordBase> cb2{ cb1->convert(newtype) };
		cb1.swap(cb2);
		copy((cb1->get_nv()).begin(), (cb1->get_nv()).end(), nv.begin());
	} else {
		nv.attr("class") = "coords";
	}
	nv.attr("fmt") = fmt;
	cb1->validate_tmpl();
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
	Rcout << *newconstCoordBase(nv, get_coordtype(nv)) << endl;
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

/// __________________________________________________
/// __________________________________________________