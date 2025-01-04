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

class CoordBase;
ostream& operator<<(ostream&, const CoordBase&);

template<class T> 
string format_coord(const T&, double);
class DecDeg;

class DegMin;
template<> 
string format_coord<DegMin>(const DegMin&, double);

class DegMinSec;
template<> 
string format_coord<DegMinSec>(const DegMinSec&, double);

class FormatBase;
class Format_DD;
class Format_DM;
class Format_DMS;
class FormatLL_DM_S;
class FormatLL_DD;

template<class T>
unique_ptr<const CoordBase> newconstCoordBase(const T&, const CoordType);

class WayPoint;
template <class T>
unique_ptr<const WayPoint> newconstWayPoint(const T&);
ostream& operator<<(ostream&, const WayPoint&);

template<class T>
inline int get_fmt_attribute(const T&);
template<class T>
inline void checkinherits(T&, const char *);
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
DataFrame waypoints_replace(DataFrame &df, int value);
DataFrame printwaypoint(DataFrame&);
const DataFrame validatewaypoint(DataFrame&);

/// __________________________________________________
/// __________________________________________________
/// Development and Debugging functions

/// Report object construction and destruction
void _ctrsgn(const type_info &obj, bool destruct = false)
{
///§	cout << (destruct ? "Destroying " : "Constructing ") << obj.name() << endl;
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
//	cout << "@factory(A1& a1)\n";
	return unique_ptr<T>(new T(std::forward<A1>(a1)));
}

template <class T, class A1, class A2>
inline unique_ptr<T>
factory(A1&& a1, A2&& a2)   // two argument version
{
//	cout << "@factory(A1& a1, A2& a2)\n";
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
inline const CoordType get_coordtype(const T &t)
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
/// Coordinate abstract base class

class CoordBase {
	protected:
		vector<double> nv;
		const vector<bool> valid { false };
		const vector<bool> latlon;
		const vector<string> names;
		const bool llgt1 = false;
		bool validator(double, bool) const;
		bool all_valid() const;
		bool waypoint = false;

	public:
		CoordBase(const NumericVector&);
		explicit CoordBase(const CoordBase&);
		CoordBase(const vector<double>, const vector<bool>&, const vector<string>&);
		CoordBase& operator=(const CoordBase&) = delete;

		virtual ~CoordBase() = 0;
		virtual const CoordType getfmt() const = 0;
		virtual int get_deg(double) const = 0;
		virtual double get_decdeg(double) const = 0;
		virtual int get_min(double) const = 0;
		virtual double get_decmin(double) const = 0;
		virtual double get_sec(double) const = 0;
		virtual vector<string> fmt_fctr_tmpl() const = 0;

		const vector<double> &get_nv() const;
		unique_ptr<const CoordBase> convert(const CoordType) const;
		void validate(bool) const;
		const vector<bool> &get_valid() const;
		const vector<string> &get_names() const;
		void warn_invalid() const;
		void set_waypoint() const;
		template <typename FunctObj, typename FunctObj2>
		vector<string> format() const;
		void print(ostream&) const;

		friend class DecDeg;
		friend class DegMin;
		friend class DegMinSec;
		friend class FormatLL_DM_S;
		friend class FormatLL_DD;
		friend ostream& operator<<(ostream&, const CoordBase&);
};

CoordBase::CoordBase(const NumericVector& nv) :
	CoordBase(
		as<vector<double>>(nv),
		nv.hasAttribute("latlon") ? as<vector<bool>>(nv.attr("latlon")) : vector<bool>(),
		nv.hasAttribute("names") ? as<vector<string>>(nv.attr("names")) : vector<string>()
	)
{
///§	cout << "@CoordBase::CoordBase(const NumericVector&) "; _ctrsgn(typeid(*this));
}


CoordBase::CoordBase(const CoordBase &c) :
	CoordBase(vector<double>(c.nv.size()), vector<bool>{ c.latlon }, vector<string>{ c.names })
{
///§	cout << "@CoordBase::CoordBase(const CoordBase&) "; _ctrsgn(typeid(*this));
}


CoordBase::CoordBase(const vector<double> n, const vector<bool>& ll, const vector<string>& _names) :
	nv(std::move(n)), latlon{ ll }, names{ std::move(_names) }, llgt1(latlon.size() > 1)
{
///§	cout << "@CoordBase::CoordBase(const vector<double>, const LogicalVector&, const vector<string>&) "; _ctrsgn(typeid(*this));
}


CoordBase::~CoordBase()
{
///§	cout << "@CoordBase::~CoordBase() "; _ctrsgn(typeid(*this), true);
}


/// __________________________________________________
/// Get const reference to nv
inline const vector<double> &CoordBase::get_nv() const
{
//	cout << "@CoordBase::get_nv()\n";
	return nv;
}


/// __________________________________________________
/// Get const reference to valid
inline const vector<bool> &CoordBase::get_valid() const
{
	return valid;
}


/// __________________________________________________
/// Get const reference to names
inline const vector<string> &CoordBase::get_names() const
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
/// Validate coord value
inline bool CoordBase::validator(double n, bool lat = false) const
{
//	cout << "@CoordBase::validator() " << " n: " << setw(9) << setfill(' ') << n << " lat: ";
//  if (latlon.size()) cout << boolalpha << lat << endl; else cout << "NA\n";

	return !((abs(get_decdeg(n)) > (latlon.size() && lat ? 90 : 180)) ||
		(abs(get_decmin(n)) >= 60) ||
		(abs(get_sec(n)) >= 60));
}


/// __________________________________________________
/// Validate coords vector
void CoordBase::validate(bool warn = true) const
{
//	cout << "@CoordBase::validate() " << typeid(*this).name() << " latlon " << LogicalVector(wrap(latlon)) << endl;
	vector<bool>& non_const_valid { const_cast<vector<bool>&>(valid) };
	non_const_valid.assign(nv.size(), {false});
	if (latlon.size()) {
		vector<bool>::const_iterator ll_it(latlon.begin());
		transform(nv.begin(), nv.end(), non_const_valid.begin(),
		 [this, &ll_it](double n) { return validator(n, llgt1 ? *ll_it++ : *ll_it); });
	} else
		transform(nv.begin(), nv.end(), non_const_valid.begin(),
		 [this](double n) { return validator(n); });
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
	bool &wpt = const_cast<bool&>(waypoint);
	wpt = true;
}


/// __________________________________________________
/// Formatted coordinate strings for printing
template <typename FunctObj, typename FunctObj2>
vector<string> CoordBase::format() const
{
//	cout << "CoordBase::format<typename FunctObj, FunctObj2>()\n";
	vector<string> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), FunctObj(*this));
	transform(out.begin(), out.end(), nv.begin(), out.begin(), FunctObj2(*this));
	return out;
}


/// __________________________________________________
/// Print coords vector
void CoordBase::print(ostream& stream) const
{
//	cout << "@CoordBase::print() type " << typeid(*this).name() << endl;
	vector<string> sv(fmt_fctr_tmpl()); 
	if (names.size()) {
		vector<string>::const_iterator nm_it(names.begin());
		for_each(sv.begin(), sv.end(),
			[&stream, &nm_it](const string &s) { stream << s << " " << *nm_it++ << "\n"; });
	} else
		for_each(sv.begin(), sv.end(), [&stream](const string &s) { stream << s << "\n"; });
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
		int get_deg(double) const;
		double get_decdeg(double x) const;
		int get_min(double x) const;
		double get_decmin(double x) const;
		double get_sec(double x) const;
		vector<string> fmt_fctr_tmpl() const;
};


DecDeg::DecDeg(const NumericVector &nv) : CoordBase(nv)
{
///§	cout << "@DecDeg::DecDeg(NumericVector&) "; _ctrsgn(typeid(*this));
}


DecDeg::DecDeg(const CoordBase &c) : CoordBase(c)
{
///§	cout << "@DecDeg::DecDeg(const CoordBase&) "; _ctrsgn(typeid(*this));
	transform(c.nv.begin(), c.nv.end(), nv.begin(), [&c](double n) { return c.get_decdeg(n); });
}


DecDeg::~DecDeg()
{
///§	cout << "@DecDeg::~DecDeg() "; _ctrsgn(typeid(*this), true);
}

inline int DecDeg::get_deg(double x) const
{
//	cout << "DecDeg.get_deg()\n";
	return int(x);
}

inline double DecDeg::get_decdeg(double x) const
{
//	cout << "DecDeg.get_decdeg()\n";
	return x;
}

inline int DecDeg::get_min(double x) const
{
//	cout << "DecDeg.get_min()\n";
	return (int(x * 1e6) % int(1e6)) * 6e-5;
}

inline double DecDeg::get_decmin(double x) const
{
//	cout << "DecDeg.get_decmin()\n";
	return polish(mod1by60(x));
}

inline double DecDeg::get_sec(double x) const
{
//	cout << "DecDeg.get_sec()\n";
	return mod1by60(get_decmin(x));
}


/// __________________________________________________
/// Instantiate functor template for formatting decimal degrees
inline vector<string> DecDeg::fmt_fctr_tmpl() const
{
//	cout << "DecDeg::fmt_fctr_tmpl()\n";
	return format<Format_DD, FormatLL_DD>();
}


/// __________________________________________________
/// Degrees and minutes derived class
class DegMin : public CoordBase {
	public:
		DegMin(const NumericVector&);
		DegMin(const CoordBase&);
		~DegMin();

		const CoordType getfmt() const { return CoordType::degmin; }
		int get_deg(double) const;
		double get_decdeg(double x) const;
		int get_min(double x) const;
		double get_decmin(double x) const;
		double get_sec(double x) const;
		vector<string> fmt_fctr_tmpl() const;
};


DegMin::DegMin(const NumericVector &nv) : CoordBase(nv)
{
///§	cout << "@DegMin::DegMin(NumericVector&) "; _ctrsgn(typeid(*this));
}

DegMin::DegMin(const CoordBase &c) : CoordBase(c)
{
///§	cout << "@DegMin::DegMin(const CoordBase&) "; _ctrsgn(typeid(*this));
	transform(c.nv.begin(), c.nv.end(), nv.begin(),
		[&c](double n) { return c.get_deg(n) * 1e2 + c.get_decmin(n); });
}

DegMin::~DegMin()
{
///§	cout << "@DegMin::~DegMin() "; _ctrsgn(typeid(*this), true);
}

inline int DegMin::get_deg(double x) const
{
//	cout << "DegMin.get_deg()\n";
	return int(x / 1e2);
}

inline double DegMin::get_decdeg(double x) const
{
//	cout << "DegMin.get_decdeg()\n";
	return int(x / 1e2) + mod1e2(x) / 60;
}

inline int DegMin::get_min(double x) const
{
//	cout << "DegMin.get_min()\n";
	return int(x) % int(1e2);
}

inline double DegMin::get_decmin(double x) const
{
//	cout << "DegMin.get_decmin()\n";
	return polish(mod1e2(x));
}

inline double DegMin::get_sec(double x) const
{
//	cout << "DegMin.get_sec()\n";
	return mod1by60(get_decmin(x));
}


/// __________________________________________________
/// Instantiate functor template for formatting degrees and minutes
inline vector<string> DegMin::fmt_fctr_tmpl() const
{
//	cout << "DegMin::fmt_fctr_tmpl()\n";
	return format<Format_DM, FormatLL_DM_S>();
}


/// __________________________________________________
/// Degrees minutes and seconds derived class
class DegMinSec : public CoordBase {
	public:
		DegMinSec(const NumericVector&);
		DegMinSec(const CoordBase&);
		~DegMinSec();

		const CoordType getfmt() const { return CoordType::degminsec; }
		int get_deg(double) const;
		double get_decdeg(double x) const;
		int get_min(double x) const;
		double get_decmin(double x) const;
		double get_sec(double x) const;
		vector<string> fmt_fctr_tmpl() const;
};


DegMinSec::DegMinSec(const NumericVector &nv) : CoordBase(nv)
{
///§	cout << "@DegMinSec::DegMinSec(NumericVector&) "; _ctrsgn(typeid(*this));
}

DegMinSec::DegMinSec(const CoordBase &c) : CoordBase(c)
{
///§	cout << "@DegMinSec::DegMinSec(const CoordBase&) "; _ctrsgn(typeid(*this));
	transform(c.nv.begin(), c.nv.end(), nv.begin(),
		[&c](double n) { return c.get_deg(n) * 1e4 + c.get_min(n) * 1e2 + c.get_sec(n); });
}

DegMinSec::~DegMinSec()
{
///§	cout << "@DegMinSec::~DegMinSec() "; _ctrsgn(typeid(*this), true);
}

inline int DegMinSec::get_deg(double x) const
{
//	cout << "DegMinSec.get_deg()\n";
	return int(x / 1e4);
}

inline double DegMinSec::get_decdeg(double x) const
{
//	cout << "DegMinSec.get_decdeg()\n";
	return int(x / 1e4) + (double)int(fmod(x, 1e4) / 1e2) / 60 + mod1e2(x) / 3600;
}

inline int DegMinSec::get_min(double x) const
{
//	cout << "DegMinSec.get_min()\n";
	return (int(x) % int(1e4)) / 1e2;
}

inline double DegMinSec::get_decmin(double x) const
{
//	cout << "DegMinSec.get_decmin()\n";
	return int(fmod(x, 1e4) / 1e2) + mod1e2(x) / 60;
}

inline double DegMinSec::get_sec(double x) const
{
//	cout << "DegMinSec.get_sec()\n";
	return mod1e2(x);
}


/// __________________________________________________
/// Instantiate functor template for formatting degrees, minutes and seconds
inline vector<string> DegMinSec::fmt_fctr_tmpl() const
{
//	cout << "DegMinSec::fmt_fctr_tmpl()\n";
	return format<Format_DMS, FormatLL_DM_S>();
}


/// __________________________________________________
/// __________________________________________________
/// Formatting functors

/// __________________________________________________
/// Format coords vector functor base class
class FormatBase {
	protected:
		const CoordBase& cb; 
		ostringstream outstrstr;
	public:
		FormatBase(const CoordBase& _cb) : cb(_cb)
		{
//			cout << "@FormatBase(const CoordBase& _cb) "; _ctrsgn(typeid(*this));
		}
		virtual ~FormatBase() = 0;
};

inline FormatBase::~FormatBase() {}


/// __________________________________________________
/// Format coords vector functor for decimal degrees
class Format_DD : public FormatBase {
public:
	Format_DD(const CoordBase& cb) : FormatBase(cb)
	{
//		cout << "@Format_DD(const CoordBase& cb) "; _ctrsgn(typeid(*this));
	}
	string operator()(double n)
	{
//		cout << "@Format_DD::operator()\n";
		outstrstr.str("");
		outstrstr << setw(11) << setfill(' ') << fixed << setprecision(6) << cb.get_decdeg(n) << "\u00B0";
		return outstrstr.str();
	}
};


/// __________________________________________________
/// Format coords vector functor for degrees and minutes
class Format_DM : public FormatBase{
public:
	Format_DM(const CoordBase& cb) : FormatBase(cb)
	{
//		cout << "@Format_DM(const CoordBase& cb) "; _ctrsgn(typeid(*this));
	}
	string operator()(double n)
	{
//		cout << "@Format_DM::operator()\n";
		outstrstr.str("");
		outstrstr << setw(3) << setfill(' ') << abs(cb.get_deg(n)) << "\u00B0"
				  << setw(7) << setfill('0') << fixed << setprecision(4) << abs(cb.get_decmin(n)) << "'";
		return outstrstr.str();
	}
};


/// __________________________________________________
/// Format coords vector functor for degrees, minutes and seconds 
class Format_DMS : public FormatBase{
public:
	Format_DMS(const CoordBase& cb) : FormatBase(cb)
	{
//		cout << "@Format_DMS(const CoordBase& cb) "; _ctrsgn(typeid(*this));
	}
	string operator()(double n)
	{
//		cout << "@Format_DMS::operator()\n";
		outstrstr.str("");
		outstrstr << setw(3) << setfill(' ') << abs(cb.get_deg(n)) << "\u00B0"
				  << setw(2) << setfill('0') << abs(cb.get_min(n)) << "'"
				  << setw(5) << fixed << setprecision(2) << abs(cb.get_sec(n)) << "\"";
		return outstrstr.str();
	}
};


/// __________________________________________________
/// Format functor for latitude and longitude strings for decimal degrees
class FormatLL_DD : public FormatBase {
	vector<bool>::const_iterator ll_it;
public:
	FormatLL_DD(const CoordBase& cb) : FormatBase(cb), ll_it(cb.latlon.begin()) {}
	string operator()(string ostr, double n)
	{
//		cout << "@FormatLL_DD::operator() cb.waypoint " << boolalpha << cb.waypoint << endl;
		if (cb.latlon.size() && !cb.waypoint)
			return ostr += ((cb.llgt1 ? *ll_it++ : *ll_it) ? " lat" : " lon");
		else
			return ostr;
	}
};


/// __________________________________________________
/// Format functor for latitude and longitude strings for degrees, minutes (and seconds)
class FormatLL_DM_S : public FormatBase {
	vector<bool>::const_iterator ll_it;
public:
	FormatLL_DM_S(const CoordBase& cb) : FormatBase(cb), ll_it(cb.latlon.begin()) {}
	string operator()(string ostr, double n)
	{
//		cout << "@FormatLL_DM_S::operator()\n";
		return ostr += cb.latlon.size() ? cardpoint(cb.get_decmin(n) < 0, cb.llgt1 ? *ll_it++ : *ll_it) : cardi_b(cb.get_decmin(n) < 0);
	}
};


/// __________________________________________________
/// __________________________________________________
/// Create unique_ptr<CoordBase> to new DecDeg, DegMin or DegMinSec object
template<class T>
unique_ptr<const CoordBase> newconstCoordBase(const T &t, const CoordType type)
{
//	cout << "@newconstCoordBase<T>(const T&, const CoordType) of type " << coordtype_to_int(type) + 1 << endl;

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
		const vector<bool> &validlat;
		const vector<bool> &validlon;
	public:
		explicit WayPoint(unique_ptr<const CoordBase>, unique_ptr<const CoordBase>);
		explicit WayPoint(const WayPoint&, CoordType);
		~WayPoint();

		const CoordBase &get_cbp(bool) const;
		unique_ptr<const WayPoint> convert(const CoordType) const;
		void validate(bool) const;
		const vector<bool> &get_valid(bool) const;
		void warn_invalid() const;
		void print(ostream& stream) const;
		vector<string> format() const;
		friend ostream& operator<<(ostream&, const WayPoint&);
};


WayPoint::WayPoint(unique_ptr<const CoordBase> _cbp_lat, unique_ptr<const CoordBase> _cbp_lon) :
	cbp_lat{std::move(_cbp_lat)}, cbp_lon{std::move(_cbp_lon)},
	validlat(cbp_lat->get_valid()), validlon(cbp_lon->get_valid())
{
///§	cout << "@WayPoint(unique_ptr<const CoordBase>, unique_ptr<const CoordBase>) "; _ctrsgn(typeid(*this));
	cbp_lat->set_waypoint();
	cbp_lon->set_waypoint();
}

WayPoint::WayPoint(const WayPoint &wp, CoordType type) :
	WayPoint{ wp.get_cbp(true).convert(type), wp.get_cbp(false).convert(type) }
{
///§	cout << "@WayPoint(const WayPoint&) "; _ctrsgn(typeid(*this));
}


WayPoint::~WayPoint()
{
///§	cout << "@WayPoint::~WayPoint() "; _ctrsgn(typeid(*this), true);
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
inline const CoordBase &WayPoint::get_cbp(bool latlon) const
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
		[](string &latstr, string &lonstr) { return latstr + "  " + lonstr; }
	);
	if (cbp_lat->get_names().size())
		transform(
			out.begin(), out.end(), cbp_lat->get_names().begin(), out.begin(),
			[](string &lls, const string &name) { return lls + "  " + name; }
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
	for_each(sv.begin(), sv.end(), [&stream](const string &s) { stream << s << "\n"; });
}


/// __________________________________________________
/// Validate WayPoint
void WayPoint::validate(bool warn = true) const
{
//	cout << "@WayPoint::validate(bool)\n";
	cbp_lat->validate(warn);
	cbp_lon->validate(warn);
}


/// __________________________________________________
/// WayPoint validity
const vector<bool> &WayPoint::get_valid(bool latlon) const
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
unique_ptr<const WayPoint> newconstWaypoint(const T &t)
{
//	cout << "@newconstWaypoint(const T&) fmt " << get_fmt_attribute(t) << endl;
	return factory<const WayPoint>(
		newconstCoordBase(as<NumericVector>(t[1]), get_coordtype(t)),
		newconstCoordBase(as<NumericVector>(t[2]), get_coordtype(t))
	);
}


/// __________________________________________________
/// Output WayPoint to ostream
ostream& operator<<(ostream& stream, const WayPoint &wp)
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
inline int get_fmt_attribute(const T &t)
{
//	cout << "@get_fmt_attribute<T>(const T&) " << as<int>(t.attr("fmt")) << endl;
	return as<int>(t.attr("fmt"));
}


/// __________________________________________________
/// Does object inherit given class?
template<class T>
inline void checkinherits(T &t, const char *classname)
{
//	cout << "checkinherits(T &t, const char *classname) t " << typeid(t).name() << " classname " << classname << endl;
	if (!t.inherits(classname)) stop("Argument must be a \"%s\" object", classname);
}


/// __________________________________________________
/// set attributes for vector column within object
template<class T, class V>
void colattrset(const T &t, int col, const char *attrib, V &&val)
{
//	cout << "@colattrset(const T&, int, const char*, V&&) attrib " << attrib << ", col " << col << endl;
	as<NumericVector>(t[col]).attr(attrib) = std::forward<V>(val);
}


/// __________________________________________________
/// __________________________________________________
/// Validation functions

/// __________________________________________________
/// Has R coords object been validated?
inline bool validcoord(NumericVector &nv)
{
//	cout << "@validcoord(NumericVector&)\n";
	LogicalVector lv { as<LogicalVector>(nv.attr("valid")) };
	return 1 == lv.size() && lv[0];
}


/// __________________________________________________
/// Return "valid" attribute or empty LogicalVector    !!!!!!! Generalise with template & specialisation !!!!!!!
inline LogicalVector get_valid(const NumericVector &nv)
{
//	cout << "@get_valid(const NumericVector&) has attr \"valid\" " << boolalpha << nv.hasAttribute("valid") << endl;
	return (nv.hasAttribute("valid") ? LogicalVector(nv.attr("valid")) : LogicalVector());
}


/// __________________________________________________
/// Check first two columns have attribute "valid" all true
template<class T>
bool check_valid(const T &t)
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
bool check_valid<NumericVector>(const NumericVector &nv)
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
vector<bool> validatecoord(const NumericVector &nv)
{
//	cout << "@validatecoord()\n";
	unique_ptr<const CoordBase> cb{ newconstCoordBase(nv, get_coordtype(nv)) };
	cb->validate();
	const_cast<NumericVector&>(nv).attr("valid") = cb->get_valid();
	return cb->get_valid();
}


/// __________________________________________________
/// __________________________________________________
/// Exported functions


/// __________________________________________________
/// Set R vector object class to coords and return,
/// or convert format of R coords object and return  !!!!!!! Template & Specialise [with waypoints()] !!!!!!!
// [[Rcpp::export]]
NumericVector coords(NumericVector &nv, int fmt = 1)
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
	cb1->validate();
	nv.attr("valid") = cb1->get_valid();
	return nv;
}


/// __________________________________________________
/// coords() as replacement function
// [[Rcpp::export(name = "`coords<-`")]]
NumericVector coords_replace(NumericVector &nv, int value)
{
//	cout << "——Rcpp::export——`coords_replace()<-`\n";
	return coords(nv, value);
}


/// __________________________________________________
/// Set latlon attribute on "coords" NumericVector and revalidate
// [[Rcpp::export(name = "`latlon<-`")]]
NumericVector latlon(NumericVector &nv, LogicalVector &value)
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
NumericVector printcoord(NumericVector &nv)
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
vector<bool> Rvalidatecoord(NumericVector &nv)
{
//	cout << "——Rcpp::export——Rvalidatecoord()\n";
	checkinherits(nv, "coords");
	return validatecoord(nv);
}


/// __________________________________________________
/// Format coords vector - S3 method format.coords()
// [[Rcpp::export(name = "format.coords")]]
vector<string> formatcoord(NumericVector &nv)
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
vector<int> get_deg(NumericVector &nv)
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
vector<double> get_decdeg(NumericVector &nv)
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
vector<int> get_min(NumericVector &nv)
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
vector<double> get_decmin(NumericVector &nv)
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
vector<double> get_sec(NumericVector &nv)
{
//	cout << "——Rcpp::export——get_sec()\n";
	checkinherits(nv, "coords");
	unique_ptr<const CoordBase> c{newconstCoordBase(nv, get_coordtype(nv))};
	vector<double> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), [&c](double n) { return c->get_sec(n); });
	return out;
}


/// __________________________________________________
/// Add "waypoints" to R data.frame object class and validate. !!!!!!! Template & Specialise !!!!!!!
// [[Rcpp::export]]
DataFrame waypoints(DataFrame &df, int fmt = 1)
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
DataFrame waypoints_replace(DataFrame &df, int value)
{
//	cout << "——Rcpp::export——`waypoints_replace()<-`\n";
	return waypoints(df, value);
}


/// __________________________________________________
/// Print waypoints vector - S3 method print.waypoints()      /////// "invisible" not working ///////
// [[Rcpp::export(name = "print.waypoints", invisible = true)]]
DataFrame printwaypoint(DataFrame &df)
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
const DataFrame validatewaypoint(DataFrame &df)
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