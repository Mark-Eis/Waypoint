/// __________________________________________________
/// CoordBase.h
/// __________________________________________________

#ifndef COORDBASE_H_
#define COORDBASE_H_

#define FMT_HEADER_ONLY
#include "fmt/base.h"		// …fmt/*.h copied to …/R/Packages/Waypoint/src.
#include <concepts>

/// __________________________________________________
/// __________________________________________________
/// Development and debugging

#define DEBUG 0

#if DEBUG > 0

void _ctrsgn(const std::type_info&, bool = true);

/// __________________________________________________
/// Check address for unnecessary copying
inline void* address(const auto& t)
{
	return (void*)&t;
}

/// __________________________________________________
/// Format strings for debugging code
constexpr auto padstr { "— — — — — — — — — "sv };
constexpr auto exportstr { "——Rcpp::export——"sv };

#endif	// if DEBUG > 0

const string demangle(const std::type_info&);


/// __________________________________________________
/// __________________________________________________
/// Type Traits

/// NumericVector
template <typename T>
struct isNumericVector : public std::false_type {};

template <>
struct isNumericVector<NumericVector> : public std::true_type {};

template<typename T>
constexpr bool isNumericVector_v = isNumericVector<T>::value;

/// DataFrame
template <typename T>
struct isDataFrame : public std::false_type {};

template <>
struct isDataFrame<DataFrame> : public std::true_type {};

template<typename T>
constexpr bool isDataFrame_v = isDataFrame<T>::value;

/// __________________________________________________
/// Concepts

/// Concept —— NumericVector
template<typename T>
concept Is_NumericVector = isNumericVector_v<T>;

/// Concept —— DataFrame
template<typename T>
concept Is_DataFrame =
	requires(T t, const string& s, const char *c) {
		{ t.attr(s) } -> std::same_as<Rcpp::AttributeProxyPolicy<Rcpp::Vector<19>>::AttributeProxy>;
		{ t.attributeNames() } -> std::same_as<vector<string>>;
		{ t.hasAttribute(s) } -> std::same_as<bool>;
		{ t.inherits(c) } -> std::same_as<bool>;
		{ t.length() } -> std::integral;
		{ t.names() } -> std::same_as<Rcpp::NamesProxyPolicy<Rcpp::Vector<19>>::NamesProxy>;
		{ t.nrows() } -> std::integral;
	};

/// Concept —— Either NumericVector or DataFrame
template<typename T>
concept NumVec_or_DataFrame =
	Is_NumericVector<T> || Is_DataFrame<T>;


/// __________________________________________________
/// __________________________________________________
/// Class and Function declarations

/// __________________________________________________
/// __________________________________________________
/// DVecType and SVecType

/// __________________________________________________
/// VecTypeBase
template<typename T>
struct VecTypeBase : public vector<T> {
	VecTypeBase(const VecTypeBase&) = delete;											// copy constructor
	VecTypeBase(const vector<T>& vt) : vector<T>{ vt }									// copy constructor
	{
#if DEBUG > 0
		_ctrsgn(typeid(*this)); fmt::print("\t(const vector<T>&)\n"); 
#endif
	}
	VecTypeBase(const NumericVector vt) : vector<T>{ as<vector<double>>(vt) }			// copy constructor
	{
#if DEBUG > 0
		_ctrsgn(typeid(*this)); fmt::print("\t(const NumericVector&)\n");
#endif
	}

	VecTypeBase& operator=(const VecTypeBase&) = delete;									// copy assignment
	VecTypeBase& operator=(const vector<T>& vt)											// copy assignment
	{
#if DEBUG > 0
		fmt::print("@VecTypeBase& operator=(const vector<T>& vt)\n");
#endif
		vector<T>::operator= (vt);
		return *this;
	}
	VecTypeBase& operator=(const NumericVector) = delete;								// copy assignment - not defaultable

#if DEBUG == 0
	VecTypeBase(VecTypeBase&&) = default;												// move constructor
#else 
	VecTypeBase(VecTypeBase&& dv) : vector<T>{ std::move(static_cast<vector<T>&&>(dv)) }
	{
		_ctrsgn(typeid(*this)); fmt::print("\t(VecTypeBase&&)\n");
	}
#endif

	VecTypeBase(vector<T>&& vt) : vector<T>{ std::move(vt) }								// move constructor
	{
#if DEBUG > 0
		_ctrsgn(typeid(*this)); fmt::print("\t(vector<T>&&)\n");
#endif
	}
	VecTypeBase(NumericVector&& vt) = delete;											// move constructor - not defaultable

#if DEBUG == 0
VecTypeBase& operator=(VecTypeBase&&) = default;											// move assignment
#else 
	VecTypeBase& operator=(VecTypeBase&& dv)
	{
		_ctrsgn(typeid(*this)); fmt::print("\tVecTypeBase& operator=(VecTypeBase&&)\n");
		vector<T>::operator=(std::move(static_cast<vector<T>&&>(dv)));
		return *this;
	}
#endif

	VecTypeBase& operator=(vector<T>&& vt)											// move assignment
	{
#if DEBUG > 0
		fmt::print("@VecTypeBase& operator=(vector<T>&& vt)\n");
#endif
		vector<T>::operator=(std::move(vt));
		return *this;
	}
	VecTypeBase& operator=(NumericVector&& vt)	 = delete;								// move assignment - not defaultable

	virtual ~VecTypeBase() = 0;
};

template<typename T>
VecTypeBase<T>::~VecTypeBase()
{
#if DEBUG > 0
	_ctrsgn(typeid(*this), false);
#endif
}

/// __________________________________________________
/// DecDegVec
template<typename T>
struct DecDegVec final : public VecTypeBase<T> {
	using VecTypeBase<T>::VecTypeBase;
};

/// __________________________________________________
/// DegMinVec
template<typename T>
struct DegMinVec final : public VecTypeBase<T> {
	using VecTypeBase<T>::VecTypeBase;
};

/// __________________________________________________
/// DegMinSecVec
template<typename T>
struct DegMinSecVec final : public VecTypeBase<T> {
	using VecTypeBase<T>::VecTypeBase;
};


/// __________________________________________________
/// Template aliases
using DecDegVecDouble = DecDegVec<double>;
using DegMinVecDouble = DegMinVec<double>;
using DegMinSecVecDouble = DegMinSecVec<double>;

/// __________________________________________________
/// Type Traits

/// DecDegVecDouble
template <typename T>
struct isDecDegVecDouble : public std::false_type {};

template <>
struct isDecDegVecDouble<DecDegVecDouble> : public std::true_type {};

template<typename T>
constexpr bool isDecDegVecDouble_v = isDecDegVecDouble<T>::value;

/// DegMinVecDouble
template <typename T>
struct isDegMinVecDouble : public std::false_type {};

template <>
struct isDegMinVecDouble<DegMinVecDouble> : public std::true_type {};

template<typename T>
constexpr bool isDegMinVecDouble_v = isDegMinVecDouble<T>::value;

/// DegMinSecVecDouble
template <typename T>
struct isDegMinSecVecDouble : public std::false_type {};

template <>
struct isDegMinSecVecDouble<DegMinSecVecDouble> : public std::true_type {};

template<typename T>
constexpr bool isDegMinSecVecDouble_v = isDegMinSecVecDouble<T>::value;

/// __________________________________________________
/// Concept —— DVecType
template <typename T>
concept DVecType = 
	isDecDegVecDouble_v<T> ||
	isDegMinVecDouble_v<T> ||
	isDegMinSecVecDouble_v<T>;

/// __________________________________________________
/// Template aliases
using DecDegVecString = DecDegVec<string>;
using DegMinVecString = DegMinVec<string>;
using DegMinSecVecString = DegMinSecVec<string>;

/// __________________________________________________
/// Type Traits

/// DecDegVecString
template <typename T>
struct isDecDegVecString : public std::false_type {};

template <>
struct isDecDegVecString<DecDegVecString> : public std::true_type {};

template<typename T>
constexpr bool isDecDegVecString_v = isDecDegVecString<T>::value;

/// DegMinVecString
template <typename T>
struct isDegMinVecString : public std::false_type {};

template <>
struct isDegMinVecString<DegMinVecString> : public std::true_type {};

template<typename T>
constexpr bool isDegMinVecString_v = isDegMinVecString<T>::value;

/// DegMinSecVecString
template <typename T>
struct isDegMinSecVecString : public std::false_type {};

template <>
struct isDegMinSecVecString<DegMinSecVecString> : public std::true_type {};

template<typename T>
constexpr bool isDegMinSecVecString_v = isDegMinSecVecString<T>::value;

/// __________________________________________________
/// Concept —— SVecType
template <typename T>
concept SVecType = 
	isDecDegVecString_v<T> ||
	isDegMinVecString_v<T> ||
	isDegMinSecVecString_v<T>;


/// __________________________________________________
/// __________________________________________________
/// Formula simplification
inline double mod1by60(double);
inline double mod1e2(double);
inline double round2(double, int);
inline double polish(double);

/// __________________________________________________
/// __________________________________________________
/// Utility
template<NumVec_or_DataFrame T, typename U> 
inline vector<U> get_vec_attr(const T&, const string);
inline int get_fmt_attribute(const NumVec_or_DataFrame auto&);
template<NumVec_or_DataFrame T>
int check_logical_attr(T t, const string attrname);
inline void checkinherits(const NumVec_or_DataFrame auto&, const string);
inline bool is_item_in_df(const DataFrame, int);
inline void stdlenstr(vector<string>&);
inline void concat_vecstr_elmnts(const vector<string>&, vector<string>&, const string = " ");
inline void concat_vecstr_elmnts(const vector<int>&, vector<string>&, const string = " ");
inline bool prefixwithnames(vector<string>&, RObject&);
inline string str_tolower(string);
int name_pos_in_df(const DataFrame, const string);
RObject getnames(const DataFrame);

/// __________________________________________________
/// __________________________________________________
/// CoordType enum
enum class CoordType : char { decdeg, degmin, degminsec };

template<>
struct fmt::formatter<CoordType>: formatter<string_view>
{
	auto format(CoordType, format_context&) const
		-> format_context::iterator;
};

/// __________________________________________________
/// CoordType access functions
inline const CoordType get_coordtype(int);
inline const CoordType get_coordtype(const NumVec_or_DataFrame auto&);
inline int coordtype_to_int(CoordType);

inline string cardpoint(bool, bool);
inline string cardi_b(bool);


/// __________________________________________________
/// __________________________________________________
/// FamousFive -- Templated

/// __________________________________________________
/// Default empty struct for SFINAE	
template<DVecType type>
struct FamousFive {};

/// __________________________________________________
/// Specialised struct for decimal degrees	
template<>
struct FamousFive<DecDegVecDouble> {
#if DEBUG > 0
	FamousFive<DecDegVecDouble>() { _ctrsgn(typeid(*this)); };
	~FamousFive<DecDegVecDouble>() { _ctrsgn(typeid(*this), false); };
#endif
	int get_deg(double x) const { return int(x); }
	double get_decdeg(double x) const { return x; }
	int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
	double get_decmin(double x) const { return polish(mod1by60(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

/// __________________________________________________
/// Specialised struct for degrees and minutes
template<>
struct FamousFive<DegMinVecDouble> {
#if DEBUG > 0
	FamousFive<DegMinVecDouble>() { _ctrsgn(typeid(*this)); };
	~FamousFive<DegMinVecDouble>() { _ctrsgn(typeid(*this), false); };
#endif
	int get_deg(double x) const { return int(x / 1e2); }
	double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
	int get_min(double x) const { return int(x) % int(1e2); }
	double get_decmin(double x) const { return polish(mod1e2(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

/// __________________________________________________
/// Specialised struct for degrees, minutes and seconds
template<>
struct FamousFive<DegMinSecVecDouble> {
#if DEBUG > 0
	FamousFive<DegMinSecVecDouble>() { _ctrsgn(typeid(*this)); };
	~FamousFive<DegMinSecVecDouble>() { _ctrsgn(typeid(*this), false); };
#endif
	int get_deg(double x) const { return int(x / 1e4); }
	double get_decdeg(double x) const { return int(x / 1e4) + (double)int(fmod(x, 1e4) / 1e2) / 60 + mod1e2(x) / 3600; }
	int get_min(double x) const { return (int(x) % int(1e4)) / 1e2; }
	double get_decmin(double x) const { return int(fmod(x, 1e4) / 1e2) + mod1e2(x) / 60; }
	double get_sec(double x) const { return mod1e2(x); }
};

/// __________________________________________________
/// __________________________________________________
/// Convertidor -- functors for converting formats

/// __________________________________________________
/// Default empty struct for SFINAE	
template<DVecType T, DVecType U>
struct Convertidor{
};

/// __________________________________________________
/// Specialised struct for decimal degrees	
template<DVecType T>
struct Convertidor<T, DecDegVecDouble>{
	FamousFive<T> ff {};
	Convertidor() {}
	double operator()(double n) const { return ff.get_decdeg(n); }
};

/// __________________________________________________
/// Specialised struct for degrees and minutes
template<DVecType T>
struct Convertidor<T, DegMinVecDouble>{
	FamousFive<T> ff {};
	Convertidor() {}
	double operator()(double n) const { return ff.get_deg(n) * 1e2 + ff.get_decmin(n); }
};

/// __________________________________________________
/// Specialised struct for degrees, minutes and seconds
template<DVecType T>
struct Convertidor<T, DegMinSecVecDouble>{
	FamousFive<T> ff {};
	Convertidor() {}
	double operator()(double n) const { return ff.get_deg(n) * 1e4 + ff.get_min(n) * 1e2 + ff.get_sec(n); }
};


/// __________________________________________________
/// __________________________________________________
/// Formateador -- functors for converting formats

/// __________________________________________________
/// Default empty struct for SFINAE	
template<DVecType T, SVecType U>
struct Formateador{
};

/// __________________________________________________
/// Specialised struct for decimal degrees	
template<DVecType T>
struct Formateador<T, DecDegVecString>{
	FamousFive<T> ff {};
	Formateador() {}
	string operator()(double n) const { return fmt::format("{:>{}.{}f}\u00B0", ff.get_decdeg(n), 11, 6); }
};

/// __________________________________________________
/// Specialised struct for degrees and minutes
template<DVecType T>
struct Formateador<T, DegMinVecString>{
	FamousFive<T> ff {};
	Formateador() {}
	string operator()(double n) const { return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) + 
					 fmt::format("{:0>{}.{}f}\u2032", fabs(ff.get_decmin(n)), 7, 4); }
};

/// __________________________________________________
/// Specialised struct for degrees, minutes and seconds
template<DVecType T>
struct Formateador<T, DegMinSecVecString>{
	FamousFive<T> ff {};
	Formateador() {}
	string operator()(double n) const { return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) +
					 fmt::format("{:0>{}}\u2032", abs(ff.get_min(n)), 2) +
					 fmt::format("{:0>{}.{}f}\u2033", fabs(ff.get_sec(n)), 5, 2); }
};


/// __________________________________________________
/// __________________________________________________
/// Coordlet class
template<DVecType T>
class Coordlet {
		T dv;
		const vector<bool> latlon;

	public:
		explicit Coordlet(T&&, const vector<bool>);
		Coordlet(const Coordlet&) = delete;							// Disallow copying
		Coordlet& operator=(const Coordlet&) = delete;				//  ——— ditto ———
		Coordlet(Coordlet&&) = delete;								// Disallow transfer ownership
		Coordlet& operator=(Coordlet&&) = delete;					// Disallow moving
#if DEBUG == 0
		virtual ~Coordlet() = default;
#elif DEBUG > 0
		virtual ~Coordlet() { _ctrsgn(typeid(*this), false); }
#endif
		template<DVecType U>
		vector<double> convert() const;
		template<SVecType U>
		vector<string> format() const;
		const vector<bool> validate() const;
};


/// __________________________________________________
/// __________________________________________________
/// Coords class
template<DVecType T>
class Coords {
		const Coordlet<T> cdlt;
	public:
		explicit Coords(T, const vector<bool>);
		Coords(const Coords&) = delete;								// Disallow copying
		Coords& operator=(const Coords&) = delete;					//  ——— ditto ———
		Coords(Coords&&) = delete;									// Disallow transfer ownership
		Coords& operator=(Coords&&) = delete;						// Disallow moving
#if DEBUG == 0
		~Coords() = default;
#elif DEBUG > 0
		~Coords() { _ctrsgn(typeid(*this), false); }
#endif

		vector<double> convert(CoordType) const;						// Don't make return type const—otherwise makes unnecessary copy
		vector<string> format(CoordType) const;						// Don't make return type const—otherwise makes unnecessary copy
		const vector<bool> validate() const;
};

/// __________________________________________________
/// Template aliases
using CoordsDecDeg = Coords<DecDegVecDouble>;
using CoordsDegMin = Coords<DegMinVecDouble>;
using CoordsDegMinSec = Coords<DegMinSecVecDouble>;

/// __________________________________________________
/// Type Traits

/// coordsdecdeg
template <typename t>
struct is_coordsdecdeg : public std::false_type {};

template <>
struct is_coordsdecdeg<CoordsDecDeg> : public std::true_type {};

template<typename t>
constexpr bool is_coordsdecdeg_v = is_coordsdecdeg<t>::value;

/// coordsdegmin
template <typename t>
struct is_coordsdegmin : public std::false_type {};

template <>
struct is_coordsdegmin<CoordsDegMin> : public std::true_type {};

template<typename t>
constexpr bool is_coordsdegmin_v = is_coordsdegmin<t>::value;

/// coordsdegminsec
template <typename t>
struct is_coordsdegminsec : public std::false_type {};

template <>
struct is_coordsdegminsec<CoordsDegMinSec> : public std::true_type {};

template<typename t>
constexpr bool is_coordsdegminsec_v = is_coordsdegminsec<t>::value;

/// __________________________________________________
/// Concept —— coords_t
template <typename T>
concept coords_t = 
	is_coordsdecdeg_v<T> ||
	is_coordsdegmin_v<T> ||
	is_coordsdegminsec_v<T>;

/// __________________________________________________
/// Instantiate Coords<DVecType> object
template<DVecType T>
inline coords_t auto coordsmaker(NumericVector, vector<bool> = vector<bool>{});

/// __________________________________________________
/// __________________________________________________
/// Switches for Coords<DVecType>
vector<double> convert_switch(const NumericVector, CoordType); 
vector<string> format_switch(const NumericVector, CoordType); 
const vector<bool> validate_switch(const NumericVector); 


/// __________________________________________________
/// __________________________________________________
/// Type aliases
template<class T>
using bisvec = array<vector<T>, 2>;
template<class T>
using bisconstvec = array<const vector<T>, 2>;

/// __________________________________________________
/// __________________________________________________
/// Waypoints class
template<DVecType T>
class Waypoints {
		const Coords<T> crdlat;
		const Coords<T> crdlon;
	public:
		explicit Waypoints(NumericVector, NumericVector);
		Waypoints(const Waypoints&) = delete;					// Disallow copying
		Waypoints& operator=(const Waypoints&) = delete;			//  ——— ditto ———
		Waypoints(Waypoints&&) = delete;						// Disallow transfer ownership
		Waypoints& operator=(Waypoints&&) = delete;				// Disallow moving
#if DEBUG == 0
		~Waypoints() = default;
#elif DEBUG > 0
		~Waypoints() { _ctrsgn(typeid(*this), false); }
#endif
		const bisvec<double> convert(CoordType) const;
		const bisvec<string> format(CoordType) const;
		const bisconstvec<bool> validate() const;
};

/// __________________________________________________
/// Template aliases
using WaypointsDecDeg = Waypoints<DecDegVecDouble>;
using WaypointsDegMin = Waypoints<DegMinVecDouble>;
using WaypointsDegMinSec = Waypoints<DegMinSecVecDouble>;

/// __________________________________________________
/// Type Traits

/// waypointsdecdeg
template <typename t>
struct is_waypointsdecdeg : public std::false_type {};

template <>
struct is_waypointsdecdeg<WaypointsDecDeg> : public std::true_type {};

template<typename t>
constexpr bool is_waypointsdecdeg_v = is_waypointsdecdeg<t>::value;

/// waypointsdegmin
template <typename t>
struct is_waypointsdegmin : public std::false_type {};

template <>
struct is_waypointsdegmin<WaypointsDegMin> : public std::true_type {};

template<typename t>
constexpr bool is_waypointsdegmin_v = is_waypointsdegmin<t>::value;

/// waypointsdegminsec
template <typename t>
struct is_waypointsdegminsec : public std::false_type {};

template <>
struct is_waypointsdegminsec<WaypointsDegMinSec> : public std::true_type {};

template<typename t>
constexpr bool is_waypointsdegminsec_v = is_waypointsdegminsec<t>::value;

/// __________________________________________________
/// Concept —— waypoints_t
template <typename T>
concept waypoints_t = 
	is_waypointsdecdeg_v<T> ||
	is_waypointsdegmin_v<T> ||
	is_waypointsdegminsec_v<T>;

/// __________________________________________________
/// Instantiate Waypoints<DVecType> object
template<DVecType T>
inline waypoints_t auto waypointsmaker(DataFrame);

/// __________________________________________________
/// __________________________________________________
/// Switches for Waypoints<DVecType>
const bisvec<double> convert_switch(const DataFrame, CoordType);
const bisvec<string> format_switch(const DataFrame, CoordType);
const bisconstvec<bool> validate_switch(const DataFrame);


/// __________________________________________________
/// __________________________________________________
/// Validation
bool check_valid(const NumericVector, bool = false);
bool check_valid(const DataFrame, bool = false);
// bool validate(const NumericVector, bool = false);
// bool validate(const DataFrame, bool = false);

// bool validate(const NumVec_or_DataFrame auto, bool = false);

template<NumVec_or_DataFrame T>
bool validate(const T, bool = false);

bool valid_ll(const DataFrame);

/// __________________________________________________
/// __________________________________________________
/// Exported functions
NumericVector as_coords(NumericVector, int);
NumericVector convertcoords(const NumericVector, int);
NumericVector latlon(NumericVector, LogicalVector);
NumericVector validatecoords(const NumericVector, const bool);
CharacterVector formatcoords(const NumericVector, bool, bool, int);
DataFrame as_waypointsdefault(DataFrame, int);
DataFrame convertwaypoints(DataFrame, int);
DataFrame validatewaypoints(DataFrame, bool);
CharacterVector formatwaypoints(DataFrame, bool, bool, int);
CharacterVector ll_headers(int, int);
NumericVector as_coordswaypoints(DataFrame, bool);


#endif  // COORDBASE_H_