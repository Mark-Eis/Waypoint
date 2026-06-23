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
	VecTypeBase(const vector<T>& vt) : vector<T>{ vt }								// copy constructor
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

	VecTypeBase& operator=(const VecTypeBase&) = delete;								// copy assignment
	VecTypeBase& operator=(const vector<T>& vt)										// copy assignment
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

	VecTypeBase(vector<T>&& vt) : vector<T>{ std::move(vt) }							// move constructor
	{
#if DEBUG > 0
		_ctrsgn(typeid(*this)); fmt::print("\t(vector<T>&&)\n");
#endif
	}
	VecTypeBase(NumericVector&& vt) = delete;											// move constructor - not defaultable

#if DEBUG == 0
	VecTypeBase& operator=(VecTypeBase&&) = default;										// move assignment
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
/// Class forward declarations
enum class CoordType : char;
class Waypoints;


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
/// FamousFive -- Templated and OO

/// __________________________________________________
/// Abstract base class with pure virtual functions	
struct FamousFive0 {
#if DEBUG > 0
	FamousFive0() { _ctrsgn(typeid(*this)); };
#endif
	virtual ~FamousFive0() = 0;
	virtual int get_deg(double x) const = 0;
	virtual double get_decdeg(double x) const = 0;
	virtual int get_min(double x) const = 0;
	virtual double get_decmin(double x) const = 0;
	virtual double get_sec(double x) const = 0;
};

/// __________________________________________________
/// Default empty derived struct for SFINAE	
template<DVecType type>
struct FamousFive final : FamousFive0 {};

/// __________________________________________________
/// Specialised derived struct for decimal degrees	
template<>
struct FamousFive<DecDegVecDouble> final : FamousFive0 {
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
/// Specialised derived struct for degrees and minutes
template<>
struct FamousFive<DegMinVecDouble> final : FamousFive0 {
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
/// Specialised derived struct for degrees, minutes and seconds
template<>
struct FamousFive<DegMinSecVecDouble> final : FamousFive0 {
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
/// Default struct for SFINAE	
template<DVecType T, DVecType U>
struct Convertidor{
	FamousFive<T> ff {};
	Convertidor() {}
	double operator()(double n) const { stop("Convertidor<DVecType, DVecType>::operator()(double) const my bad"); return -1; }
};

/// __________________________________________________
/// Specialised derived struct for decimal degrees	
template<DVecType T>
struct Convertidor<T, DecDegVecDouble>{
	FamousFive<T> ff {};
	Convertidor() {}
	double operator()(double n) const { return ff.get_decdeg(n); }
};

/// __________________________________________________
/// Specialised derived struct for degrees and minutes
template<DVecType T>
struct Convertidor<T, DegMinVecDouble>{
	FamousFive<T> ff {};
	Convertidor() {}
	double operator()(double n) const { return ff.get_deg(n) * 1e2 + ff.get_decmin(n); }
};

/// __________________________________________________
/// Specialised derived struct for degrees, minutes and seconds
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
/// Default struct for SFINAE	
template<DVecType T, SVecType U>
struct Formateador{
	FamousFive<T> ff {};
	Formateador() {}
	string operator()(double n) const { stop("Formateador<DVecType, DVecType>::operator()(double) const my bad"); return ""s; }
};

/// __________________________________________________
/// Specialised derived struct for decimal degrees	
template<DVecType T>
struct Formateador<T, DecDegVecString>{
	FamousFive<T> ff {};
	Formateador() {}
	string operator()(double n) const { return fmt::format("{:>{}.{}f}\u00B0", ff.get_decdeg(n), 11, 6); }
};

/// __________________________________________________
/// Specialised derived struct for degrees and minutes
template<DVecType T>
struct Formateador<T, DegMinVecString>{
	FamousFive<T> ff {};
	Formateador() {}
	string operator()(double n) const { return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) + 
					 fmt::format("{:0>{}.{}f}\u2032", fabs(ff.get_decmin(n)), 7, 4); }
};

/// __________________________________________________
/// Specialised derived struct for degrees, minutes and seconds
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
		void report() const;											// Temporary —— delete
};


/// __________________________________________________
/// __________________________________________________
/// CrdWptBase class
class CrdWptBase {
	public:
		explicit CrdWptBase();
		CrdWptBase(const CrdWptBase&) = delete;						// Disallow copying
		CrdWptBase& operator=(const CrdWptBase&) = delete;			//  ——— ditto ———
		CrdWptBase(CrdWptBase&&) = delete;							// Disallow transfer ownership
		CrdWptBase& operator=(CrdWptBase&&) = delete;				// Disallow moving
		virtual ~CrdWptBase() = 0;

		virtual const vector<double> convert(CoordType) const = 0;
		virtual vector<string> format(CoordType) const = 0;
		virtual const vector<bool> validate() const = 0;
		virtual void report() const = 0;							// Temporary —— delete
};


/// __________________________________________________
/// __________________________________________________
/// Coords class
template<DVecType T>
class Coords : public CrdWptBase {
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

		const vector<double> convert(CoordType) const;
		vector<string> format(CoordType) const;
		const vector<bool> validate() const;
		void report() const;											// Temporary —— delete
};


/// __________________________________________________
/// __________________________________________________
/// Instantiate Coords<DVecType> object with unique_ptr to base
unique_ptr<CrdWptBase> coordsmaker(NumericVector);

/*
/// __________________________________________________
/// TBC  !!!!!!!!!!!!!!
/// __________________________________________________
/// Waypoints class
class Waypoints : public CrdWptBase {
		DataFrame df;
		NumericVector nvlat;
		NumericVector nvlon;
		vector<bool> validlat { false };
		vector<bool> validlon { false };

		void suffix_nesw(vector<string>&, bool) const;
	public:

		explicit Waypoints(DataFrame);
		Waypoints(const Waypoints&) = delete;					// Disallow copying
		Waypoints& operator=(const Waypoints&) = delete;			//  ——— ditto ———
		Waypoints(Waypoints&&) = delete;						// Disallow transfer ownership
		Waypoints& operator=(Waypoints&&) = delete;				// Disallow moving
		~Waypoints();

		void convert(CoordType);
		vector<string> format(CoordType) const;
		const bool validate() const;
};
*/

/// __________________________________________________
/// __________________________________________________
/// Validation
bool check_valid(const NumericVector);
bool check_valid(const DataFrame);
bool revalidate(const NumericVector);
bool revalidate(const DataFrame);
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