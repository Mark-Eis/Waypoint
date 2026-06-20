/// __________________________________________________
/// CoordBase.h
/// __________________________________________________

#ifndef COORDBASE_H_
#define COORDBASE_H_

#define FMT_HEADER_ONLY
#include "fmt/base.h"		// …fmt/*.h copied to …/R/Packages/Waypoint/src.


/// __________________________________________________
/// __________________________________________________
/// Development and debugging


#define DEBUG 0

#if DEBUG > 0

void _ctrsgn(const std::type_info&, bool = false);

/// __________________________________________________
/// Check address for unnecessary copying
inline void* address(const auto& t)
{
	return (void*)&t;
}

/// __________________________________________________
/// Format string for debugging code
constexpr auto exportstr { "——Rcpp::export——"sv };

#endif	// if DEBUG > 0

const string demangle(const std::type_info&);


/// __________________________________________________
/// __________________________________________________
/// Class and Function declarations

/// __________________________________________________
/// Class forward declarations
enum class CoordType : char;
class Coordlet;
class CrdWptBase;
class Coords;
class Waypoints;

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
/// Concept
template<typename T>
concept Coords_or_Waypoints =
	requires (T t) {
		{ t.convert(CoordType::decdeg) } -> std::same_as<void>;
		{ t.convert(CoordType::degmin) } -> std::same_as<void>;
		{ t.convert(CoordType::degminsec) } -> std::same_as<void>;
		{ t.format(CoordType::decdeg) } -> std::same_as<vector<string>>;
		{ t.format(CoordType::degmin) } -> std::same_as<vector<string>>;
		{ t.format(CoordType::degminsec) } -> std::same_as<vector<string>>;
		{ t.validate() } -> std::same_as<bool>;
	};


/// __________________________________________________
/// CoordType access functions
inline const CoordType get_coordtype(int);
inline const CoordType get_coordtype(const NumVec_or_DataFrame auto&);
inline int coordtype_to_int(CoordType);

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
// template<NumVec_or_DataFrame T, typename U> 
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
inline string cardpoint(bool, bool);
inline string cardi_b(bool);


/// __________________________________________________
/// __________________________________________________
/// FamousFive -- Templated and OO

/// __________________________________________________
/// Abstract base class with pure virtual functions	
struct FamousFive0 {
	virtual int get_deg(double x) const = 0;
	virtual double get_decdeg(double x) const = 0;
	virtual int get_min(double x) const = 0;
	virtual double get_decmin(double x) const = 0;
	virtual double get_sec(double x) const = 0;
	virtual ~FamousFive0() = 0;
};

/// __________________________________________________
/// Default empty derived struct for SFINAE	
template<CoordType type>
struct FamousFive final : FamousFive0 {};

/// __________________________________________________
/// Specialised derived struct for decimal degrees	
template<>
struct FamousFive<CoordType::decdeg> final : FamousFive0 {
	int get_deg(double x) const { return int(x); }
	double get_decdeg(double x) const { return x; }
	int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
	double get_decmin(double x) const { return polish(mod1by60(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

/// __________________________________________________
/// Specialised derived struct for degrees and minutes
template<>
struct FamousFive<CoordType::degmin> final : FamousFive0 {
	int get_deg(double x) const { return int(x / 1e2); }
	double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
	int get_min(double x) const { return int(x) % int(1e2); }
	double get_decmin(double x) const { return polish(mod1e2(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

/// __________________________________________________
/// Specialised derived struct for degrees, minutes and seconds
template<>
struct FamousFive<CoordType::degminsec> final : FamousFive0 {
	int get_deg(double x) const { return int(x / 1e4); }
	double get_decdeg(double x) const { return int(x / 1e4) + (double)int(fmod(x, 1e4) / 1e2) / 60 + mod1e2(x) / 3600; }
	int get_min(double x) const { return (int(x) % int(1e4)) / 1e2; }
	double get_decmin(double x) const { return int(fmod(x, 1e4) / 1e2) + mod1e2(x) / 60; }
	double get_sec(double x) const { return mod1e2(x); }
};


/// __________________________________________________
/// __________________________________________________
/// Coordlet class
// template<CoordType current_type>
class Coordlet {
		unique_ptr<FamousFive0> ff;
		NumericVector nv;
		const vector<bool> latlon;

		unique_ptr<FamousFive0> switch_ff(NumericVector);
	public:
		explicit Coordlet(NumericVector);
		Coordlet(const Coordlet&) = delete;						// Disallow copying
		Coordlet& operator=(const Coordlet&) = delete;			//  ——— ditto ———
		Coordlet(Coordlet&&) = delete;							// Disallow transfer ownership
		Coordlet& operator=(Coordlet&&) = delete;				// Disallow moving
		virtual ~Coordlet() = default;
//		virtual ~Coordlet() { fmt::print("§Coordlet::~Coordlet() "); _ctrsgn(typeid(*this), true); }

		void convert(CoordType);
		vector<string> format(CoordType) const;
		const vector<bool> validate() const;
};


/// __________________________________________________
/// CrdWptBase class
class CrdWptBase {
	protected:
		const CoordType ct;
	public:
		explicit CrdWptBase(CoordType);
		CrdWptBase(const CrdWptBase&) = delete;						// Disallow copying
		CrdWptBase& operator=(const CrdWptBase&) = delete;			//  ——— ditto ———
		CrdWptBase(CrdWptBase&&) = delete;							// Disallow transfer ownership
		CrdWptBase& operator=(CrdWptBase&&) = delete;				// Disallow moving
		virtual ~CrdWptBase() = 0;

		virtual void convert(CoordType) = 0;
		virtual vector<string> format(CoordType) const = 0;
		virtual const bool validate() const = 0;
};


/// __________________________________________________
/// Coords class
class Coords : public CrdWptBase {
		NumericVector nv;
		vector<bool> valid { false };
		void suffix_nesw(vector<string>&) const;
		void suffix_latlon(vector<string>&) const;
	public:
		explicit Coords(NumericVector);
		Coords(const Coords&) = delete;						// Disallow copying
		Coords& operator=(const Coords&) = delete;			//  ——— ditto ———
		Coords(Coords&&) = delete;							// Disallow transfer ownership
		Coords& operator=(Coords&&) = delete;				// Disallow moving
		~Coords() = default;
//		~Coords() { fmt::print("§Coords::~Coords(); {}; ", ct); _ctrsgn(typeid(*this), true); }

		void convert(CoordType);
		vector<string> format(CoordType) const;
		const bool validate() const;
};


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


/// __________________________________________________
/// __________________________________________________
/// Validation
bool check_valid(const NumericVector);
bool check_valid(const DataFrame);
template<NumVec_or_DataFrame T>
bool revalidate(const T);

bool valid_ll(const DataFrame);

/// __________________________________________________
/// __________________________________________________
/// Exported functions
NumericVector as_coords(NumericVector, int);
NumericVector convertcoords(NumericVector, int);
NumericVector latlon(NumericVector, LogicalVector);
NumericVector validatecoords(const NumericVector, const bool);
CharacterVector formatcoords(NumericVector, bool, bool, int);
DataFrame as_waypointsdefault(DataFrame, int);
DataFrame convertwaypoints(DataFrame, int);
DataFrame validatewaypoints(DataFrame, bool);
CharacterVector formatwaypoints(DataFrame, bool, bool, int);
CharacterVector ll_headers(int, int);
NumericVector as_coordswaypoints(DataFrame, bool);


#endif  // COORDBASE_H_