/// __________________________________________________
/// CoordBase.h
/// __________________________________________________

#ifndef COORDBASE_H_
#define COORDBASE_H_

#define FMT_HEADER_ONLY
// #include "fmt/base.h"		// …fmt/*.h copied to …/R/Packages/Waypoint/src.  Works, but not in pkgdown
#include "/Users/frzmce/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/R/Packages/Waypoint/src/fmt/format.h"

/// __________________________________________________
/// __________________________________________________
/// Class and Function declarations

/// __________________________________________________
/// Class forward declarations
enum class CoordType : char;
template<CoordType current_type>
class Coordlet;
class Coords;
class Waypoints;

/// __________________________________________________
/// Concept
template <typename T>
concept NumericVector_or_DataFrame = 
	std::is_same_v<NumericVector, T> || std::is_same_v<const NumericVector, T> ||
	std::is_same_v<DataFrame, T> || std::is_same_v<const DataFrame, T>;

template <typename T>
concept List_or_DataFrame = std::is_same_v<List, T> || std::is_same_v<DataFrame, T>;


/// __________________________________________________
/// __________________________________________________
/// Development and debugging
void _ctrsgn(const std::type_info&, bool = false);
const string demangle(const std::type_info&);

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
template<NumericVector_or_DataFrame T, class U> 
inline vector<U> get_vec_attr(const T&, const char*);
template<NumericVector_or_DataFrame T>
inline int get_fmt_attribute(const T&);
template<NumericVector_or_DataFrame T>
int check_logical_attr(T t, const char* attrname);
template<NumericVector_or_DataFrame T>
inline void checkinherits(T&, const char*);
template<class T>
inline bool is_item_in_obj(const T, int);
inline void stdlenstr(vector<string>&);
inline void concat_vecstr_elmnts(const vector<string>&, vector<string>&, const string = " ");
inline void concat_vecstr_elmnts(const vector<int>&, vector<string>&, const string = " ");
inline bool prefixwithnames(vector<string>&, RObject&);
inline string str_tolower(string);
template<List_or_DataFrame T>
int nameinobj(const T, const char*);
RObject getnames(const DataFrame);

/// __________________________________________________
/// __________________________________________________
/// CoordType enum
enum class CoordType : char { decdeg, degmin, degminsec };

template <>
struct fmt::formatter<CoordType>: formatter<string_view>
{
	auto format(CoordType, format_context&) const
		-> format_context::iterator;
};

/// __________________________________________________
/// __________________________________________________
/// CoordType Type Traits

/// __________________________________________________
/// CoordType::decdeg
template <auto T>
struct isDecDeg : public std::false_type {};

template <>
struct isDecDeg<CoordType::decdeg> : public std::true_type {};

template<auto T>
constexpr bool isDecDeg_v = isDecDeg<T>::value;

/// __________________________________________________
/// CoordType::degmin
template <auto T>
struct isDegMin : public std::false_type {};

template <>
struct isDegMin<CoordType::degmin> : public std::true_type {};

template<auto T>
constexpr bool isDegMin_v = isDegMin<T>::value;

/// __________________________________________________
/// CoordType::degminsec
template <auto T>
struct isDegMinSec : public std::false_type {};

template <>
struct isDegMinSec<CoordType::degminsec> : public std::true_type {};

template<auto T>
constexpr bool isDegMinSec_v = isDegMinSec<T>::value;

/// __________________________________________________
/// CoordType access functions
inline const CoordType get_coordtype(int);
template<NumericVector_or_DataFrame T>
inline const CoordType get_coordtype(const T&);
inline int coordtype_to_int(CoordType);

inline string cardpoint(bool, bool);
inline string cardi_b(bool);


/// __________________________________________________
/// Concept
template <typename T>
concept Coords_or_Waypoints =
	requires (T t) {
		t.template convert<CoordType::decdeg>();
		t.template convert<CoordType::degmin>();
		t.template convert<CoordType::degminsec>();
		t.validate();
//		t.format(CoordType);
//		t.format();
	};


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
template<CoordType current_type>
class Coordlet {

		FamousFive<current_type> ff;
		NumericVector nv;
		const vector<bool> latlon;

		template<CoordType>
		void convert0();
		template<CoordType> 
		vector<string> format0() const;
	public:
		explicit Coordlet(NumericVector);
		Coordlet(const Coordlet&) = delete;						// Disallow copying
		Coordlet& operator=(const Coordlet&) = delete;			//  ——— ditto ———
		Coordlet(Coordlet&&) = delete;							// Disallow transfer ownership
		Coordlet& operator=(Coordlet&&) = delete;				// Disallow moving
		virtual ~Coordlet() = default;
//		virtual ~Coordlet() { fmt::print("§Coordlet::~Coordlet(); {}; ", current_type); _ctrsgn(typeid(*this), true); }
		vector<string> format_switch(CoordType) const;
		void convert_switch(CoordType);
		const vector<bool> validate() const;
};


/// __________________________________________________
/// Coords class
class Coords {
	protected:
		const CoordType ct;
		NumericVector nv;
		vector<bool> valid { false };
//		template<CoordType>
//		void format_suffix(vector<string>&) const;  // Possibly?
		void format_suffix(vector<string>&, CoordType) const;
	public:
		explicit Coords(NumericVector);
		Coords(const Coords&) = delete;						// Disallow copying
		Coords& operator=(const Coords&) = delete;			//  ——— ditto ———
		Coords(Coords&&) = delete;							// Disallow transfer ownership
		Coords& operator=(Coords&&) = delete;				// Disallow moving
		~Coords() = default;
//		virtual ~Coords() { fmt::print("§Coords::~Coords(); {}; ", ct); _ctrsgn(typeid(*this), true); }

		template<CoordType>
		void convert();
		const bool validate() const;
		vector<string> format(CoordType) const;
//		void format_suffix_switch(CoordType) const;  // Possibly?
};


/// __________________________________________________
/// Waypoints class
class Waypoints {
	protected:
		const CoordType ct;
		DataFrame df;
		NumericVector nvlat;
		NumericVector nvlon;
		vector<bool> validlat { false };
		vector<bool> validlon { false };
		void format_suffix(vector<string>&, const bool) const;
	public:
		explicit Waypoints(DataFrame);
		Waypoints(const Waypoints&) = delete;					// Disallow copying
		Waypoints& operator=(const Waypoints&) = delete;			//  ——— ditto ———
		Waypoints(Waypoints&&) = delete;						// Disallow transfer ownership
		Waypoints& operator=(Waypoints&&) = delete;				// Disallow moving
		~Waypoints();

		template<CoordType>
		void convert();
		const bool validate() const;
		vector<string> format(CoordType) const;
};


/// __________________________________________________
/// __________________________________________________
/// CoordType switches
vector<string> format_switch_current(NumericVector, const CoordType);
template<CoordType current_type> 
vector<string> format_dispatch(NumericVector, const CoordType);

void convert_switch_current(NumericVector, const CoordType);
template<CoordType current_type> 
inline void convert_dispatch(NumericVector, const CoordType);

const vector<bool> validate_switch_current(const NumericVector);
template<CoordType current_type> 
const vector<bool> validate_dispatch(const NumericVector);


/// __________________________________________________
/// __________________________________________________
/// Validation
bool check_valid(const NumericVector);
bool check_valid(const DataFrame);

template<NumericVector_or_DataFrame T, Coords_or_Waypoints U>
bool revalidate(const T);

constexpr auto revalid_Coords = &revalidate<NumericVector, Coords>;
constexpr auto revalid_WayPoints = &revalidate<DataFrame, Waypoints>;

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