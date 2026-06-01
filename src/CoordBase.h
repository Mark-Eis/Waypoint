/// __________________________________________________
/// Waypoint.h
/// __________________________________________________

#ifndef Waypoint_H_
#define Waypoint_H_

#define FMT_HEADER_ONLY
// #include "fmt/base.h"		// …fmt/*.h copied to …/R/Packages/Waypoint/src.  Works, but not in pkgdown
#include "/Users/frzmce/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/R/Packages/Waypoint/src/fmt/format.h"

/// __________________________________________________
/// Class and Function declarations

/// Concept
template <typename T>
concept NumericVector_or_DataFrame = std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value;

template <typename T>
concept List_or_DataFrame = std::is_same<List, T>::value || std::is_same<DataFrame, T>::value;

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
template<class T>
inline void prefixvecstr(vector<string>&, const vector<T>&);
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

inline const CoordType get_coordtype(int);
template<NumericVector_or_DataFrame T>
inline const CoordType get_coordtype(const T&);
inline int coordtype_to_int(CoordType);

inline string cardpoint(bool, bool);
inline string cardi_b(bool);

/// __________________________________________________
/// __________________________________________________
/// FamousFiveT Templated and OO
struct FamousFive {
	virtual int get_deg(double x) const = 0;
	virtual double get_decdeg(double x) const = 0;
	virtual int get_min(double x) const = 0;
	virtual double get_decmin(double x) const = 0;
	virtual double get_sec(double x) const = 0;
};

/// __________________________________________________
/// Default struct for decimal degrees	
template<CoordType type>
struct FamousFiveT final : FamousFive {
	int get_deg(double x) const { return int(x); }
	double get_decdeg(double x) const { return x; }
	int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
	double get_decmin(double x) const { return polish(mod1by60(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

/// __________________________________________________
/// Specialised struct for degrees and minutes
template<>
struct FamousFiveT<CoordType::degmin> final : FamousFive {
	int get_deg(double x) const { return int(x / 1e2); }
	double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
	int get_min(double x) const { return int(x) % int(1e2); }
	double get_decmin(double x) const { return polish(mod1e2(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

/// __________________________________________________
/// Specialised struct for degrees, minutes and seconds
template<>
struct FamousFiveT<CoordType::degminsec> final : FamousFive {
	int get_deg(double x) const { return int(x / 1e4); }
	double get_decdeg(double x) const { return int(x / 1e4) + (double)int(fmod(x, 1e4) / 1e2) / 60 + mod1e2(x) / 3600; }
	int get_min(double x) const { return (int(x) % int(1e4)) / 1e2; }
	double get_decmin(double x) const { return int(fmod(x, 1e4) / 1e2) + mod1e2(x) / 60; }
	double get_sec(double x) const { return mod1e2(x); }
};


/// __________________________________________________
/// __________________________________________________
/// Class forward declarations
template<CoordType current_type>
class Coordlet;
template<CoordType current_type>
class WayPoint;

/// __________________________________________________
/// Concept
template <typename T>
concept Coord_or_WayPoint =
	requires (T t) {
		t.template convert<CoordType::decdeg>();
		t.template convert<CoordType::degmin>();
		t.template convert<CoordType::degminsec>();
		t.template format<CoordType::decdeg>();
		t.template format<CoordType::degmin>();
		t.template format<CoordType::degminsec>();
		t.get_coordtype();
		t.validate();
	};


/// __________________________________________________
/// __________________________________________________
/// Coordlet class
template<CoordType current_type>
class Coordlet {

		FamousFiveT<current_type> ff;
		NumericVector nv;
		vector<bool> valid { false };
		const vector<bool> latlon;

		template<CoordType required_type>
		void convert0();
		template<CoordType required_type> 
		vector<string> format0(bool wpt) const;
	public:
		explicit Coordlet(NumericVector nv);
		Coordlet(const Coordlet&) = delete;						// Disallow copying
		Coordlet& operator=(const Coordlet&) = delete;			//  ——— ditto ———
		Coordlet(Coordlet&&) = delete;							// Disallow transfer ownership
		Coordlet& operator=(Coordlet&&) = delete;				// Disallow moving
		virtual ~Coordlet() = default;
		vector<string> format_switch(CoordType required_type, bool wpt) const;
		void convert_switch(CoordType required_type);
		void validate(bool = true, const char* = "");
};


/// __________________________________________________
/// __________________________________________________
/// CoordType switches
vector<string> format_switch_current(NumericVector, CoordType, CoordType, bool = false);
template<CoordType current_type> 
vector<string> format_required(NumericVector, CoordType, bool);

void convert_switch_current(NumericVector, CoordType, CoordType);
template<CoordType current_type> 
inline void convert_required(NumericVector, CoordType);

void validate_switch_current(NumericVector, CoordType, bool, const char* = "");


/// __________________________________________________
/// __________________________________________________
/// Validation
bool check_valid(const NumericVector);
bool check_valid(const DataFrame);

template<NumericVector_or_DataFrame T /*, Coord_or_WayPoint U */>
bool revalidate(const T);

// constexpr auto revalid_Coord = &revalidate<NumericVector, Coord>;
// constexpr auto revalid_WayPoint = &revalidate<DataFrame, WayPoint>;

template<NumericVector_or_DataFrame T /*, Coord_or_WayPoint U */>
inline const T validate(const T);

bool valid_ll(const DataFrame);

/// __________________________________________________
/// __________________________________________________
/// Exported functions
NumericVector as_coords(NumericVector, int);
NumericVector convertcoords(NumericVector, int);
NumericVector latlon(NumericVector, LogicalVector);
NumericVector validatecoords(NumericVector, bool);
CharacterVector formatcoords(NumericVector, bool, bool, int);
DataFrame as_waypointsdefault(DataFrame, int);
DataFrame convertwaypoints(DataFrame, int);
DataFrame validatewaypoints(DataFrame, bool);
CharacterVector formatwaypoints(DataFrame, bool, bool, int);
CharacterVector ll_headers(int, int);
NumericVector as_coordswaypoints(DataFrame, bool);


#endif  // Waypoint_H_