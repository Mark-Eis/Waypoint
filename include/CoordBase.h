/// __________________________________________________
/// CoordBase.h
/// __________________________________________________

#ifndef COORDBASE_H_
#define COORDBASE_H_

/// __________________________________________________
/// Class and Function declarations

// Development and debugging
void _ctrsgn(const type_info&, bool);
class Demangler;
ostream& operator<< (ostream&, const Demangler&);

// Formula simplification
inline double mod1by60(double);
inline double mod1e2(double);
inline double round2(double, int);
inline double polish(double);

// Utility
template<class T, class U> 
inline vector<U> get_vec_attr(const T&, const char*);
template<class T>
inline int get_fmt_attribute(const T&);
template<class T>
inline void checkinherits(T&, const char*);
template<class T>
inline bool is_item_in_obj(const T, const int);
inline void stdlenstr(vector<string>&);
template<class T>
inline void prefixvecstr(vector<string>&, const vector<T>&);
inline bool prefixwithnames(vector<string>&, RObject&);
inline string str_tolower(string);
template<class T>
int nameinobj(const T, const char*);
RObject getnames(const DataFrame);

//CoordType
enum class CoordType : char { decdeg, degmin, degminsec };

inline const CoordType get_coordtype(const int);
template<class T>
inline const CoordType get_coordtype(const T&);
inline const int coordtype_to_int(CoordType);

inline string cardpoint(bool, bool);
inline string cardi_b(bool);

//FamousFive
struct FamousFive;
struct FF_decdeg;
struct FF_degmin;
struct FF_degminsec;

//Convertor
template<CoordType type>
class Convertor;

template<CoordType type>
class Format;

template<class T, CoordType type>
class FormatLL;

class Validator;

//CoordType switches
template<class T, class U>
void convert_switch(T, CoordType);
template<class T>
vector<string> format_switch(const T&, CoordType);

// Coordbase
class Coordbase;

// Coord
class Coord;

// WayPoint
class WayPoint;

// Validation
bool check_valid(const NumericVector);
bool check_valid(const DataFrame);

template<class T>
bool validated(T, const char*, bool&);

template<class T, class U>
const T revalidate(const T);

constexpr auto revalid_Coord = &revalidate<NumericVector, Coord>;
constexpr auto revalid_WayPoint = &revalidate<DataFrame, WayPoint>;

template<class T, class U>
inline const T validate(const T);

bool valid_ll(const DataFrame);

// Exported
NumericVector as_coords(NumericVector, const int);
NumericVector convertcoords(NumericVector, const int);
NumericVector latlon(NumericVector, LogicalVector);
NumericVector validatecoords(NumericVector, const bool);
CharacterVector formatcoords(NumericVector, bool);
DataFrame as_waypointsdefault(DataFrame, const int);
DataFrame convertwaypoints(DataFrame, const int);
DataFrame validatewaypoints(DataFrame, const bool);
CharacterVector formatwaypoints(DataFrame, bool);
CharacterVector ll_headers(int, const int);
NumericVector as_coordswaypoints(DataFrame, bool);



/*
#define FMT_HEADER_ONLY
#include </opt/homebrew/Cellar/fmt/11.1.4/include/fmt/base.h>		// verbose path needs sorting!

/// Either use format_as() function or specialise the formatter struct template.
// #define USE_FORMAT_AS

int fmtdemo();
int fmtdemo2(int);
std::string fmtdemo3(int);

/// __________________________________________________
/// CoordType enum class

namespace CoordType_enum {

	enum class CoordType : char { decdeg, degmin, degminsec };
	
	inline const CoordType get_coordtype(int i);


#ifdef USE_FORMAT_AS
// Use format_as()

	auto format_as(CoordType ct) -> std::string;
#endif  // ifdef USE_FORMAT_AS

}

#ifndef USE_FORMAT_AS
// Use formatter struct template specialisation (can't be within namespace CoordType_enum)

template <> struct fmt::formatter<CoordType_enum::CoordType>: formatter<string_view> {
	// parse is inherited from formatter<string_view>.

	auto format(CoordType_enum::CoordType, format_context&) const
		-> format_context::iterator;
};
#endif  // ifndef USE_FORMAT_AS
*/

#endif  // COORDBASE_H_
