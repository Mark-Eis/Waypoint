/// __________________________________________________
/// CoordBase.h
/// __________________________________________________

#ifndef COORDBASE_H_
#define COORDBASE_H_

#define FMT_HEADER_ONLY
#include </opt/homebrew/Cellar/fmt/11.1.4/include/fmt/base.h>		// verbose path needs sorting!

/// __________________________________________________
/// Class and Function declarations

// Development and debugging
void _ctrsgn(const std::type_info&, bool);
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
template <>
//struct fmt::formatter<CoordType_enum::CoordType>: formatter<string_view>
struct fmt::formatter<CoordType>: formatter<string_view>
{
//	auto format(CoordType_enum::CoordType, format_context&) const
	auto format(CoordType, format_context&) const
		-> format_context::iterator;
};

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
vector<string> format_switch(const T&);

// Coordbase class
class Coordbase {
	protected:
		CoordType ct;
		const FamousFive& ff;

	public:
		Coordbase(CoordType _ct);
		Coordbase(const Coordbase&) = delete;						// Disallow copying
		Coordbase& operator=(const Coordbase&) = delete;			//  ——— ditto ———
		Coordbase(Coordbase&&) = delete;							// Disallow transfer ownership
		Coordbase& operator=(Coordbase&&) = delete;					// Disallow moving
		virtual ~Coordbase() = 0;
		CoordType get_coordtype() const;
};


// Coordinate derived class
class Coord : public Coordbase {
	protected:
		const NumericVector nv;
		const vector<bool> valid { false };
		const vector<bool> latlon;

	public:
		Coord(CoordType, const NumericVector);
		~Coord() = default;
//		~Coord() { cout << "§Coord::~Coord() "; _ctrsgn(typeid(*this), true); }

		template<CoordType type>
		void convert() const;
		void validate(bool warn = true) const;
		template<CoordType type>
		vector<string> format_ct() const;
		vector<string> format(bool usenames) const;
};


// Waypoint derived class
class WayPoint : public Coordbase {
	protected:
		const DataFrame df;
		const NumericVector nvlat;
		const NumericVector nvlon;
		const vector<bool> validlat { false };
		const vector<bool> validlon { false };
	public:
		explicit WayPoint(CoordType, const DataFrame);
		~WayPoint() = default;
//		~WayPoint() { cout << "§WayPoint::~WayPoint() "; _ctrsgn(typeid(*this), true); }

		template<CoordType type>
		void convert() const;
		void validate(bool = true) const;
		template<CoordType type>
		vector<string> format_ct() const;
		vector<string> format(bool usenames) const;
};


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
