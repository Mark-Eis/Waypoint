/// __________________________________________________
/// CoordBase.h
/// __________________________________________________

#ifndef COORDBASE_H_
#define COORDBASE_H_

#define FMT_HEADER_ONLY
#include </opt/homebrew/Cellar/fmt/11.1.4/include/fmt/base.h>		// verbose path needs sorting!

/// __________________________________________________
/// Class and Function declarations

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

inline const CoordType get_coordtype(const int);
template<class T>
inline const CoordType get_coordtype(const T&);
inline const int coordtype_to_int(CoordType);

inline string cardpoint(bool, bool);
inline string cardi_b(bool);

/// __________________________________________________
/// __________________________________________________
/// FamousFive Class and Derived Classes
struct FamousFive {
//	FamousFive() { fmt::print("§{} ", "FamousFive()"); _ctrsgn(typeid(*this)); }
	virtual ~FamousFive() = 0;	
	virtual int get_deg(double x) const = 0;
	virtual double get_decdeg(double x) const = 0;
	virtual int get_min(double x) const = 0;
	virtual double get_decmin(double x) const = 0;
	virtual double get_sec(double x) const = 0;
};

inline FamousFive::~FamousFive()
{
//	fmt::print("§{} ", "~FamousFive()"); _ctrsgn(typeid(*this), true);
//	std::fflush(nullptr);
}	

/// __________________________________________________
/// Derived class for decimal degrees	
struct FF_decdeg : public FamousFive {
//	FF_decdeg() { fmt::print("§{} ", "FF_decdeg()"); _ctrsgn(typeid(*this)); }
	~FF_decdeg() = default;
//	~FF_decdeg() { fmt::print("§{} ", "~FF_decdeg()"); _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x); }
	double get_decdeg(double x) const { return x; }
	int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
	double get_decmin(double x) const { return polish(mod1by60(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
} ff_decdeg;

/// __________________________________________________
/// Derived class for degrees and minutes
struct FF_degmin : public FamousFive {
//	FF_degmin() { fmt::print("§{} ", "FF_degmin()"); _ctrsgn(typeid(*this)); }
	~FF_degmin() = default;
//	~FF_degmin() { fmt::print("§{} ", "~FF_degmin()"); _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x / 1e2); }
	double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
	int get_min(double x) const { return int(x) % int(1e2); }
	double get_decmin(double x) const { return polish(mod1e2(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
} ff_degmin;

/// __________________________________________________
/// Derived class for degrees, minutes and seconds
struct FF_degminsec : public FamousFive {
//	FF_degminsec() { fmt::print("§{} ", "FF_degminsec()"); _ctrsgn(typeid(*this)); }
	~FF_degminsec() = default;
//	~FF_degminsec() { fmt::print("§{} ", "~FF_degminsec()"); _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x / 1e4); }
	double get_decdeg(double x) const { return int(x / 1e4) + (double)int(fmod(x, 1e4) / 1e2) / 60 + mod1e2(x) / 3600; }
	int get_min(double x) const { return (int(x) % int(1e4)) / 1e2; }
	double get_decmin(double x) const { return int(fmod(x, 1e4) / 1e2) + mod1e2(x) / 60; }
	double get_sec(double x) const { return mod1e2(x); }
} ff_degminsec;


/// __________________________________________________
/// __________________________________________________
/// Class forward declarations
class Coordbase;
class Coord;
class WayPoint;


/// __________________________________________________
/// __________________________________________________
/// Templated coord type conversion functors

template<CoordType type>
class Convertor {
	protected:
		const FamousFive& ff; 
	public:
		Convertor(const FamousFive& _ff) : ff(_ff)
		{
//			fmt::print("§Convertor<CoordType::{}>::Convertor()(const FamousFive&) ", type); _ctrsgn(typeid(*this));
			std::fflush(nullptr);}
		~Convertor() = default;
//		~Convertor() { fmt::print("§Convertor<CoordType::{}>::~Convertor() ", type); _ctrsgn(typeid(*this), true); }
		double operator()(double n);
};


/// __________________________________________________
/// Default operator(), for decimal degrees
template<CoordType type>
inline double Convertor<type>::operator()(double n)
{
//	fmt::print("@Convertor<CoordType::{}>::operator() [default]\n", type);
	return ff.get_decdeg(n);
}


/// __________________________________________________
/// Specialised operator() for degrees and minutes
template<>
inline double Convertor<CoordType::degmin>::operator()(double n)
{
//	fmt::print("@{}\n", "Convertor<CoordType::degmin>::operator()");
	return ff.get_deg(n) * 1e2 + ff.get_decmin(n);
}


/// __________________________________________________
/// Specialised operator() for degrees, minutes and seconds
template<>
inline double Convertor<CoordType::degminsec>::operator()(double n)
{
//	fmt::print("@{}\n", "Convertor<CoordType::degminsec>::operator()");
	return ff.get_deg(n) * 1e4 + ff.get_min(n) * 1e2 + ff.get_sec(n);
}


/// __________________________________________________
/// __________________________________________________
/// Templated coord formatting functors
template<CoordType type>
class Format {
	protected:
		const FamousFive& ff;
	public:
		Format(const FamousFive& _ff) : ff(_ff)
		{
//			fmt::print("§Format<CoordType::{}>::Format()(const FamousFive&) ", type); _ctrsgn(typeid(*this));
//			std::fflush(nullptr);
		}
		~Format() = default;
//		~Format() { fmt::print("§Format<CoordType::{}>::~Format() ", type); _ctrsgn(typeid(*this), true); }
		string operator()(double n) const;
};


/// __________________________________________________
/// Default operator(), for decimal degrees
template<CoordType type>
inline string Format<type>::operator()(double n) const
{
//	fmt::print("@Format<CoordType::{}>::operator() [default]\n", type);
	return fmt::format("{:>{}.{}f}\u00B0", ff.get_decdeg(n), 11, 6);
}

/// __________________________________________________
/// Specialised operator() for degrees and minutes
template<>
inline string Format<CoordType::degmin>::operator()(double n) const
{
//	fmt::print("@Format<CoordType::{}>::operator()\n", CoordType::degmin);
	return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) +
		   fmt::format("{:0>{}.{}f}\u2032", abs(ff.get_decmin(n)), 7, 4);
}

/// __________________________________________________
/// Specialised operator() for degrees, minutes and seconds
template<>
inline string Format<CoordType::degminsec>::operator()(double n) const
{
//	fmt::print("@Format<CoordType::{}>::operator()\n", CoordType::degminsec);
	return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) +
		   fmt::format("{:0>{}}\u2032", abs(ff.get_min(n)), 2) +
		   fmt::format("{:0>{}.{}f}\u2033", abs(ff.get_sec(n)), 5, 2);
}


/// __________________________________________________
/// __________________________________________________
/// Formatting functors for latitude and longitude

/// Default functor for degrees, minutes (and seconds)
template<class T, CoordType type>
class FormatLL {
		const FamousFive& ff; 
		vector<bool>::const_iterator ll_it;
		int ll_size;
	public:
		FormatLL(const FamousFive& _ff, const vector<bool>& ll) : ff(_ff), ll_it(ll.begin()), ll_size(ll.size())
		{
			static_assert(std::is_same<Coord, T>::value || std::is_same<WayPoint, T>::value, "T must be Coord or WayPoint");
//			fmt::print("§FormatLL<{}, CoordType::{}>::FormatLL(const FamousFive&, vector<bool>&) ", "Coord or WayPoint", type); _ctrsgn(typeid(*this));
//			std::fflush(nullptr);
		}
		~FormatLL() = default;
//		~FormatLL() { fmt::print("§FormatLL<{}, CoordType::{}>::~FormatLL() ", "Coord or WayPoint", type); _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
//			fmt::print("@FormatLL<{}, CoordType::{}>::operator(string, double) [default]\n", "Coord or WayPoint", type);
			return ostr += ll_size ? cardpoint(ff.get_decmin(n) < 0, ll_size > 1 ? *ll_it++ : *ll_it) : cardi_b(ff.get_decmin(n) < 0);
		}
};

/// __________________________________________________
/// Specialised functor for decimal degrees Coord
template<>
class FormatLL<Coord, CoordType::decdeg> {
		vector<bool>::const_iterator ll_it;
		int ll_size;
	public:
		FormatLL(const FamousFive& _ff, const vector<bool>& ll) : ll_it(ll.begin()), ll_size(ll.size())
		{
//			fmt::print("§FormatLL<{}, CoordType::{}>::FormatLL(const FamousFive&, vector<bool>&) ", "Coord", CoordType::decdeg ); _ctrsgn(typeid(*this));
//			std::fflush(nullptr);
		}
		~FormatLL() = default;
//		~FormatLL() { fmt::print("§FormatLL<{}, CoordType::{}>::~FormatLL() ", "Coord", CoordType::decdeg); _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
//			fmt::print("@FormatLL<{}, CoordType::{}>::operator(string, double)\n", "Coord", CoordType::decdeg);
			if (ll_size)
				return ostr += ((ll_size > 1 ? *ll_it++ : *ll_it) ? " lat" : " lon");
			else
				return ostr;
		}
};

/// __________________________________________________
/// Specialised functor for decimal degrees WayPoint—simply returns its argument
template<>
class FormatLL<WayPoint, CoordType::decdeg> {
		vector<bool>::const_iterator ll_it;
		int ll_size;
	public:
		FormatLL(const FamousFive& _ff, const vector<bool>& ll) : ll_it(ll.begin()), ll_size(ll.size())
		{
//			fmt::print("§FormatLL<{}, CoordType::{}>::FormatLL(const FamousFive&, vector<bool>&) ", "WayPoint", CoordType::decdeg); _ctrsgn(typeid(*this));
//			std::fflush(nullptr);
		}
		~FormatLL() = default;
//		~FormatLL() { fmt::print("§FormatLL<{}, CoordType::{}>::~FormatLL() ", "WayPoint", CoordType::decdeg); _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
//			fmt::print("@FormatLL<{}, CoordType::{}>::operator(string, double)\n", "WayPoint", CoordType::decdeg);
			return ostr;
		}
};


/// __________________________________________________
/// __________________________________________________
/// Validate functor

class Validator {
		const FamousFive& ff;
		vector<bool>::const_iterator ll_it;
		const int ll_size;
	public:
		Validator(const FamousFive& _ff, const vector<bool>& ll) : ff(_ff), ll_it(ll.begin()), ll_size(ll.size())
		{
//			fmt::print("§{} ", "Validator::Validator(const FamousFive&, vector<bool>&)"); _ctrsgn(typeid(*this));
		}
		~Validator() = default;
//		~Validator() { fmt::print("§{} ", "Validator::~Validator()"); _ctrsgn(typeid(*this), true); }
		bool operator()(double n)
		{
//			fmt::print("@{} validating: {: {}f}\n", "Validator()", n, 9);
			return !((abs(ff.get_decdeg(n)) > (ll_size && (ll_size > 1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(abs(ff.get_decmin(n)) >= 60) ||
				(abs(ff.get_sec(n)) >= 60));
		}
};


/// __________________________________________________
/// __________________________________________________
///CoordType switches
template<class T, class U>
void convert_switch(T, CoordType);
template<class T>
vector<string> format_switch(const T&);


/// __________________________________________________
/// __________________________________________________
/// Coordbase class
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

/// __________________________________________________
/// Coordinate derived class
class Coord : public Coordbase {
	protected:
		const NumericVector nv;
		const vector<bool> valid { false };
		const vector<bool> latlon;

	public:
		Coord(CoordType, const NumericVector);
		~Coord() = default;
//		~Coord() { fmt::print("§{} {} ", "Coord::~Coord()", ct); _ctrsgn(typeid(*this), true); }

		template<CoordType type>
		void convert() const;
		void validate(bool warn = true) const;
		template<CoordType type>
		vector<string> format_ct() const;
		vector<string> format(bool usenames) const;
};

/// __________________________________________________
/// Waypoint derived class
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
//		~WayPoint() { fmt::print("§{} {} ", "WayPoint::~WayPoint()", ct); _ctrsgn(typeid(*this), true); }

		template<CoordType type>
		void convert() const;
		void validate(bool = true) const;
		template<CoordType type>
		vector<string> format_ct() const;
		vector<string> format(bool usenames) const;
};


/// __________________________________________________
/// __________________________________________________
/// Validation
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

/// __________________________________________________
/// __________________________________________________
/// Exported functions
NumericVector as_coords(NumericVector, const int);
NumericVector convertcoords(NumericVector, const int);
NumericVector latlon(NumericVector, LogicalVector);
NumericVector validatecoords(NumericVector, const bool);
CharacterVector formatcoords(NumericVector, bool, bool);
DataFrame as_waypointsdefault(DataFrame, const int);
DataFrame convertwaypoints(DataFrame, const int);
DataFrame validatewaypoints(DataFrame, const bool);
CharacterVector formatwaypoints(DataFrame, bool, bool);
CharacterVector ll_headers(int, const int);
NumericVector as_coordswaypoints(DataFrame, bool);


#endif  // COORDBASE_H_
