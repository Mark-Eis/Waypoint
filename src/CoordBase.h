/// __________________________________________________
/// CoordBase.h
/// __________________________________________________

#ifndef COORDBASE_H_
#define COORDBASE_H_

#include <concepts>

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
	explicit VecTypeBase( vector<T>::size_type count ) : vector<T>(count) {}			// ≈ "default"
	VecTypeBase(const VecTypeBase&) = delete;											// copy constructor
	VecTypeBase(const vector<T>& vt) : vector<T>{ vt } {}								// copy constructor

	VecTypeBase& operator=(const VecTypeBase&) = delete;									// copy assignment
	VecTypeBase& operator=(const vector<T>& vt)											// copy assignment
	{
		vector<T>::operator=(vt);
		return *this;
	}

	VecTypeBase(VecTypeBase&&) = default;												// move constructor

	VecTypeBase(vector<T>&& vt) : vector<T>{ std::move(vt) } {}							// move constructor

	VecTypeBase& operator=(VecTypeBase&&) = default;									// move assignment

	VecTypeBase& operator=(vector<T>&& vt)												// move assignment
	{
		vector<T>::operator=(std::move(vt));
		return *this;
	}

	virtual ~VecTypeBase() = 0;
};

template<typename T>
VecTypeBase<T>::~VecTypeBase() {}


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
/// Concept —— vectype
template <typename T>
concept vectype = 
	DVecType<T> || SVecType<T>;

/// __________________________________________________
/// __________________________________________________
/// Formula simplification
inline double mod1by60(double);
inline double mod1e2(double);
inline double round2(double, int = 2);
inline double polish(double);

/// __________________________________________________
/// __________________________________________________
/// Utility
template<typename U> 
inline vector<U> get_vec_attr(const NumVec_or_DataFrame auto&, const string);
inline int get_fmt_attribute(const NumVec_or_DataFrame auto&);
int check_logical_attr(NumVec_or_DataFrame auto, const string);
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
	int get_deg(double x) const { return int(x); }
	double get_decdeg(double x) const { return x; }
	int get_min(double x) const { return int(get_decmin(x)); }
	double get_decmin(double x) const { return polish(mod1by60(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
};

/// __________________________________________________
/// Specialised struct for degrees and minutes
template<>
struct FamousFive<DegMinVecDouble> {
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
	double operator()(double n) const { return ff.get_decdeg(n); }
};

/// __________________________________________________
/// Specialised struct for degrees and minutes
template<DVecType T>
struct Convertidor<T, DegMinVecDouble>{
	FamousFive<T> ff {};
	double operator()(double n) const { return ff.get_deg(n) * 1e2 + ff.get_decmin(n); }
};

/// __________________________________________________
/// Specialised struct for degrees, minutes and seconds
template<DVecType T>
struct Convertidor<T, DegMinSecVecDouble>{
	FamousFive<T> ff {};
	double operator()(double n) const { return ff.get_deg(n) * 1e4 + ff.get_min(n) * 1e2 + ff.get_sec(n); }
};

/// __________________________________________________
/// Template aliases
template<DVecType T>
using ConvertidorDecDegVec = Convertidor<DecDegVecDouble, T>;
template<DVecType T>
using ConvertidorDegMinVec = Convertidor<DegMinVecDouble, T>;
template<DVecType T>
using ConvertidorDegMinSecVec = Convertidor<DegMinSecVecDouble, T>;


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
	ostringstream ostrstr;
	string operator()(double n)
	{
		ostrstr.str(""s);
		ostrstr << setw(11) << setfill(' ')  << fixed << setprecision(6) << ff.get_decdeg(n) << "\u00B0";
		return ostrstr.str();
	}
};

/// __________________________________________________
/// Specialised struct for degrees and minutes
template<DVecType T>
struct Formateador<T, DegMinVecString>{
	FamousFive<T> ff {};
	ostringstream ostrstr;
	string operator()(double n)
	{
		auto deg {abs(ff.get_deg(n))};
		auto min {fabs(ff.get_decmin(n))};
		if (round2(min) > 59.99995) {
			++deg;
			min = 0;
		} 
		ostrstr.str(""s);
		ostrstr << setw(3) << setfill(' ') << deg << "\u00B0"
				<< setw(7) << setfill('0') << fixed << setprecision(4) << min << "\u2032";
		return ostrstr.str();
	}
};

/// __________________________________________________
/// Specialised struct for degrees, minutes and seconds
template<DVecType T>
struct Formateador<T, DegMinSecVecString>{
	FamousFive<T> ff {};
	ostringstream ostrstr;
	string operator()(double n)
	{
		auto min {abs(ff.get_min(n))};
		auto sec {fabs(ff.get_sec(n))};
		if (round2(sec) > 59.995) {
			++min;
			sec = 0;
		} 
		ostrstr.str(""s);
		ostrstr << setw(3) << setfill(' ') << abs(ff.get_deg(n)) << "\u00B0"
				<< setw(2) << setfill('0') << min << "\u2032"
				<< setw(5) << fixed << setprecision(2) << sec << "\u2033";
		return ostrstr.str();
	}
};

/// __________________________________________________
/// Template aliases
template<SVecType T>
using FormateadorDecDegVec = Formateador<DecDegVecDouble, T>;
template<SVecType T>
using FormateadorDegMinVec = Formateador<DegMinVecDouble, T>;
template<SVecType T>
using FormateadorDegMinSecVec = Formateador<DegMinSecVecDouble, T>;


/// __________________________________________________
/// Concept —— float_or_string
template<typename T>
concept float_or_string =
	std::floating_point<T> || std::same_as<T, string>;

/// __________________________________________________
/// Concept —— functador
template<typename T>
concept functador =
	requires (T t, double n) {
		{ t.operator()(n) } -> float_or_string;
	};

/// __________________________________________________
/// __________________________________________________
template<DVecType T, typename S>
class Coords;

/// Concept —— sufijo
template<typename T>
concept sufijo =
	requires (T t, DecDegVecString& vt) {
		{ t.suffix(vt) };
	} ||
	requires (T t, DegMinVecString& vt) {
		{ t.suffix(vt) };
	} ||
	requires (T t, DegMinSecVecString& vt) {
		{ t.suffix(vt) };
	};

/// __________________________________________________
/// Concept —— coords_t
template<typename T>
concept coords_t =
	requires (T t, CoordType ct, DecDegVecString& vt) {
		{ t.validate() } -> std::same_as<const vector<bool>>;
		{ t.add_suffix(vt) };
	} ||
	requires (T t, CoordType ct, DegMinVecString& vt) {
		{ t.validate() } -> std::same_as<const vector<bool>>;
		{ t.add_suffix(vt) };
	} ||
	requires (T t, CoordType ct, DegMinSecVecString& vt) {
		{ t.validate() } -> std::same_as<const vector<bool>>;
		{ t.add_suffix(vt) };
	};


/// __________________________________________________
/// Coords class
template<DVecType T, typename S>
class Coords {
	protected:
		T dv;
		const vector<bool> latlon;

		template<vectype U, functador V>
		inline U conform0() const;
	public:
		explicit Coords(NumericVector);
		Coords(const Coords&) = delete;								// Disallow copying
		Coords& operator=(const Coords&) = delete;					//  ——— ditto ———
		Coords(Coords&&) = delete;									// Disallow transfer ownership
		Coords& operator=(Coords&&) = delete;						// Disallow moving
		virtual ~Coords() = 0;

        	template<typename U, template <typename V> typename F>
		vector<U> conform(CoordType) const;							// Non-const return type avoids making unnecessary copy
		const vector<bool> validate() const;
		void add_suffix(vectype auto&) const;
};


/// __________________________________________________
/// SufijoCoords class
template<DVecType T>
class SufijoCoords final : public Coords<T, SufijoCoords<T>> {
	public:
		using Coords<T, SufijoCoords<T>>::Coords;
		using Coords<T, SufijoCoords<T>>::latlon;
		using Coords<T, SufijoCoords<T>>::dv;
		void suffix(vectype auto&) const;
};


/// __________________________________________________
/// SufijoWaypoints class
template<DVecType T>
class SufijoWaypoints final : public Coords<T, SufijoWaypoints<T>> {
	public:
		using Coords<T, SufijoWaypoints<T>>::Coords;
		using Coords<T, SufijoWaypoints<T>>::latlon;
		using Coords<T, SufijoWaypoints<T>>::dv;
		void suffix(vectype auto&) const;
};


/// __________________________________________________
/// __________________________________________________
/// Switches for Coords<DVecType>
vector<double> convert_switch(const NumericVector, CoordType); 
vector<string> format_switch_c(const NumericVector, CoordType); 
vector<string> format_switch_w(const NumericVector, CoordType); 
const vector<bool> validate_switch(const NumericVector); 

/// __________________________________________________
/// __________________________________________________
/// Type aliases
template<typename T>
using bisvec = array<vector<T>, 2>;
template<typename T>
using bisconstvec = array<const vector<T>, 2>;

/// __________________________________________________
/// __________________________________________________
/// Waypoints class
class Waypoints {
		NumericVector nv_lat;
		NumericVector nv_lon;
		const vector<int> llcols;
	public:
		Waypoints(const DataFrame&);
		Waypoints(const Waypoints&) = delete;					// Disallow copying
		Waypoints& operator=(const Waypoints&) = delete;			//  ——— ditto ———
		Waypoints(Waypoints&&) = delete;						// Disallow transfer ownership
		Waypoints& operator=(Waypoints&&) = delete;				// Disallow moving
		~Waypoints();

		vector<double> convert(CoordType, bool) const;
		vector<string> format(CoordType, bool) const;
		const vector<bool> validate(bool) const;
};

/// __________________________________________________
/// __________________________________________________
/// for Waypoints
inline const bisconstvec<bool> validate_switch(const DataFrame);

/// __________________________________________________
/// __________________________________________________
/// Validation
bool check_valid(const NumericVector, bool = false);
bool check_valid(const DataFrame, bool = false);
bool validate(const NumVec_or_DataFrame auto, bool = false);
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