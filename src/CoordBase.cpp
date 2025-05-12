#include <Rcpp.h>
#include <cxxabi.h>
using namespace Rcpp;

using std::vector; 
using std::string;
using std::transform;
using std::ostream;
using std::cout;	// Deprecate ?
using std::endl;	// Deprecate ?

#include "CoordBase.h"

#define FMT_HEADER_ONLY
#include </opt/homebrew/Cellar/fmt/11.1.4/include/fmt/format.h>		// verbose path needs sorting!
// #include <fmt/format.h>		// verbose path needs sorting!

/// __________________________________________________
/// __________________________________________________
/// Development and Debugging functions

/// Report object construction and destruction
void _ctrsgn(const std::type_info& obj, bool destruct = false)
{
	fmt::print("{}ing ", destruct ? "Destroy" : "Construct");
	std::fflush(nullptr);
	string s = obj.name();
	system(("c++filt -t " + s).data());
}

/// Demangle object names functor
class Demangler {
	char* p;
	int status = 0;
public:
	Demangler(const std::type_info& obj) : p(abi::__cxa_demangle(obj.name(), NULL, NULL, &status)) {}
	~Demangler() { std::free(p); }
	operator string() const { return string("\"") + p + "\" (status " + std::to_string(status) + ")"; }
};

ostream& operator<< (ostream& stream, const Demangler& d)
{
//  cout << "ostream& operator<< (ostream&, const Demangler&) ";
  return stream << string(d);
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
/// Utility functions

/// __________________________________________________
/// Return named attribute as vector<U> or empty vector<U>
template<class T, class U> 
inline vector<U> get_vec_attr(const T& t, const char* attrname)
{
//	cout << "@get_vec_attr<T, U>(const T&, const char*) attr \"" << attrname << "\" " << boolalpha << t.hasAttribute(attrname) << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	return t.hasAttribute(attrname) ? as<vector<U>>(t.attr(attrname)) : vector<U>();
}


/// __________________________________________________
/// Return "fmt" attribute as int
template<class T>
inline int get_fmt_attribute(const T& t)
{
	fmt::print("@{} fmt={}\n", "get_fmt_attribute<T>(const T&)", as<int>(t.attr("fmt")));
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	return as<int>(t.attr("fmt"));
}


/// __________________________________________________
/// Does object inherit given class?
template<class T>
inline void checkinherits(T& t, const char* classname)
{
//	cout << "@checkinherits<T>(T& t, const char* classname) t " << Demangler(typeid(t))  << " classname \"" << classname << "\"" << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	if (!t.inherits(classname)) stop("Argument must be a \"%s\" object", classname);
}


/// __________________________________________________
/// Is item number present in object? (Using C++ numbering)
template<class T>
inline bool is_item_in_obj(const T t, const int item)
{
//	cout << "@is_item_in_obj(T, int)\n";
	if (NA_INTEGER == item)
		return false;
	else
		return !(item < 0) && item < t.size();
}


/// __________________________________________________
/// Standarise width of strings in vector to that of the longest
inline void stdlenstr(vector<string>& sv)
{
	fmt::print("@{}\n", "stdlenstr(vector<string>&)");
	int maxwdth = max_element(sv.begin(), sv.end(), [](const string& a, const string& b){ return a.size() < b.size(); })->size();
	transform(sv.begin(), sv.end(), sv.begin(), [maxwdth](const string& s) { return fmt::format("{:<{}}", s, maxwdth); });
}


/// __________________________________________________
/// Prefix vector<string> elements with elements of vector<T>—default for vector<string> prefix
template<class T>
inline void prefixvecstr(vector<string>& sv, const vector<T>& prefix)
{
	fmt::print("@{}\n", "prefixvecstr<T>(vector<string>&, const vector<T>&)");
	transform(sv.begin(), sv.end(), prefix.begin(), sv.begin(), [](string& lls, const string& name) { return name + "  " + lls; });	
}


/// __________________________________________________
/// Specialisation for vector<int> prefix
template<>
inline void prefixvecstr(vector<string>& sv, const vector<int>& prefix)
{
//	cout << "@prefixvecstr<>(vector<string>&, const vector<int>&)\n";
	transform(sv.begin(), sv.end(), prefix.begin(), sv.begin(), [](string& lls, const int name) { return std::to_string(name) + "  " + lls; });	
}


/// __________________________________________________
/// Prefix vector<string> elements with elements of RObject 
inline bool prefixwithnames(vector<string>& sv, RObject& namesobj)
{
//	cout << "@prefixwithnames(vector<string>&, RObject&)\n";
	if (is<CharacterVector>(namesobj)) {
		vector<string>&& names = as<vector<string>>(namesobj);
		stdlenstr(names);
		prefixvecstr(sv, names);
	} else if(is<IntegerVector>(namesobj))
		prefixvecstr(sv, as<vector<int>>(namesobj));
	else
		return false;
	return true;
}


/// __________________________________________________
/// string to lower case (see cppreference.com std::tolower)
inline string str_tolower(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return tolower(c); });
    return s;
}


/// __________________________________________________
/// Find position of name within object names
template<class T>
int nameinobj(const T t, const char* name)
{
//	cout << "@nameinobj<T>(const T, const char*) name is " << name << endl;
	static_assert(std::is_same<List, T>::value || std::is_same<DataFrame, T>::value, "T must be List or DataFrame");
	vector<string> names { get_vec_attr<T, string>(t, "names") };
	if (!names.size())
		return -1;
	typedef decltype(names.size()) Tmp;
	Tmp i = 0;
	for (auto str : names ) {
//		cout << "@nameinobj<T>(const T, const char*) testing " << str << endl;
		if (!str_tolower(str).compare(name)) {
//			cout << "@nameinobj<T>(const T, const char*) found " << str << endl;
			break;
		}
		i++;
	}
	if (i == names.size())
		i = -1;
	return i;
}


/// __________________________________________________
/// Retrieve names column or row.names from DataFrame as Robject
RObject getnames(const DataFrame df)
{
//	cout << "@getnames(const DataFrame)\n";
	vector<int> namescolvec { get_vec_attr<DataFrame, int>(df, "namescol") };
	if (1 == namescolvec.size()) {
		int namescol = namescolvec[0] - 1;
		if (is_item_in_obj(df, namescol))
			return df[namescol];
		else
			stop("Invalid \"namescol\" attribute! (item not in object)");
	} else
		if (df.hasAttribute("row.names"))
			return df.attr("row.names");
		else
			stop("Missing row.names!");
}


/// __________________________________________________
/// __________________________________________________
/// CoordType enum class

/// __________________________________________________
/// Formatter struct template specialisation
 
// auto fmt::formatter<CoordType_enum::CoordType>::format(CoordType_enum::CoordType ct, format_context& ctx) const
auto fmt::formatter<CoordType>::format(CoordType ct, format_context& ctx) const
	-> format_context::iterator
{
//	using enum CoordType_enum::CoordType;  // see: using enum (https://en.cppreference.com/w/cpp/language/enum#Using-enum-declaration)
	using enum CoordType;
	string_view name = "unknown";
	switch (ct) {
		case decdeg:
			name = "DecDeg";
			break;

		case degmin:
			name = "DegMin";
			break;

		case degminsec:
			name = "DegMinSec";
			break;

			default:
//				stop("fmt::formatter<CoordType_enum::CoordType>::format(CoordType_enum::CoordType, format_context&) my bad");
				stop("fmt::formatter<CoordType>::format(CoordType, format_context&) my bad");
	}
	return formatter<string_view>::format(name, ctx);
}


/// __________________________________________________
/// Convert int to CoordType enum
inline const CoordType get_coordtype(const int i)
{
//	cout << "@get_coordtype(int) " << i << endl;
	fmt::print("@{} {}\n", "get_coordtype(int)" , i);
	if (i < 1 || i > 3)
		stop("\"fmt\" must be between 1 and 3");
	return vector<CoordType>{ CoordType::decdeg, CoordType::degmin, CoordType::degminsec }[i - 1];
}


/// __________________________________________________
/// Convert "fmt" attribute to CoordType enum
template<class T>
inline const CoordType get_coordtype(const T& t)
{
//	cout << "@get_coordtype<T>(const T&) " << get_fmt_attribute(t) << endl;
	fmt::print("@{} fmt={}\n", "get_coordtype<T>(const T&)", get_fmt_attribute(t));
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	return get_coordtype(get_fmt_attribute(t));
}


/// __________________________________________________
/// Convert CoordType enum to int
inline const int coordtype_to_int(CoordType ct)
{
//	cout << "@coordtype_to_int(CoordType ct) " << static_cast<char>(ct) + 1 << endl;
	fmt::print("@{} ct={}\n", "coordtype_to_int(CoordType)", ct);
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
/// FamousFive Class and Derived Classes

struct FamousFive {
	FamousFive() { fmt::print("§{} ", "FamousFive()"); _ctrsgn(typeid(*this)); }
	virtual ~FamousFive() = 0;	
	virtual int get_deg(double x) const = 0;
	virtual double get_decdeg(double x) const = 0;
	virtual int get_min(double x) const = 0;
	virtual double get_decmin(double x) const = 0;
	virtual double get_sec(double x) const = 0;
};

inline FamousFive::~FamousFive()
{
	fmt::print("§{} ", "~FamousFive()"); _ctrsgn(typeid(*this), true);
	std::fflush(nullptr);
}	

/// __________________________________________________
/// Derived class for decimal degrees	
struct FF_decdeg : public FamousFive {
	FF_decdeg() { fmt::print("§{} ", "FF_decdeg()"); _ctrsgn(typeid(*this)); }
//	~FF_decdeg() = default;
	~FF_decdeg() { fmt::print("§{} ", "~FF_decdeg()"); _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x); }
	double get_decdeg(double x) const { return x; }
	int get_min(double x) const { return (int(x * 1e6) % int(1e6)) * 6e-5; }
	double get_decmin(double x) const { return polish(mod1by60(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
} ff_decdeg;

/// __________________________________________________
/// Derived class for degrees and minutes
struct FF_degmin : public FamousFive {
	FF_degmin() { fmt::print("§{} ", "FF_degmin()"); _ctrsgn(typeid(*this)); }
//	~FF_degmin() = default;
	~FF_degmin() { fmt::print("§{} ", "~FF_degmin()"); _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x / 1e2); }
	double get_decdeg(double x) const { return int(x / 1e2) + mod1e2(x) / 60; }
	int get_min(double x) const { return int(x) % int(1e2); }
	double get_decmin(double x) const { return polish(mod1e2(x)); }
	double get_sec(double x) const { return mod1by60(get_decmin(x)); }
} ff_degmin;

/// __________________________________________________
/// Derived class for degrees, minutes and seconds
struct FF_degminsec : public FamousFive {
	FF_degminsec() { fmt::print("§{} ", "FF_degminsec()"); _ctrsgn(typeid(*this)); }
//	~FF_degminsec() = default;
	~FF_degminsec() { fmt::print("§{} ", "~FF_degminsec()"); _ctrsgn(typeid(*this), true); }
	int get_deg(double x) const { return int(x / 1e4); }
	double get_decdeg(double x) const { return int(x / 1e4) + (double)int(fmod(x, 1e4) / 1e2) / 60 + mod1e2(x) / 3600; }
	int get_min(double x) const { return (int(x) % int(1e4)) / 1e2; }
	double get_decmin(double x) const { return int(fmod(x, 1e4) / 1e2) + mod1e2(x) / 60; }
	double get_sec(double x) const { return mod1e2(x); }
} ff_degminsec;

vector<FamousFive*> vff { &ff_decdeg, &ff_degmin, &ff_degminsec };


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
//			cout << "§Convertor<CoordType>::Convertor(const FamousFive&) "; _ctrsgn(typeid(*this));
		}
		~Convertor() = default;
//		~Convertor() { cout << "§Convertor<type>::~Convertor() "; _ctrsgn(typeid(*this), true); }
		double operator()(double n);
};


/// __________________________________________________
/// Default operator(), for decimal degrees
template<CoordType type>
inline double Convertor<type>::operator()(double n)
{
//	cout << "@Convertor<CoordType>::operator() [default for CoordType::decdeg]\n";
	return ff.get_decdeg(n);
}


/// __________________________________________________
/// Specialised operator() for degrees and minutes
template<>
inline double Convertor<CoordType::degmin>::operator()(double n)
{
//	cout << "@Convertor<CoordType::degmin>::operator()\n";
	return ff.get_deg(n) * 1e2 + ff.get_decmin(n);
}


/// __________________________________________________
/// Specialised operator() for degrees, minutes and seconds
template<>
inline double Convertor<CoordType::degminsec>::operator()(double n)
{
//	cout << "@Convertor<CoordType::degminsec>::operator()\n";
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
			fmt::print("§Format<CoordType::{}>::Format()(const FamousFive&) ", type); _ctrsgn(typeid(*this));
			std::fflush(nullptr);
		}
//		~Format() = default;
		~Format() { fmt::print("§Format<CoordType::{}>::~Format() ", type); _ctrsgn(typeid(*this), true); }
		string operator()(double n) const;
};


/// __________________________________________________
/// Default operator(), for decimal degrees
template<CoordType type>
inline string Format<type>::operator()(double n) const
{
	fmt::print("@Format<CoordType::{}>::operator() [default]\n", type);
	return fmt::format("{:>{}.{}f}\u00B0", ff.get_decdeg(n), 11, 6);
}

/// __________________________________________________
/// Specialised operator() for degrees and minutes
template<>
inline string Format<CoordType::degmin>::operator()(double n) const
{
	fmt::print("@Format<CoordType::{}>::operator()\n", CoordType::degmin);
	return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) +
		   fmt::format("{:0>{}.{}f}\u2032", abs(ff.get_decmin(n)), 7, 4);
}

/// __________________________________________________
/// Specialised operator() for degrees, minutes and seconds
template<>
inline string Format<CoordType::degminsec>::operator()(double n) const
{
	fmt::print("@Format<CoordType::{}>::operator()\n", CoordType::degminsec);
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
			fmt::print("§FormatLL<{}, CoordType::{}>::FormatLL(const FamousFive&, vector<bool>&) ", "Coord or WayPoint", type); _ctrsgn(typeid(*this));
			std::fflush(nullptr);
		}
		~FormatLL() { fmt::print("§FormatLL<{}, CoordType::{}>::~FormatLL() ", "Coord or WayPoint", type); _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
			fmt::print("@FormatLL<{}, CoordType::{}>::operator(string, double) [default]\n", "Coord or WayPoint", type);
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
			fmt::print("§FormatLL<{}, CoordType::{}>::FormatLL(const FamousFive&, vector<bool>&) ", "Coord", CoordType::decdeg ); _ctrsgn(typeid(*this));
			std::fflush(nullptr);
		}
//		~FormatLL() = default;
		~FormatLL() { fmt::print("§FormatLL<{}, CoordType::{}>::~FormatLL() ", "Coord", CoordType::decdeg); _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
			fmt::print("@FormatLL<{}, CoordType::{}>::operator(string, double)\n", "Coord", CoordType::decdeg);
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
			fmt::print("§FormatLL<{}, CoordType::{}>::FormatLL(const FamousFive&, vector<bool>&) ", "WayPoint", CoordType::decdeg); _ctrsgn(typeid(*this));
			std::fflush(nullptr);
		}
//		~FormatLL() = default;
		~FormatLL() { fmt::print("§FormatLL<{}, CoordType::{}>::~FormatLL() ", "WayPoint", CoordType::decdeg); _ctrsgn(typeid(*this), true); }
		string operator()(string ostr, double n)
		{
			fmt::print("@FormatLL<{}, CoordType::{}>::operator(string, double)\n", "WayPoint", CoordType::decdeg);
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
//			cout << "§Validator::Validator(const FamousFive&, vector<bool>&) "; _ctrsgn(typeid(*this));
		}
		~Validator() = default;
//		~Validator() { cout << "§Validator::~Validator() "; _ctrsgn(typeid(*this), true); }
		bool operator()(double n)
		{
//			cout << "@Validator() " << " validating: " << std::setw(9) << std::setfill(' ') << n << endl;
			return !((abs(ff.get_decdeg(n)) > (ll_size && (ll_size > 1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(abs(ff.get_decmin(n)) >= 60) ||
				(abs(ff.get_sec(n)) >= 60));
		}
};


/// __________________________________________________
/// __________________________________________________
/// CoordType switches

/// __________________________________________________
/// Convert coords or waypoints format CoordType switch 
template<class T, class U>
void convert_switch(T t, CoordType newtype)
{
//	cout << "@convert_switch<T&, U>(T, CoordType) t " << Demangler(typeid(t)) << " newtype " << coordtype_to_int(newtype) << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	CoordType type = get_coordtype(t);
	U u(type, t);
	u.validate();

	if (type != newtype) {
		switch (newtype)
		{
			case CoordType::decdeg:
				u.template convert<CoordType::decdeg>();
				break;

			case CoordType::degmin:
				u.template convert<CoordType::degmin>();
				break;

			case CoordType::degminsec:
				u.template convert<CoordType::degminsec>();
				break;

			default:
				stop("convert_switch<T&, U>(const T&, U) my bad");
		}
		t.attr("fmt") = coordtype_to_int(newtype) + 1;
	}
}


/// __________________________________________________
/// Format coords or waypoints vector<string> CoordType switch 
template<class T>
vector<string> format_switch(const T& t)
{
//	cout << "@format_switch<T>(const T&, CoordType) " << Demangler(typeid(t)) << " ct " << coordtype_to_int(ct) << endl;
	static_assert(std::is_same<Coord, T>::value || std::is_same<WayPoint, T>::value, "T must be Coord or WayPoint");
	switch (t.get_coordtype())
	{
		case CoordType::decdeg:
			return t.template format_ct<CoordType::decdeg>();

		case CoordType::degmin:
			return t.template format_ct<CoordType::degmin>();

		case CoordType::degminsec:
			return t.template format_ct<CoordType::degminsec>();

		default:
			stop("format_switch(const T&, CoordType) my bad");
	}
}


/// __________________________________________________
/// __________________________________________________
/// Coordbase class

Coordbase::Coordbase(CoordType _ct) :
	ct(_ct), ff(*vff[coordtype_to_int(ct)])
{
//	cout << "§Coordbase::Coordbase(CoordType) "; _ctrsgn(typeid(*this));
}


Coordbase::~Coordbase()
{
//	cout << "§Coordbase::~Coordbase() "; _ctrsgn(typeid(*this), true);
}


CoordType Coordbase::get_coordtype() const
{
//	cout << "@Coordbase::get_coordtype() ct " << coordtype_to_int(ct) << endl;
	fmt::print("@{} ct={}\n", "Coordbase::get_coordtype()", coordtype_to_int(ct));
	return ct;
}


/// __________________________________________________
/// Coordinate derived class

Coord::Coord(CoordType ct, const NumericVector nv) :
	Coordbase(ct), nv(nv),
	latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") } //,
{
//	cout << "§Coord::Coord(CoordType, const NumericVector) "; _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Convert NumericVector coordinate format
template<CoordType newtype>
inline void Coord::convert() const
{
//	cout << "@Coord::convert<CoordType>() newtype " << coordtype_to_int(newtype) + 1 << endl;
	transform(nv.begin(), nv.end(), const_cast<NumericVector&>(nv).begin(), Convertor<newtype>(ff));
}


/// __________________________________________________
/// Validate coords vector
void Coord::validate(bool warn) const
{
//	cout << "@Coord::validate() " << Demangler(typeid(*this)) << " latlon " << LogicalVector(wrap(latlon)) << endl;
	vector<bool>& non_const_valid { const_cast<vector<bool>&>(valid) };
	non_const_valid.assign(nv.size(), {false});
	transform(nv.begin(), nv.end(), non_const_valid.begin(), Validator(ff, latlon));
	if (all_of(valid.begin(), valid.end(), [](bool v) { return v;}))
		non_const_valid.assign({true});
	else
		if (warn)
			warning("Validation failed!");
	const_cast<NumericVector&>(nv).attr("valid") = valid;
}


/// __________________________________________________
/// Format coordinates as vector<string> of CoordType
template<CoordType type>
vector<string> Coord::format_ct() const
{
//	cout << "@Coord::format_ct<CoordType>() " << Demangler(typeid(*this)) << endl;
	vector<string> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), Format<type>(ff));
	transform(out.begin(), out.end(), nv.begin(), out.begin(), FormatLL<Coord, type>(ff, latlon));
	return out;
}


/// __________________________________________________
/// Format coords vector<string> with names
vector<string> Coord::format(bool usenames) const
{
//	cout << "@Coord::format(bool) " << Demangler(typeid(*this)) << endl;
	vector<string>&& sv = format_switch(*this);
	vector<string> names { get_vec_attr<NumericVector, string>(nv, "names") };
	if (names.size() && usenames) {
		stdlenstr(names);
		prefixvecstr(sv, names);
	}
	return sv;
}


/// __________________________________________________
/// __________________________________________________
/// Waypoint class

WayPoint::WayPoint(CoordType ct, const DataFrame df) :
	Coordbase(ct), df(df),
	nvlat(df[get_vec_attr<DataFrame, int>(df, "llcols")[0] - 1]), 
	nvlon(df[get_vec_attr<DataFrame, int>(df, "llcols")[1] - 1])
{
//	cout << "§WayPoint::WayPoint(CoordType ct, const DataFrame) "; _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Convert DataFrame coordinate format
template<CoordType newtype>
inline void WayPoint::convert() const
{
// 	cout << "@WayPoint::convert<CoordType>() newtype " << coordtype_to_int(newtype) + 1 << endl;
	transform(nvlat.begin(), nvlat.end(), const_cast<NumericVector&>(nvlat).begin(), Convertor<newtype>(ff));
	transform(nvlon.begin(), nvlon.end(), const_cast<NumericVector&>(nvlon).begin(), Convertor<newtype>(ff));
}


/// __________________________________________________
/// Validate WayPoint
void WayPoint::validate(bool warn) const
{
//	cout << "@WayPoint::validate(bool) " << Demangler(typeid(*this)) << endl;

	vector<bool>& non_const_validlat { const_cast<vector<bool>&>(validlat) };
	non_const_validlat.assign(nvlat.size(), {false});
	transform(nvlat.begin(), nvlat.end(), non_const_validlat.begin(), Validator(ff, vector<bool>{ true }));

	vector<bool>& non_const_validlon { const_cast<vector<bool>&>(validlon) };
	non_const_validlon.assign(nvlon.size(), {false});
	transform(nvlon.begin(), nvlon.end(), non_const_validlon.begin(), Validator(ff, vector<bool>{ false }));

	if (all_of(validlat.begin(), validlat.end(), [](bool v) { return v;}))
		non_const_validlat.assign({true});
	else
		if (warn)
			warning("Validation of latitude failed!");
	const_cast<DataFrame&>(df).attr("validlat") = validlat;

	if (all_of(validlon.begin(), validlon.end(), [](bool v) { return v;}))
		non_const_validlon.assign({true});
	else
		if (warn)
			warning("Validation of longitude failed!");
	const_cast<DataFrame&>(df).attr("validlon") = validlon;
}


/// __________________________________________________
/// Format waypoints as vector<string> of CoordType
template<CoordType type>
vector<string> WayPoint::format_ct() const
{
//	cout << "@WayPoint::format_ct() " << Demangler(typeid(*this)) << endl;
	vector<string> sv_lat(nvlat.size());
	transform(nvlat.begin(), nvlat.end(), sv_lat.begin(), Format<type>(ff));
	transform(sv_lat.begin(), sv_lat.end(), nvlat.begin(), sv_lat.begin(), FormatLL<WayPoint, type>(ff, vector<bool>{ true }));
	vector<string> sv_lon(nvlon.size());
	transform(nvlon.begin(), nvlon.end(), sv_lon.begin(), Format<type>(ff));
	transform(sv_lon.begin(), sv_lon.end(), nvlon.begin(), sv_lon.begin(), FormatLL<WayPoint, type>(ff, vector<bool>{ false }));

	vector<string> out(sv_lat.size());
	transform(sv_lat.begin(), sv_lat.end(), sv_lon.begin(), out.begin(), [](string& latstr, string& lonstr) { return latstr + "  " + lonstr; });
	return out;
}


/// __________________________________________________
/// Format waypoints vector<string> with names
vector<string> WayPoint::format(bool usenames) const
{
//	cout << "@WayPoint::format(bool) " << Demangler(typeid(*this)) << endl;
	vector<string>&& sv = format_switch(*this);
	if (usenames) {
		RObject names = getnames(df);
		if (!prefixwithnames(sv, names))
			stop("Invalid \"namescol\" attribute!");
	}
	return sv;
}


/// __________________________________________________
/// __________________________________________________
/// Validation functions

/// __________________________________________________
/// Check "valid" attribute of NumericVector all true

bool check_valid(const NumericVector nv)
{
//	cout << "@check_valid(const NumericVector)" << endl;
	bool unvalidated = false;
	bool valid = validated(nv, "valid", unvalidated);
	if (unvalidated)
		revalid_Coord(nv);
	return valid;
}


/// __________________________________________________
/// Check "lat_valid" and "lon_valid attributes of DataFrame are all true
bool check_valid(const DataFrame df)
{
//	cout << "@check_valid(const DataFrame)\n";
	bool unvalidated = false;

	bool latvalid = validated(df, "validlat", unvalidated);
	if (unvalidated)
		return revalid_WayPoint(df);
	bool lonvalid = validated(df, "validlon", unvalidated);
	if (unvalidated)
		return revalid_WayPoint(df);

	if (!latvalid)
		warning("Invalid latitude!");
	if (!lonvalid)
		warning("Invalid longitude!");
	return latvalid || lonvalid;
}


/// __________________________________________________
/// Check NumericVector or DataFrame has been validated and valid vector attribute all true
template<class T>
bool validated(T t, const char* attrname, bool& unvalidated)
{
//	cout << "@validated<T>(T, const char*, bool&)" << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	const vector<bool>&& validvec = get_vec_attr<T, bool>(t, attrname);
	bool valid = all_of(validvec.begin(), validvec.end(), [](bool v) { return v;});
	unvalidated = (validvec.size()) ? false : true;
	return valid;
}


/// __________________________________________________
/// Revalidate NumericVector or DataFrame
template<class T, class U>
const T revalidate(const T t)
{
//	cout << "@revalidate<T, U>(const T) t" << Demangler(typeid(t))  << endl;
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	warning("Revalidating %s…!", Demangler(typeid(t)));
	validate<T, U>(t);	
	return check_valid(t);
}


/// __________________________________________________
/// Validate NumericVector or DataFrame
template<class T, class U>
inline const T validate(const T t)
{
//	cout << "@validate<T, U>(const T)\n";
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	U(get_coordtype(t), t).validate();
	return t;	
}


/// __________________________________________________
/// Check df has valid "llcols" attribute
bool valid_ll(const DataFrame df)
{
//	cout << "@valid_ll(const DataFrame)\n";
	bool valid = false;
	vector<int> llcols { get_vec_attr<DataFrame, int>(df, "llcols") };
	if (2 == llcols.size()) {
		transform(llcols.begin(), llcols.end(), llcols.begin(), [](int x){ return --x; });
		if (is_item_in_obj(df, llcols[0]) && is_item_in_obj(df, llcols[1]) && llcols[0] != llcols[1])
			if (is<NumericVector>(df[llcols[0]]) && is<NumericVector>(df[llcols[1]]))
				valid = true;
	}
	return valid;
}


/// __________________________________________________
/// __________________________________________________
/// Exported functions

/// __________________________________________________
/// Create coords
//' @rdname coords 
// [[Rcpp::export(name = "as_coords.default")]]
NumericVector as_coords(NumericVector object, const int fmt = 1)
{
//	cout << "——Rcpp::export——coords(NumericVector)\n";
	object.attr("fmt") = fmt;
	convert_switch<NumericVector, Coord>(object, get_coordtype(fmt));
	object.attr("class") = "coords";
	return object;
}


/// __________________________________________________
/// Convert coords format
//' @rdname convert
// [[Rcpp::export(name = "convert.coords")]]
NumericVector convertcoords(NumericVector x, const int fmt)
{
//	cout << "——Rcpp::export——convertcoords(NumericVector, const int) from " << get_fmt_attribute(x) << " to " << fmt << endl;
	checkinherits(x, "coords");
	CoordType type = get_coordtype(x);
	CoordType newtype = get_coordtype(fmt);
	if (newtype == type) {
//		cout << "——fmt out == fmt in!——" << endl;
		if (!check_valid(x))
			stop("Invalid coords!");
	} else 
		convert_switch<NumericVector, Coord>(x, newtype);
	return x;
}


/// __________________________________________________
/// Set latlon attribute on "coords" NumericVector and revalidate
//' @rdname coords
// [[Rcpp::export(name = "`latlon<-`")]]
NumericVector latlon(NumericVector cd, LogicalVector value)
{
//	cout << "——Rcpp::export——latlon(NumericVector, LogicalVector)\n";
	checkinherits(cd, "coords");
	if (value.size() != cd.size() && value.size() != 1)
		stop("value must be either length 1 or length(cd)");
	else
		cd.attr("latlon") = value;
	validate<NumericVector, Coord>(cd);
	return cd;
}


/// __________________________________________________
/// Validate coords vector
//' @rdname validate
// [[Rcpp::export(name = "validate.coords")]]
NumericVector validatecoords(NumericVector x, const bool force = true)
{
//	cout << "——Rcpp::export——validatecoords(NumericVector, const bool) format " << get_fmt_attribute(x) << endl;
	checkinherits(x, "coords");
	if (force)
		return validate<NumericVector, Coord>(x);
	else {
		if (!check_valid(x))
			warning("Invalid coords!");
		return x;
	}
}


/// __________________________________________________
/// Format coords vector - S3 method format.coords()
//' @rdname format
// [[Rcpp::export(name = "format.coords")]]
CharacterVector formatcoords(NumericVector x, bool usenames = true)
{
//	cout << "——Rcpp::export——formatcoords(NumericVector)\n";
	checkinherits(x, "coords");
	if(!x.size())
		stop("x has 0 length!");
	if (!check_valid(x))
		warning("Formatting invalid coords!");
	return wrap(Coord(get_coordtype(x), x).format(usenames));
}


/// __________________________________________________
/// Create waypoints
//' @rdname waypoints
// [[Rcpp::export(name = "as_waypoints.default")]]
DataFrame as_waypoints(DataFrame object, const int fmt = 1)
{
//	cout << "——Rcpp::export——as_waypoints(DataFrame, const int)\n";
	checkinherits(object, "data.frame");
	CoordType type = get_coordtype(fmt);
	object.attr("fmt") = fmt;
	int namescol = 0;
	if (!object.hasAttribute("namescol")) {
		namescol = nameinobj(object, "name");
		if (++namescol)
			object.attr("namescol") = namescol;
	}
	if (!object.hasAttribute("llcols")) {
		const vector<int> llcols { namescol + 1, namescol + 2 };
		object.attr("llcols") = llcols;
	}
	if(!valid_ll(object))
		stop("Invalid llcols attribute!");
	convert_switch<DataFrame, WayPoint>(object, type);
	object.attr("class") = CharacterVector{"waypoints", "data.frame"};
	return object;
}


/// __________________________________________________
/// Convert waypoints format
//' @rdname convert
// [[Rcpp::export(name = "convert.waypoints")]]
DataFrame convertwaypoints(DataFrame x, const int fmt)
{
//	cout << "——Rcpp::export——convertwaypoints(DataFrame, int) from " << get_fmt_attribute(x) << " to " << fmt << endl;
	checkinherits(x, "waypoints");
	CoordType newtype = get_coordtype(fmt);
	CoordType type = get_coordtype(x);
	if (newtype == type) {
//		cout << "——fmt out == fmt in!——" << endl;
		if (!check_valid(x))
			stop("Invalid waypoints!");
	} else {
		if(!valid_ll(x))
			stop("Invalid llcols attribute!");
		convert_switch<DataFrame, WayPoint>(x, newtype);
	}
	return x;
}


/// __________________________________________________
/// Validate waypoints vector
//' @rdname validate
// [[Rcpp::export(name = "validate.waypoints")]]
DataFrame validatewaypoints(DataFrame x, const bool force = true)
{
//      cout << "——Rcpp::export——validatewaypoints(DataFrame, const bool) format " << get_fmt_attribute(x) << endl;
        checkinherits(x, "waypoints");
        if(!valid_ll(x))
                stop("Invalid llcols attribute!");
	if (force)
		return validate<DataFrame, WayPoint>(x);
	else {
		if (!check_valid(x))
			warning("Invalid waypoints!");
		return x;
	}
}


/// __________________________________________________
/// Format waypoints vector - S3 method format.waypoints()
//' @rdname format
// [[Rcpp::export(name = "format.waypoints")]]
CharacterVector formatwaypoints(DataFrame x, bool usenames = true)
{
//	cout << "——Rcpp::export——formatwaypoints(DataFrame)\n";
	checkinherits(x, "waypoints");
	if(!x.nrows())
		stop("x has 0 rows!");
	if(!valid_ll(x))
		stop("Invalid llcols attribute!");
	if (!check_valid(x))
		warning("Formatting invalid waypoints!");
	return wrap(WayPoint(get_coordtype(x), x).format(usenames));
}


/// __________________________________________________
/// Latitude and longitude headers for S3 print.waypoint()
//' @rdname format
// [[Rcpp::export]]
CharacterVector ll_headers(const CharacterVector aswidth, const int fmt)
{
//	cout << "@ll_headers(int, const int)  width " <<  width << " fmt " << fmt << endl;
	constexpr int spacing[] { 5,  7,  8,
							 11, 13, 14 };
	vector<string> sv {
			string("Latitude") + string(spacing[fmt - 1], ' ') + "Longitude ",
			string(spacing[fmt + 2], '_') + string(2, ' ') + string(spacing[fmt + 2] + 1, '_')
		};

	constexpr int adjust[] = { 2, 6, 10 };
	const int width = (as<vector<string>>(aswidth)[0]).size() - adjust[fmt - 1];
	std::ostringstream ostrstr;
	transform(sv.begin(), sv.end(), sv.begin(), [&ostrstr, width](const string& s)
		{ ostrstr.str(""); ostrstr << std::setw(width) << s; return ostrstr.str(); });
	return wrap(sv);
}


/// __________________________________________________
/// Clone coords object from waypoints vector
//' @rdname coords
// [[Rcpp::export(name = "as_coords.waypoints")]]
NumericVector as_coordswaypoints(DataFrame object, bool which)
{
//	cout << "——Rcpp::export——as_coord(DataFrame)\n";
	checkinherits(object, "waypoints");
	NumericVector nv = object[get_vec_attr<DataFrame, int>(object, "llcols")[which ? 0 : 1] - 1];
	nv = clone(nv);
	nv.attr("class") = "coords";
	nv.attr("fmt") = object.attr("fmt");
	nv.attr("valid") = object.attr(which ? "validlat" : "validlon");
	nv.attr("latlon") = which;
	nv.attr("names") = getnames(object);
	return nv;
}


/// __________________________________________________
/// __________________________________________________