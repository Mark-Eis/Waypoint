/// __________________________________________________
/// CoordBase.cpp
/// __________________________________________________

// [[Rcpp::plugins(cpp23)]]

#include <Rcpp.h>
#include <cxxabi.h>
#include <array>

using namespace Rcpp;

using std::array;
using std::vector;
using std::string;
using std::string_view;
using namespace std::string_view_literals;
using std::transform;

#include "CoordBase.h"

#define FMT_HEADER_ONLY
#include "fmt/format.h"		// …fmt/*.h copied to ~/Documents/R/Packages/Waypoint/src/fmt. Works, but not always in pkgdown
#include "fmt/ranges.h"		// …fmt/*.h copied to ~/Documents/R/Packages/Waypoint/src/fmt. Works, but not always in pkgdown

#define DEBUG 1


/// __________________________________________________
/// __________________________________________________
/// Development and Debugging functions

#if DEBUG > 0

/// Report object construction and destruction
void _ctrsgn(const std::type_info& obj, bool destruct /* = false */)		// default arg should be in header
{
//	fmt::print("{}ing ", destruct ? "Destroy" : "Construct");
	std::fflush(nullptr);
	string s = obj.name();
	system(("c++filt -t " + s).data());
}

/// Format string for debugging code
constexpr auto exportstr { "——Rcpp::export——"sv };

#endif

/// Demangle object names
const string demangle(const std::type_info& obj)
{
	int status = 0;
	char* p { abi::__cxa_demangle(obj.name(), NULL, NULL, &status) };
	string str { fmt::format("\"{}\" (status {})", p, std::to_string(status)) };
	std::free(p);
	return str;
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
template<NumericVector_or_DataFrame T, class U> 
inline vector<U> get_vec_attr(const T& t, const char* attrname)
{
//	fmt::print("@{} attr=\"{}\" {}\n", "get_vec_attr<T, U>(const T&, const char*)", attrname, t.hasAttribute(attrname) ? true : false);
	return t.hasAttribute(attrname) ? as<vector<U>>(t.attr(attrname)) : vector<U>();
}


/// __________________________________________________
/// Return "fmt" attribute as int
template<NumericVector_or_DataFrame T>
inline int get_fmt_attribute(const T& t)
{
//	fmt::print("@{} fmt={}\n", "get_fmt_attribute<T>(const T&)", as<int>(t.attr("fmt")));
	return as<int>(t.attr("fmt"));
}


/// __________________________________________________
/// Check whether a NumericVector or DataFrame has a specified logical vector attribute and whether all true
template<NumericVector_or_DataFrame T>
int check_logical_attr(T t, const char* attrname)
{
//	fmt::print("@check_logical_attr<NumericVector_or_DataFrame>(T, const char*); T: {}; attrname {}\n", demangle(typeid(t)), attrname);
	const vector vec_attr{ get_vec_attr<T, bool>(t, attrname) };
	if (vec_attr.size()) {
		return all_of(vec_attr.begin(), vec_attr.end(), [](bool v) { return v;}) ? 0b11 : 0b01;
	} else {
		return 0b00;
	}
}


/// __________________________________________________
/// Does object inherit given class?
template<NumericVector_or_DataFrame T>
inline void checkinherits(T& t, const char* classname)
{
//	fmt::print("@{} T {} classname \"{}\"\n", "checkinherits<T>(T&, const char*)", demangle(typeid(t)), classname);
	if (!t.inherits(classname)) stop("Argument must be a \"%s\" object", classname);
}


/// __________________________________________________
/// Is item number present in object? (Using C++ numbering)
template<class T>
inline bool is_item_in_obj(const T t, int item)
{
//	fmt::print("@{} T {} item={}\n", "is_item_in_obj<T>(T, int)", demangle(typeid(t)), item);
	if (NA_INTEGER == item)
		return false;
	else
		return !(item < 0) && item < t.size();
}


/// __________________________________________________
/// Standarise width of strings in vector to that of the longest
inline void stdlenstr(vector<string>& sv)
{
//	fmt::print("@{}\n", "stdlenstr(vector<string>&)");
	int maxwdth = max_element(sv.begin(), sv.end(), [](const string& a, const string& b){ return a.size() < b.size(); })->size();
	transform(sv.begin(), sv.end(), sv.begin(), [maxwdth](const string& s) { return fmt::format("{:<{}}", s, maxwdth); });
}


/// __________________________________________________
/// Prefix vector<string> elements with elements of vector<T>—default for vector<string> prefix
template<class T>
inline void prefixvecstr(vector<string>& sv, const vector<T>& prefix)
{
//	fmt::print("@{} T {}\n", "prefixvecstr<T>(vector<string>&, const vector<T>&) [default]", "vector<string>");
	transform(sv.begin(), sv.end(), prefix.begin(), sv.begin(), [](string& lls, const string& name) { return name + "  " + lls; });	
}


/// __________________________________________________
/// Specialisation for vector<int> prefix
template<>
inline void prefixvecstr(vector<string>& sv, const vector<int>& prefix)
{
//	fmt::print("@{}\n", "prefixvecstr<>(vector<string>&, const vector<int>&)");
	transform(sv.begin(), sv.end(), prefix.begin(), sv.begin(), [](string& lls, int name) { return std::to_string(name) + "  " + lls; });	
}


/// __________________________________________________
/// Prefix vector<string> elements with elements of RObject 
inline bool prefixwithnames(vector<string>& sv, RObject& namesobj)
{
//	fmt::print("@{}\n", "prefixwithnames(vector<string>&, RObject&)");
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
//	fmt::print("@{}\n", "str_tolower(string)");
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return tolower(c); });
    return s;
}


/// __________________________________________________
/// Find position of name within object names
template<List_or_DataFrame T>
int nameinobj(const T t, const char* name)
{
//	fmt::print("@{} name={}\n", "nameinobj<T>(const T, const char*)", name);
	vector names{ get_vec_attr<T, string>(t, "names") };
	if (!names.size())
		return -1;
	typedef decltype(names.size()) Tmp;
	Tmp i = 0;
	for (auto str : names ) {
//		fmt::print("@{} testing {}\n", "nameinobj<T>(const T, const char*)", str);
		if (!str_tolower(str).compare(name)) {
//			fmt::print("@{} found {}\n", "nameinobj<T>(const T, const char*)", str);
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
//	fmt::print("@{}\n", "getnames(const DataFrame)");
	vector namescolvec{ get_vec_attr<DataFrame, int>(df, "namescol") };
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
#if DEBUG > 0

auto fmt::formatter<CoordType>::format(CoordType ct, format_context& ctx) const
	-> format_context::iterator
{
	constexpr array<const char*, 3> names {"DecDeg", "DegMin", "DegMinSec"};
	return formatter<string_view>::format(names[fmt::underlying(ct)], ctx);
}

#endif

/// __________________________________________________
/// Convert int to CoordType enum
inline const CoordType get_coordtype(int i)
{
//	fmt::print("@{} {}\n", "get_coordtype(int)" , i);
	if (i < 1 || i > 3)
		stop("\"fmt\" must be between 1 and 3");
	using enum CoordType;
	constexpr array<CoordType, 3> coordtypes{ decdeg, degmin, degminsec };
	return coordtypes[i - 1];
}


/// __________________________________________________
/// Convert "fmt" attribute to CoordType enum
template<NumericVector_or_DataFrame T>
inline const CoordType get_coordtype(const T& t)
{
//	fmt::print("@{} t {}\n", "get_coordtype<T>(const T&)", demangle(typeid(t)));
	return get_coordtype(get_fmt_attribute(t));
}


/// __________________________________________________
/// Convert CoordType enum to int
inline int coordtype_to_int(CoordType ct)
{
//	fmt::print("@{} ct={}\n", "coordtype_to_int(CoordType)", ct);
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


/// ——— ¡¡Nothing new above this line!! ———


/// __________________________________________________
/// __________________________________________________
/// Coordlet class

/// __________________________________________________
/// Constructor of Coordlet
template<CoordType current_type>
Coordlet<current_type>::Coordlet(NumericVector nv) :
	nv{ nv },
	latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") }
{
//	fmt::print("§{} {} ", fmt::format("Coordlet<CoordType::{}>::Coordlet(NumericVector)", current_type), "");  _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Format Coordlet::nv as vector<string> of CoordType
template<CoordType current_type> template<CoordType required_type>
vector<string> Coordlet<current_type>::format0(bool wpt) const
{
//	fmt::print("@Coordlet<CoordType::{}>::format0<CoordType::{}>() const; wpt: {}\n", current_type, required_type, wpt);

	vector<bool>::const_iterator ll_it { latlon.begin() };
	const auto ll_size { latlon.size() };
	vector out_sv{ vector<string>(nv.size()) };

	if constexpr (CoordType::decdeg == required_type)
		transform(nv.begin(), nv.end(), out_sv.begin(), [this](double n){
				return fmt::format("{:>{}.{}f}\u00B0", ff.get_decdeg(n), 11, 6);
			});
	else if constexpr (CoordType::degmin == required_type)
		transform(nv.begin(), nv.end(), out_sv.begin(), [this](double n){
				return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) +
					   fmt::format("{:0>{}.{}f}\u2032", fabs(ff.get_decmin(n)), 7, 4);
			});
	else
		transform(nv.begin(), nv.end(), out_sv.begin(), [this](double n){
				return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) +
					   fmt::format("{:0>{}}\u2032", abs(ff.get_min(n)), 2) +
					   fmt::format("{:0>{}.{}f}\u2033", fabs(ff.get_sec(n)), 5, 2);
			});

	if (wpt) {
	} else {
		if constexpr (CoordType::decdeg == required_type) {
			if (ll_size) 
				transform(out_sv.begin(), out_sv.end(), nv.begin(), out_sv.begin(), [&ll_it, &ll_size](string& out_sv, double n){
					return out_sv + ((ll_size > 1 ? *ll_it++ : *ll_it) ? " lat" : " lon");});
		} else
			transform(out_sv.begin(), out_sv.end(), nv.begin(), out_sv.begin(), [&ll_it, &ll_size](string& out_sv, double n){
				return out_sv + (ll_size ? cardpoint(n < 0, ll_size > 1 ? *ll_it++ : *ll_it) : cardi_b(n < 0));});
	}

	return out_sv;
}


/// __________________________________________________
/// __________________________________________________
/// CoordType switches

/// __________________________________________________
/// Switch CoordType required for Coordlet<CoordType>::format0()
template<CoordType current_type>
vector<string> Coordlet<current_type>::format_switch(CoordType required_type, bool wpt) const
{
//	fmt::print("@Coordlet<CoordType::{}>::format_switch(CoordType); required: {}; wpt: {}\n", current_type, required_type, wpt);

	switch (required_type)
	{
		case CoordType::decdeg:
			return format0<CoordType::decdeg>(wpt);

		case CoordType::degmin:
			return format0<CoordType::degmin>(wpt);

		case CoordType::degminsec:
			return format0<CoordType::degminsec>(wpt);

		default:
			stop("Coordlet<CoordType>::format_switch(CoordType) my bad");
	}
}


/// __________________________________________________
/// Convert Coordlet::nv to a new CoordType
template<CoordType current_type> template<CoordType required_type>
void Coordlet<current_type>::convert0()
{
//	fmt::print("@Coordlet<CoordType::{}>::convert0<CoordType::{}>()\n", current_type, required_type);

	if constexpr (CoordType::decdeg == required_type)
		transform(nv.begin(), nv.end(), nv.begin(), [this](double n){
				return ff.get_decdeg(n);
			});
	else if constexpr (CoordType::degmin == required_type)
		transform(nv.begin(), nv.end(), nv.begin(), [this](double n){
				return ff.get_deg(n) * 1e2 + ff.get_decmin(n);
			});
	else
		transform(nv.begin(), nv.end(), nv.begin(), [this](double n){
				return ff.get_deg(n) * 1e4 + ff.get_min(n) * 1e2 + ff.get_sec(n);
			});
}


/// __________________________________________________
/// Switch CoordType required for Coordlet<CoordType>::convert0()
template<CoordType current_type>
void Coordlet<current_type>::convert_switch(CoordType required_type)
{
//	fmt::print("@Coordlet<CoordType::{}>::convert_switch(CoordType); required_type: {}\n", current_type, required_type);
	validate(true, "coords to be converted");
	switch (required_type)
	{
		case CoordType::decdeg:
			convert0<CoordType::decdeg>();
			break;

		case CoordType::degmin:
			convert0<CoordType::degmin>();
			break;

		case CoordType::degminsec:
			convert0<CoordType::degminsec>();
			break;

		default:
			stop("Coordlet<CoordType>::convert_switch(CoordType) my bad");
	}
}


/// __________________________________________________
/// Validate Coordlet::nv
template<CoordType current_type>
void Coordlet<current_type>::validate(bool warn, const char* what)
{
//	fmt::print("@Coordlet<CoordType::{}>::validate(bool); warn: {}; what: {}; latlon: {}\n", current_type, warn, what, fmt::join(latlon, ", "));
	vector<bool>::const_iterator ll_it{ latlon.begin() };
	auto ll_size { latlon.size() };

	valid.assign(nv.size(), {false});
	transform(nv.begin(), nv.end(), valid.begin(), [this, &ll_it, &ll_size](double n){
		return !((fabs(ff.get_decdeg(n)) > (ll_size && (ll_size > 1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(fabs(ff.get_decmin(n)) >= 60) ||
				(fabs(ff.get_sec(n)) >= 60));
	});

	if (all_of(valid.begin(), valid.end(), [](bool v) { return v;}))
		valid.assign({true});
	else
		if (warn)
			warning("Validation%s%s failed!", strlen(what) ? " of " : "", what);
	nv.attr("valid") = valid;
}


/// __________________________________________________
/// Switch current CoordType to format nv
vector<string> format_switch_current(NumericVector nv, CoordType current_type, CoordType required_type, bool wpt)
{
//	fmt::print("@format_switch_current(NumericVector, CoordType, CoordType, bool); current: {}; required: {}; wpt: {}\n", current_type, required_type, wpt);
	using enum CoordType;

	switch (current_type)
	{
		case decdeg:
			return format_required<decdeg>(nv, required_type, wpt);

		case degmin:
			return format_required<degmin>(nv, required_type, wpt);

		case degminsec:
			return format_required<degminsec>(nv, required_type, wpt);

		default:
			stop("format_switch_current(NumericVector, CoordType) my bad");
	}
}


/// __________________________________________________
/// Switch required CoordType to format nv
template<CoordType current_type> 
inline vector<string> format_required(NumericVector nv, CoordType required_type, bool wpt)
{
//	fmt::print("@format_required<CoordType::{}>(NumericVector, CoordType); required: {}; wpt: {}\n", current_type, required_type, wpt);
	return Coordlet<current_type>{ nv }.format_switch(required_type, wpt);
}


/// __________________________________________________
/// Switch current CoordType to convert nv
void convert_switch_current(NumericVector nv, CoordType current_type, CoordType required_type)
{
//	fmt::print("@convert_switch_current(NumericVector, CoordType, CoordType); current_type: {}; required_type: {}\n", current_type, required_type);
	using enum CoordType;

	if (required_type != current_type) {
		switch (current_type)
		{
			case decdeg:
				convert_required<decdeg>(nv, required_type);
				break;
	
			case degmin:
				convert_required<degmin>(nv, required_type);
				break;
	
			case degminsec:
				convert_required<degminsec>(nv, required_type);
				break;
	
			default:
				stop("convert_switch_current(NumericVector, CoordType) my bad");
		}
		nv.attr("fmt") = coordtype_to_int(required_type) + 1;
	} else {
		fmt::print("——fmt out == fmt in!——\n"); std::fflush(nullptr);
		if (!check_valid(nv))
			stop("Invalid coords!");
	}
}


/// __________________________________________________
/// Switch required CoordType to convert nv
template<CoordType current_type> 
inline void convert_required(NumericVector nv, CoordType required_type)
{
//	fmt::print("@convert_required<CoordType::{}>(NumericVector, CoordType); required_type: {}\n", current_type, required_type);
	Coordlet<current_type>{ nv }.convert_switch(required_type);
}


/// __________________________________________________
/// Switch current CoordType to validate nv
void validate_switch_current(NumericVector nv, CoordType current_type, bool warn, const char* what)
{
//	fmt::print("@validate_switch_current(NumericVector, CoordType, bool, const char*); current_type: {}; warn: {}; what: {}\n", current_type, warn, what);
	using enum CoordType;

	switch (current_type)
	{
		case decdeg:
			Coordlet<decdeg>{ nv }.validate(warn, what);
			break;

		case degmin:
			Coordlet<degmin>{ nv }.validate(warn, what);
			break;

		case degminsec:
			Coordlet<degminsec>{ nv }.validate(warn, what);
			break;

		default:
			stop("validate_switch_current(NumericVector, CoordType, bool, const char*) my bad");
	}
}


/// __________________________________________________
/// __________________________________________________
/// Validation functions

/// __________________________________________________
/// Check "valid" attribute of NumericVector all true
bool check_valid(const NumericVector nv)
{
//	fmt::print("@check_valid(const NumericVector)\n");
	int validated = check_logical_attr(nv, "valid");
//	fmt::print("@check_valid(const NumericVector); validated is: {}; all valid is: {}\n", bool(validated), bool(validated >> 1));  // To be deleted
	if (!validated)
		return revalidate(nv);
	return validated >> 1;
}


/// __________________________________________________
/// Check "lat_valid" and "lon_valid attributes of DataFrame are all true
bool check_valid(const DataFrame df)
{
//	fmt::print("@check_valid(const DataFrame)\n");

	int latvalidated = check_logical_attr(df, "validlat");
	if (!latvalidated)
		return revalidate(df);

	int lonvalidated = check_logical_attr(df, "validlon");
	if (!lonvalidated)
		return revalidate(df);

	if (!(latvalidated >> 1))
		warning("Invalid latitude!");
	if (!(lonvalidated >> 1))
		warning("Invalid longitude!");
	return latvalidated >> 1 || lonvalidated >> 1;
}


/// __________________________________________________
/// Revalidate NumericVector —— Temporary solution, needs to work for DataFrame
template<NumericVector_or_DataFrame T>
bool revalidate(const T t)
{
//	fmt::print("@revalidate<T>(const T); T: {}\n", demangle(typeid(t)));
	warning("Revalidating %s…!", demangle(typeid(t)));
	if constexpr (std::is_same_v<T, NumericVector>)
		validate(t);										// Temporary solution
	else
		stop("Compiler pacifier!");							// Temporary solution
	return check_valid(t);
}


/// __________________________________________________
/// Validate NumericVector or DataFrame —— Temporary solution, needs to work for DataFrame
template<NumericVector_or_DataFrame T>
inline const T validate(const T t)
{
//	fmt::print("@validate<T>(const T); T: {}\n", demangle(typeid(t)));
	validate_switch_current(t, get_coordtype(t), true, "coords [validate<T>(const T)]");			//¡¡¡Temporary solution!!!
	return t;	
}


/// __________________________________________________
/// Check df has valid "llcols" attribute
bool valid_ll(const DataFrame df)
{
//	fmt::print("@{}\n", "valid_ll(const DataFrame)");
	bool valid = false;
	vector llcols { get_vec_attr<DataFrame, int>(df, "llcols") };
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
NumericVector as_coords(NumericVector object, int fmt = 1)
{
//	fmt::print("{1}@{0} fmt={2}\n", "as_coords(NumericVector, int)", exportstr, fmt);
	object.attr("fmt") = fmt;
	validate_switch_current(object, get_coordtype(fmt), true, "coords");
	object.attr("class") = "coords";
	return object;
}


/// __________________________________________________
/// Convert coords format
//' @rdname convert
// [[Rcpp::export(name = "convert.coords")]]
NumericVector convertcoords(NumericVector x, int fmt)
{
	checkinherits(x, "coords");
	CoordType type = get_coordtype(x);
	CoordType newtype = get_coordtype(fmt);
//	fmt::print("{1}@{0} from {2} to {3}\n", "convertcoords(NumericVector, int)", exportstr, type, newtype);
	convert_switch_current(x, type, newtype);
	return x;
}


/// __________________________________________________
/// Set latlon attribute on "coords" NumericVector and revalidate
//' @rdname coords
// [[Rcpp::export(name = "`latlon<-`")]]
NumericVector latlon(NumericVector cd, LogicalVector value)
{
//	fmt::print("{1}@{0}\n", "latlon(NumericVector, LogicalVector)", exportstr);
	checkinherits(cd, "coords");
	if (value.size() != cd.size() && value.size() != 1)
		stop("value must be either length 1 or length(cd)");
	else
		cd.attr("latlon") = value;
	validate_switch_current(cd, get_coordtype(cd), true, "coords [`latlon<-`]");
	return cd;
}


/// __________________________________________________
/// Validate coords vector
//' @rdname validate
// [[Rcpp::export(name = "validate.coords")]]
NumericVector validatecoords(NumericVector x, bool force = true)
{
//	fmt::print("{1}@{0} force: {2}\n", "validatecoords(NumericVector, bool)", exportstr, force);
	checkinherits(x, "coords");
	if (force)									
		validate_switch_current(x, get_coordtype(x), true, "coords [validatecoords]");
	else {
		if (!check_valid(x))
			warning("Invalid coords!");
	}
	return x;
}


/// __________________________________________________
/// Format coords vector - S3 method format.coords()
//' @rdname format
// [[Rcpp::export(name = "format.coords")]]
CharacterVector formatcoords(NumericVector x, bool usenames = true, bool validate = true, int fmt = 0)
{
//	fmt::print("{1}@{0} usenames: {2}, validate: {3}\n", "formatcoords(NumericVector, bool, bool)", exportstr, usenames, validate);
	checkinherits(x, "coords");
	if(!x.size())
		stop("x has 0 length!");
	if (validate)
		if (!check_valid(x))
			warning("Formatting invalid coords!");

	CoordType ct { get_coordtype(x) };
	vector sv{ format_switch_current(x, ct, fmt ? get_coordtype(fmt) : ct) };
	vector names{ get_vec_attr<NumericVector, string>(x, "names") };
	if (names.size() && usenames) {
		stdlenstr(names);
		prefixvecstr(sv, names);
	}
	return wrap(sv);
}


/// __________________________________________________
/// Create waypoints
//' @rdname waypoints
// [[Rcpp::export(name = "as_waypoints.default")]]
DataFrame as_waypoints(DataFrame object, int fmt = 1)
{
	fmt::print("{1}@{0} fmt={2}\n", "as_waypoints(DataFrame, int)", exportstr, fmt);
	object.attr("fmt") = fmt;
	int namescol = 0;
	if (!object.hasAttribute("namescol")) {
		namescol = nameinobj(object, "name");
		if (++namescol)
			object.attr("namescol") = namescol;
	}
	if (!object.hasAttribute("llcols")) {
		const vector llcols{ namescol + 1, namescol + 2 };
		object.attr("llcols") = llcols;
	}
	if(!valid_ll(object))
		stop("Invalid llcols attribute!");
//	WayPoint{get_coordtype(fmt), object}.validate();
	object.attr("class") = CharacterVector{"waypoints", "data.frame"};
	return object;
}


/// __________________________________________________
/// __________________________________________________

/*** R

as_coords <- function(object, ...) 
    UseMethod("as_coords")

convert <- function(x, ...) 
    UseMethod("convert")

validate <- function(x, ...) 
    UseMethod("validate")

as_waypoints <- function(object, ...) 
    UseMethod("as_waypoints")



## ========================================
##  S3method print.coords(x, ...)
#'
#' @rdname format
#' @export

print.coords <- function (x, ..., max = NULL) {
    n <- length(x)
    validate(x, force = FALSE)
    max <- max %||% getOption("max.print", 99999L)
    if (!is.finite(max)) 
        stop("invalid 'max' / getOption(\"max.print\"): ", max)
    omit <- (n0 <- max %/% (if (is.null(names(x))) 1L else 2L)) < n
    if (omit)
        x0 <- x[seq_len(n0)]
    else
        x0 <- x
    writeLines(format(x0, validate = FALSE, ...))
    if (omit) 
        cat(" [ reached 'max' / getOption(\"max.print\") -- omitted", n - n0, "entries ]\n")
    invisible(x)
}


## ========================================
##  S3method `[.coords`(x, i)
#'
#' @export


`[.coords` <- function(x, i) 
{
    if (missing(i)) return(x)
    if (max(i) > length(x)) stop(gettext("subscript out of bounds"), domain = NA, call. = FALSE)
    attrx <- attributes(x)
    y <- NextMethod("[")
    attributes(y) <- lapply(attrx[names(attrx) %in% c("fmt", "latlon", "names", "valid")], \(ax) if(length(ax) == length(x)) ax[i] else ax)
    class(y) <- oldClass(x)
    y
}


## ========================================
##  S3method `[.coords<-`(x, i, value)
#'
#' @rdname Extract
#' @export


`[<-.coords` <- function(x, i, value) 
{
    if (max(i) > length(x)) stop(gettext("subscript out of bounds"), domain = NA, call. = FALSE)
    vfmt <- attr(value, "fmt", exact = TRUE)
    if (!is.null(vfmt) && attr(x, "fmt", exact = TRUE) != vfmt) stop("value has different fmt attribute\n", call. = FALSE)
    y <- NextMethod("[<-")
    vll <- attr(value, "latlon", exact = TRUE)
    if (!is.null(vll)) {
        if (length(vll) != length(value) && length(vll) != 1)
            stop("\"latlon\" attribute of replacement value must be either length 1 or length(value)\n", call. = FALSE)
        xll <- attr(x, "latlon", exact = TRUE)
        if (length(xll) == length(x))
            attr(y, "latlon")[i] <- vll
        else if (length(xll) == 1 && !all(vll == xll)) {
            attr(y, "latlon") <- rep(xll, length(x))
            attr(y, "latlon")[i] <- vll
        }
    }
    class(y) <- oldClass(x)
    validate(y)
}

## ========================================

# newds <- T
newds <- F

if (newds) {

	dd <- c(51.5077650, 49.5462100, 48.1072317, 38.8894933, 0.0000000, -37.1117400, -53.1047817, -25.2401550,
			-0.1279233, 18.3985617, -122.7786717, -77.0352417, 0.0000000, -12.2886300, 73.5172833, -57.519227)
	
	dm <- c(5130.4659, 4932.7726, 4806.4339, 3853.3696, 0.0000, -3706.7044, -5306.2869, -2514.4093,
			-007.6754, 1823.9137, -12246.7203, -7702.1145, 0.0000, -1217.3178, 7331.0370, -5731.1536)
	
	dms <- c(513027.95, 493246.36, 480626.04, 385322.18, 00000.00, -370642.26, -530617.21, -251424.56,
			-00740.53, 182354.82, -1224643.22, -770206.87, 00000.00, -121719.07, 733102.22, -573109.21)
			
	names(dd) <- names(dm) <- names(dms) <- rep(c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
                   "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"), 2)

	wp1 <- data.frame(
	    name = c("Nelson's Column", "Ostravice", "Tally Ho", "Washington Monument", "Null Island",
	             "Tristan da Cunha", "Mawson Peak", "Silvio Pettirossi International Airport"),
	    lat = c(513027.95, 493246.36, 480626.04, 385322.18, 0, -370642.26, -530617.21, -251424.56),
	    lon = c(-00740.53, 182354.82, -1224643.22, -770206.87, 0, -121719.07, 733102.22, -573109.21)
	)

}
*/



/// __________________________________________________
/// "Testlets"
/// __________________________________________________
/// __________________________________________________
// [[Rcpp::export]]
void bin()
{
	Rcout << "Hello Mimiland!\n\n";

	Rcout << "sizeof char: " << sizeof(char) << "\n\n";
	Rcout << "sizeof char: " << sizeof(char) << "\n\n";
	Rcout << "sizeof char: " << sizeof(char) << "\n\n";
	Rcout << "sizeof short int: " << sizeof(short int) << "\n\n";
	Rcout << "sizeof int: " << sizeof(int) << "\n\n";
	Rcout << "sizeof long int: " << sizeof(long int) << "\n\n";
	Rcout << "sizeof vector<bool>: " << sizeof(vector<bool>) << "\n\n";
	
	Rcout << "sizeof 0b1: " << sizeof(0b1) << "\n\n";
	Rcout << "sizeof 0b10: " << sizeof(0b10) << "\n\n";
	auto Ob1_auto = 0b1;
	auto Ob2_auto = 0b10;
	Rcout << "sizeof Ob1_auto: " << sizeof(Ob1_auto) << "; "<< demangle(typeid(Ob1_auto)) << "\n\n";
	Rcout << "sizeof 0b10_auto: " << sizeof(Ob2_auto) << "; "<< demangle(typeid(Ob2_auto)) << "\n\n";
	char Ob1_char = 0b1;
	char Ob2_char = 0b10;
	Rcout << "sizeof Ob1_char: " << sizeof(Ob1_char) << "; "<< demangle(typeid(Ob1_char)) << "\n\n";
	Rcout << "sizeof 0b10_char: " << sizeof(Ob2_char) << "; "<< demangle(typeid(Ob2_char)) << "\n\n";
	Rcout << "sizeof 0b10_char: " << sizeof(Ob2_char) << "; "<< demangle(typeid(Ob2_char)) << "\n\n";

	Rcout << "0b1 is " << 0b1 << "; 0b10 is " << 0b10 << "; 0b1 & 0b10 is " << int(0b1 & 0b10) << "; 0b1 | 0b10 is " << int(0b1 | 0b10) << "\n\n";
}


/// __________________________________________________
/// __________________________________________________
// [[Rcpp::export]]
CharacterVector foo(NumericVector object, int fmt_cur = 1, int fmt_req = 1)
{
	Rcout << "Hello Mimiland!\n";
	Rcout << "object: " << object << "\n\n";

	CoordType ct_cur { get_coordtype(fmt_cur) };
	CoordType ct_req { get_coordtype(fmt_req) };
	return wrap(format_switch_current(object, ct_cur, ct_req));
}


/// __________________________________________________
/// __________________________________________________
// [[Rcpp::export]]
NumericVector bar(NumericVector x, int fmt_cur = 1, int fmt_req = 1)
{
	Rcout << "Hello Mimiland!\n";
	Rcout << "x: " << x << "\n\n";

	CoordType ct_cur { get_coordtype(fmt_cur) };
	CoordType ct_req { get_coordtype(fmt_req) };
	convert_switch_current(x, ct_cur, ct_req);
	return x;
}


/// __________________________________________________
/// __________________________________________________
// [[Rcpp::export]]
NumericVector bas(NumericVector x, int fmt_cur = 1, const char* what = "bas")
{
	Rcout << "Hello Mimiland!\n";
	Rcout << "x: " << x << "\n\n";

	CoordType ct_cur { get_coordtype(fmt_cur) };
	validate_switch_current(x, ct_cur, true, what);
	return x;
}


/// __________________________________________________
/// __________________________________________________
// [[Rcpp::export]]
int cla_test(RObject object, const char* attrname = "valid")
{
	fmt::print("Testing: -\n\n{}\n{}\n\n",
		"template<NumericVector_or_DataFrame T>",
		"int check_logical_attr(T t, const char* attrname)"
	);

	fmt::print("attrname is {}\n\n", attrname);

    if (is<NumericVector>(object)) {
    	fmt::print("object is a {}\n\n", "NumericVector");
        NumericVector nv = as<NumericVector>(object);
		return check_logical_attr(nv, attrname);
    } else if (is<DataFrame>(object)) {
    	fmt::print("object is a {}\n\n", "DataFrame");
        DataFrame df = as<DataFrame>(object);
		return check_logical_attr(df, attrname);
    } else
        stop("cla_test(): object not a NumericVector or DataFrame ");
}
