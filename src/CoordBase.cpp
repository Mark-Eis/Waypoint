
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

#define DEBUG 0


/// __________________________________________________
/// __________________________________________________
/// Development and Debugging functions

#if DEBUG > 0

/// Report object construction and destruction
void _ctrsgn(const std::type_info& obj, bool destruct)
{
	fmt::print("{}ing ", destruct ? "Destroy" : "Construct");
	std::fflush(nullptr);
	string s = obj.name();
	system(("c++filt -t " + s).data());
}

/// Format string for debugging code
//constexpr string_view exportstr { "——Rcpp::export——" };
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


/// __________________________________________________
/// __________________________________________________
/// Instantiate FamousFive Derived Classes

vector<FamousFive*> vff { &ff_decdeg, &ff_degmin, &ff_degminsec };


/// __________________________________________________
/// __________________________________________________
/// CoordType switches

/// __________________________________________________
/// Convert coords or waypoints format CoordType switch 
template<NumericVector_or_DataFrame T, Coord_or_WayPoint U>
void convert_switch(T t, CoordType newtype)
{
	CoordType type = get_coordtype(t);
//	fmt::print("@{} T: {} oldtype: {}, newtype: {}\n", "convert_switch<T&, U>(T, CoordType)", demangle(typeid(t)), type, newtype);
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
template<Coord_or_WayPoint T>
vector<string> format_switch(const T& t, CoordType ctreq)
{
//	fmt::print("@{} T: {} CoordType::{}, ctreq CoordType::{}\n", "format_switch<T>(const T&, CoordType)", demangle(typeid(t)), t.get_coordtype(), ctreq);
	switch (ctreq)
	{
		case CoordType::decdeg:
			return t.template format<CoordType::decdeg>();

		case CoordType::degmin:
			return t.template format<CoordType::degmin>();

		case CoordType::degminsec:
			return t.template format<CoordType::degminsec>();

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
//	fmt::print("§{} {} ", "Coordbase::Coordbase(CoordType)", ct); _ctrsgn(typeid(*this));
}


Coordbase::~Coordbase()
{
//	fmt::print("§{} {} ", "Coordbase::~Coordbase()", ct); _ctrsgn(typeid(*this), true);
}


CoordType Coordbase::get_coordtype() const
{
//	fmt::print("@{} ct={}\n", "Coordbase::get_coordtype()", coordtype_to_int(ct));
	return ct;
}


/// __________________________________________________
/// Convert NumericVector CoordType
template<CoordType type>
void Coordbase::convert0(NumericVector nv)
{
//	fmt::print("@Coordbase::convert0<CoordType::{}>()\n", type);

	if constexpr (CoordType::decdeg == type)
		transform(nv.begin(), nv.end(), nv.begin(), [this](double n){
				return ff.get_decdeg(n);
			});
	else if constexpr (CoordType::degmin == type)
		transform(nv.begin(), nv.end(), nv.begin(), [this](double n){
				return ff.get_deg(n) * 1e2 + ff.get_decmin(n);
			});
	else
		transform(nv.begin(), nv.end(), nv.begin(), [this](double n){
				return ff.get_deg(n) * 1e4 + ff.get_min(n) * 1e2 + ff.get_sec(n);
			});
}


/// __________________________________________________
/// Validate coordinates in NumericVector
void Coordbase::validate0(NumericVector nv, vector<bool>& valid, const vector<bool>& latlon)
{
//	fmt::print("@{} latlon: {}\n", "Coordbase::validate0()", fmt::join(latlon, ", "));
	vector<bool>::const_iterator ll_it{ latlon.begin() };
	auto ll_size { latlon.size() };

	valid.assign(nv.size(), {false});
	transform(nv.begin(), nv.end(), valid.begin(), [this, &ll_it, &ll_size](double n){
		return !((fabs(ff.get_decdeg(n)) > (ll_size && (ll_size > 1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(fabs(ff.get_decmin(n)) >= 60) ||
				(fabs(ff.get_sec(n)) >= 60));
	});
}


/// __________________________________________________
/// Format coordinates as vector<string> of CoordType
template<CoordType type>
vector<string> Coordbase::format0(NumericVector nv) const
{
//	fmt::print("@Coordbase::format0<CoordType::{}>() const\n", type);
	vector outstr{ vector<string>(nv.size()) };

	if constexpr (CoordType::decdeg == type)
		transform(nv.begin(), nv.end(), outstr.begin(), [this](double n){
				return fmt::format("{:>{}.{}f}\u00B0", ff.get_decdeg(n), 11, 6);
			});
	else if constexpr (CoordType::degmin == type)
		transform(nv.begin(), nv.end(), outstr.begin(), [this](double n){
				return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) +
					   fmt::format("{:0>{}.{}f}\u2032", fabs(ff.get_decmin(n)), 7, 4);
			});
	else
		transform(nv.begin(), nv.end(), outstr.begin(), [this](double n){
				return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) +
					   fmt::format("{:0>{}}\u2032", abs(ff.get_min(n)), 2) +
					   fmt::format("{:0>{}.{}f}\u2033", fabs(ff.get_sec(n)), 5, 2);
			});
	
	return outstr;
}


/// __________________________________________________
/// Coordinate derived class

Coord::Coord(CoordType ct, NumericVector nv) :
	Coordbase(ct), nv(nv),
	latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") }
{
//	fmt::print("§{} {} ", "Coord::Coord(CoordType, NumericVector)", ct); _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Convert Coord NumericVector CoordType
template<CoordType newtype>
inline void Coord::convert()
{
//	fmt::print("@Coord::convert<{}>() to {}\n", ct, newtype);
	convert0<newtype>(nv);
}


/// __________________________________________________
/// Validate coords vector
void Coord::validate(bool warn)
{
//	fmt::print("@{} latlon: {}\n", "Coord::validate()", fmt::join(latlon, ", "));
	validate0(nv, valid, latlon);

	if (all_of(valid.begin(), valid.end(), [](bool v) { return v;}))
		valid.assign({true});
	else
		if (warn)
			warning("Validation failed!");
	nv.attr("valid") = valid;
}


/// __________________________________________________
/// Format coordinates as vector<string> of CoordType
template<CoordType type>
vector<string> Coord::format() const
{
//	fmt::print("@Coord::format<CoordType::{}>() const; ll_size type: {}, ll_size: {}\n", type, demangle(typeid(latlon.size())), latlon.size());
	vector<bool>::const_iterator ll_it { latlon.begin() };
	const auto ll_size { latlon.size() };
	vector out_sv{ format0<type>(nv) };

	if constexpr (CoordType::decdeg == type) {
		if (ll_size) 
			transform(out_sv.begin(), out_sv.end(), nv.begin(), out_sv.begin(), [&ll_it, &ll_size](string& outstr, double n){
				return outstr + ((ll_size > 1 ? *ll_it++ : *ll_it) ? " lat" : " lon");});
	} else
		transform(out_sv.begin(), out_sv.end(), nv.begin(), out_sv.begin(), [&ll_it, &ll_size](string& outstr, double n){
			return outstr + (ll_size ? cardpoint(n < 0, ll_size > 1 ? *ll_it++ : *ll_it) : cardi_b(n < 0));});

	return out_sv;
}


/// __________________________________________________
/// __________________________________________________
/// Waypoint class

WayPoint::WayPoint(CoordType ct, DataFrame df) :
	Coordbase(ct), df(df),
	nvlat(df[get_vec_attr<DataFrame, int>(df, "llcols")[0] - 1]), 
	nvlon(df[get_vec_attr<DataFrame, int>(df, "llcols")[1] - 1])
{
//	fmt::print("§{} {} ", "WayPoint::WayPoint(CoordType, DataFrame)", ct); _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Convert both WayPoint NumericVectors CoordType
template<CoordType newtype>
inline void WayPoint::convert()
{
//	fmt::print("@WayPoint::convert<{}>() to {}\n", ct, newtype);
	convert0<newtype>(nvlat);
	convert0<newtype>(nvlon);
}


/// __________________________________________________
/// Validate WayPoint
void WayPoint::validate(bool warn)
{
//	fmt::print("@{}\n", "WayPoint::validate(bool)");

	validate0(nvlat, validlat, vector{ true });
	validate0(nvlon, validlon, vector{ false });

	if (all_of(validlat.begin(), validlat.end(), [](bool v) { return v;}))
		validlat.assign({true});
	else
		if (warn)
			warning("Validation of latitude failed!");
	df.attr("validlat") = validlat;

	if (all_of(validlon.begin(), validlon.end(), [](bool v) { return v;}))
		validlon.assign({true});
	else
		if (warn)
			warning("Validation of longitude failed!");
	df.attr("validlon") = validlon;
}


/// __________________________________________________
/// Format waypoints as vector<string> of CoordType
template<CoordType type>
vector<string> WayPoint::format() const
{
//	fmt::print("@WayPoint::format<CoordType::{}>()\n", type);
	vector sv_lat{ format2<type>(true) };
	vector sv_lon{ format2<type>(false) };

	transform(sv_lat.begin(), sv_lat.end(), sv_lon.begin(), sv_lat.begin(), [](auto& latstr, auto& lonstr){return latstr + "  " + lonstr;});
	return sv_lat;
}


/// __________________________________________________
/// Format waypoints auxillary function
template<CoordType type>
vector<string> WayPoint::format2(const bool lat) const
{
//	fmt::print("@WayPoint::format2<CoordType::{}>(const bool lat) const; {}\n", type, lat? "lat" : "lon");
	auto& nv{ lat ? nvlat : nvlon };
	vector out_sv{ format0<type>(nv) };

	if constexpr (CoordType::decdeg != type)
		transform(out_sv.begin(), out_sv.end(), nv.begin(), out_sv.begin(), [lat](string& outstr, double n){return outstr + cardpoint(n < 0, lat);});
	return out_sv;
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
	if (!validated)
		revalid_Coord(nv);
	return validated >> 1;
}


/// __________________________________________________
/// Check "lat_valid" and "lon_valid attributes of DataFrame are all true
bool check_valid(const DataFrame df)
{
//	fmt::print("@check_valid(const DataFrame)\n");

	int latvalidated = check_logical_attr(df, "validlat");
	if (!latvalidated)
		return revalid_WayPoint(df);

	int lonvalidated = check_logical_attr(df, "validlon");
	if (!lonvalidated)
		return revalid_WayPoint(df);

	if (!(latvalidated >> 1))
		warning("Invalid latitude!");
	if (!(lonvalidated >> 1))
		warning("Invalid longitude!");
	return latvalidated >> 1 || lonvalidated >> 1;
}


/// __________________________________________________
/// Revalidate NumericVector or DataFrame
template<NumericVector_or_DataFrame T, Coord_or_WayPoint U>
const T revalidate(const T t)
{
//	fmt::print("@{} T: {}\n", "revalidate<T, U>(const T)", demangle(typeid(t)));
	warning("Revalidating %s…!", demangle(typeid(t)));
	validate<T, U>(t);	
	return check_valid(t);
}


/// __________________________________________________
/// Validate NumericVector or DataFrame
template<NumericVector_or_DataFrame T, Coord_or_WayPoint U>
inline const T validate(const T t)
{
//	fmt::print("@{} T: {}\n", "validate<T, U>(const T)", demangle(typeid(t)));
	U(get_coordtype(t), t).validate();
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
	convert_switch<NumericVector, Coord>(object, get_coordtype(fmt));
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
	if (newtype == type) {
//		fmt::print("——fmt out == fmt in!——\n"); std::fflush(nullptr);
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
//	fmt::print("{1}@{0}\n", "latlon(NumericVector, LogicalVector)", exportstr);
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
NumericVector validatecoords(NumericVector x, bool force = true)
{
//	fmt::print("{1}@{0} force: {2}\n", "validatecoords(NumericVector, bool)", exportstr, force);
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
	vector sv{ format_switch(Coord(ct, x), fmt ? get_coordtype(fmt) : ct) };
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
//	fmt::print("{1}@{0} fmt={2}\n", "as_waypoints(DataFrame, int)", exportstr, fmt);
	CoordType type = get_coordtype(fmt);
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
	convert_switch<DataFrame, WayPoint>(object, type);
	object.attr("class") = CharacterVector{"waypoints", "data.frame"};
	return object;
}


/// __________________________________________________
/// Convert waypoints format
//' @rdname convert
// [[Rcpp::export(name = "convert.waypoints")]]
DataFrame convertwaypoints(DataFrame x, int fmt)
{
	checkinherits(x, "waypoints");
	CoordType type = get_coordtype(x);
	CoordType newtype = get_coordtype(fmt);
//	fmt::print("{1}@{0} from {2} to {3}\n", "convertwaypoints(DataFrame, int)", exportstr, type, newtype);
	if (newtype == type) {
//		fmt::print("——fmt out == fmt in!——\n"); std::fflush(nullptr);
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
DataFrame validatewaypoints(DataFrame x, bool force = true)
{
//	fmt::print("{1}@{0} force: {2}\n", "validatewaypoints(DataFrame, bool)", exportstr, force);
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
CharacterVector formatwaypoints(DataFrame x, bool usenames = true, bool validate = true, int fmt = 0)
{
//	fmt::print("{1}@{0} usenames: {2}, validate: {3}\n", "formatwaypoints(DataFrame, bool, bool)", exportstr, usenames, validate);
	checkinherits(x, "waypoints");
	if(!x.nrows())
		stop("x has 0 rows!");
	if(!valid_ll(x))
		stop("Invalid llcols attribute!");
	if (validate)
		if (!check_valid(x))
			warning("Formatting invalid waypoints!");
	CoordType ct { get_coordtype(x) };
	vector sv{ format_switch(WayPoint(ct, x), fmt ? get_coordtype(fmt) : ct) };
	if (usenames) {
		RObject names = getnames(x);
		if (!prefixwithnames(sv, names))
			stop("Invalid \"namescol\" attribute!");
	}
	return wrap(sv);
}


/// __________________________________________________
/// Latitude and longitude headers for S3 print.waypoint()
//' @rdname format
// [[Rcpp::export]]
CharacterVector ll_headers(int width, int fmt)
{
//	fmt::print("{1}@{0} width={2}, fmt={3}\n", "ll_headers(int, int)", exportstr, width, fmt);
	--fmt;  //      to C++ array numbering
	constexpr int spacing[][3] { {15,  17,  18}, {11, 13, 14} };
	return wrap(vector {
		fmt::format("{:>{}}{:>{}}", "Latitude", width - spacing[0][fmt], "Longitude", spacing[0][fmt] - 1), // --fmt —> C++ array numbering
		fmt::format("{:>{}}", string(spacing[1][fmt], '_') + string(2, ' ') + string(spacing[1][fmt] + 1, '_'), width),
	});
}


/// __________________________________________________
/// Clone coords object from waypoints vector
//' @rdname coords
// [[Rcpp::export(name = "as_coords.waypoints")]]
NumericVector as_coordswaypoints(DataFrame object, bool which)
{
//	fmt::print("{1}@{0} which: {2}\n", "as_coord(DataFrame)", exportstr, which ? "lat" : "lon");
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