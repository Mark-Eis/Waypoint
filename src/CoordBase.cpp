#include <Rcpp.h>
#include <cxxabi.h>
using namespace Rcpp;

using std::vector; 
using std::string;
using std::transform;
using std::ostream;

#include "CoordBase.h"

#define FMT_HEADER_ONLY
// #include </opt/homebrew/Cellar/fmt/11.1.4/include/fmt/format.h>		// verbose path not found!
// #include </opt/homebrew/Cellar/fmt/11.1.4/include/fmt/ranges.h>		// verbose path not found!
#include <fmt/format.h>		// …/fmt/11.1.4/include/fmt/ copied to /Library/R/arm64/4.5/library/Rcpp/include
#include <fmt/ranges.h>		// …/fmt/11.1.4/include/fmt/ copied to /Library/R/arm64/4.5/library/Rcpp/include


/// __________________________________________________
/// __________________________________________________
/// Development and Debugging functions

/// Report object construction and destruction
void _ctrsgn(const std::type_info& obj, bool destruct)
{
//	fmt::print("{}ing ", destruct ? "Destroy" : "Construct");
//	std::fflush(nullptr);
	string s = obj.name();
	system(("c++filt -t " + s).data());
}

/// Demangle object names
const string demangle(const std::type_info& obj)
{
	int status = 0;
	char* p { abi::__cxa_demangle(obj.name(), NULL, NULL, &status) };
	string str { fmt::format("\"{}\" (status {})", p, std::to_string(status)) };
	std::free(p);
	return str;
}

/// Format string for debugging code
const auto& exportstr { "——Rcpp::export——" };

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
//	fmt::print("@{} attr=\"{}\" {}\n", "get_vec_attr<T, U>(const T&, const char*)", attrname, t.hasAttribute(attrname) ? true : false);
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	return t.hasAttribute(attrname) ? as<vector<U>>(t.attr(attrname)) : vector<U>();
}


/// __________________________________________________
/// Return "fmt" attribute as int
template<class T>
inline int get_fmt_attribute(const T& t)
{
//	fmt::print("@{} fmt={}\n", "get_fmt_attribute<T>(const T&)", as<int>(t.attr("fmt")));
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	return as<int>(t.attr("fmt"));
}


/// __________________________________________________
/// Does object inherit given class?
template<class T>
inline void checkinherits(T& t, const char* classname)
{
//	fmt::print("@{} T {} classname \"{}\"\n", "checkinherits<T>(T&, const char*)", demangle(typeid(t)), classname);
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
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
template<class T>
int nameinobj(const T t, const char* name)
{
//	fmt::print("@{} name={}\n", "nameinobj<T>(const T, const char*)", name);
	static_assert(std::is_same<List, T>::value || std::is_same<DataFrame, T>::value, "T must be List or DataFrame");
	vector<string> names { get_vec_attr<T, string>(t, "names") };
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
 
auto fmt::formatter<CoordType>::format(CoordType ct, format_context& ctx) const
	-> format_context::iterator
{
	string_view name = "unknown";
	switch (ct) {
		case CoordType::decdeg:
			name = "DecDeg";
			break;

		case CoordType::degmin:
			name = "DegMin";
			break;

		case CoordType::degminsec:
			name = "DegMinSec";
			break;

			default:
				stop("fmt::formatter<CoordType>::format(CoordType, format_context&) my bad");
	}
	return formatter<string_view>::format(name, ctx);
}


/// __________________________________________________
/// Convert int to CoordType enum
inline const CoordType get_coordtype(int i)
{
//	fmt::print("@{} {}\n", "get_coordtype(int)" , i);
	if (i < 1 || i > 3)
		stop("\"fmt\" must be between 1 and 3");
	return vector<CoordType>{ CoordType::decdeg, CoordType::degmin, CoordType::degminsec }[i - 1];
}


/// __________________________________________________
/// Convert "fmt" attribute to CoordType enum
template<class T>
inline const CoordType get_coordtype(const T& t)
{
//	fmt::print("@{} t {}\n", "get_coordtype<T>(const T&)", demangle(typeid(t)));
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
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
template<class T, class U>
void convert_switch(T t, CoordType newtype)
{
	CoordType type = get_coordtype(t);
//	fmt::print("@{} T: {} oldtype: {} newtype: {}\n", "convert_switch<T&, U>(T, CoordType)", demangle(typeid(t)), type, newtype);
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	static_assert(std::is_same<Coord, U>::value || std::is_same<WayPoint, U>::value, "T must be Coord or WayPoint");
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
//	fmt::print("@{} T: {} CoordType::{}\n", "format_switch<T>(const T&)", demangle(typeid(t)), t.get_coordtype());
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
/// Coordinate derived class

Coord::Coord(CoordType ct, const NumericVector nv) :
	Coordbase(ct), nv(nv),
	latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") } //,
{
//	fmt::print("§{} {} ", "Coord::Coord(CoordType, const NumericVector)", ct); _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Convert NumericVector coordinate format
template<CoordType newtype>
inline void Coord::convert() const
{
//	fmt::print("@Coord::convert<{}>() to {}\n", ct, newtype);
	transform(nv.begin(), nv.end(), const_cast<NumericVector&>(nv).begin(), Convertor<newtype>(ff));
}


/// __________________________________________________
/// Validate coords vector
void Coord::validate(bool warn) const
{
//	fmt::print("@{} latlon: {}\n", "Coord::validate()", fmt::join(latlon, ", "));
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
//	fmt::print("@Coord::format_ct<CoordType::{}>()\n", type);
	vector<string> out(nv.size());
	transform(nv.begin(), nv.end(), out.begin(), Format<type>(ff));
	transform(out.begin(), out.end(), nv.begin(), out.begin(), FormatLL<Coord, type>(ff, latlon));
	return out;
}


/// __________________________________________________
/// Format coords vector<string> with names
vector<string> Coord::format(bool usenames) const
{
//	fmt::print("@{}\n", "Coord::format(bool)");
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
//	fmt::print("§{} {} ", "WayPoint::WayPoint(WayPointType, const DataFrame)", ct); _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Convert DataFrame coordinate format
template<CoordType newtype>
inline void WayPoint::convert() const
{
//	fmt::print("@WayPoint::convert<{}>() to {}\n", ct, newtype);
	transform(nvlat.begin(), nvlat.end(), const_cast<NumericVector&>(nvlat).begin(), Convertor<newtype>(ff));
	transform(nvlon.begin(), nvlon.end(), const_cast<NumericVector&>(nvlon).begin(), Convertor<newtype>(ff));
}


/// __________________________________________________
/// Validate WayPoint
void WayPoint::validate(bool warn) const
{
//	fmt::print("@{}\n", "WayPoint::validate(bool)");

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
//	fmt::print("@WayPoint::format_ct<CoordType::{}>()\n", type);
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
//	fmt::print("@{}\n", "WayPoint::format(bool)");
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
//	fmt::print("@{}\n", "check_valid(const NumericVector)");
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
//	fmt::print("@{}\n", "check_valid(const DataFrame)");
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
//	fmt::print("@{} T: {} attrname: {} \n", "validated<T>(T, const char*, bool&)", demangle(typeid(t)), attrname);
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
//	fmt::print("@{} T: {}\n", "revalidate<T, U>(const T)", demangle(typeid(t)));
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	static_assert(std::is_same<Coord, U>::value || std::is_same<WayPoint, U>::value, "T must be Coord or WayPoint");
	warning("Revalidating %s…!", demangle(typeid(t)));
	validate<T, U>(t);	
	return check_valid(t);
}


/// __________________________________________________
/// Validate NumericVector or DataFrame
template<class T, class U>
inline const T validate(const T t)
{
//	fmt::print("@{} T: {}\n", "validate<T, U>(const T)", demangle(typeid(t)));
	static_assert(std::is_same<NumericVector, T>::value || std::is_same<DataFrame, T>::value, "T must be NumericVector or DataFrame");
	static_assert(std::is_same<Coord, U>::value || std::is_same<WayPoint, U>::value, "T must be Coord or WayPoint");
	U(get_coordtype(t), t).validate();
	return t;	
}


/// __________________________________________________
/// Check df has valid "llcols" attribute
bool valid_ll(const DataFrame df)
{
//	fmt::print("@{}\n", "valid_ll(const DataFrame)");
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
//		fmt::print("——fmt out == fmt in!——\n");
//		std::fflush(nullptr);
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
CharacterVector formatcoords(NumericVector x, bool usenames = true, bool validate = true)
{
//	fmt::print("{1}@{0} usenames: {2}, validate: {3}\n", "formatcoords(NumericVector, bool, bool)", exportstr, usenames, validate);
	checkinherits(x, "coords");
	if(!x.size())
		stop("x has 0 length!");
	if (validate)
		if (!check_valid(x))
			warning("Formatting invalid coords!");
	return wrap(Coord(get_coordtype(x), x).format(usenames));
}


/// __________________________________________________
/// Create waypoints
//' @rdname waypoints
// [[Rcpp::export(name = "as_waypoints.default")]]
DataFrame as_waypoints(DataFrame object, int fmt = 1)
{
//	fmt::print("{1}@{0} fmt={2}\n", "as_waypoints(DataFrame, int)", exportstr, fmt);
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
DataFrame convertwaypoints(DataFrame x, int fmt)
{
	checkinherits(x, "waypoints");
	CoordType type = get_coordtype(x);
	CoordType newtype = get_coordtype(fmt);
//	fmt::print("{1}@{0} from {2} to {3}\n", "convertwaypoints(DataFrame, int)", exportstr, type, newtype);
	if (newtype == type) {
//		fmt::print("——fmt out == fmt in!——\n");
//		std::fflush(nullptr);
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
CharacterVector formatwaypoints(DataFrame x, bool usenames = true, bool validate = true)
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
	return wrap(WayPoint(get_coordtype(x), x).format(usenames));
}


/// __________________________________________________
/// Latitude and longitude headers for S3 print.waypoint()
//' @rdname format
// [[Rcpp::export]]
CharacterVector ll_headers(const CharacterVector aswidth, int fmt)
{
//	fmt::print("{1}@{0} fmt={2}\n", "ll_headers(int, int)", exportstr, fmt);
	constexpr int spacing[] { 5,  7,  8,
							 11, 13, 14 };
	vector<string> sv {
			string("Latitude") + string(spacing[fmt - 1], ' ') + "Longitude ",
			string(spacing[fmt + 2], '_') + string(2, ' ') + string(spacing[fmt + 2] + 1, '_')
		};

	constexpr int adjust[] = { 2, 6, 10 };
	int width = (as<vector<string>>(aswidth)[0]).size() - adjust[fmt - 1];
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