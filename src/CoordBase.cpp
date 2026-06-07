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
void _ctrsgn(const std::type_info& obj, bool destruct)
{ /*
*/	fmt::print("{}ing ", destruct ? "Destroy" : "Construct");
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
	return t.hasAttribute(attrname) ? as<vector<U>>(t.attr(attrname)) : vector<U>{};
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
	auto maxwdth = max_element(sv.begin(), sv.end(), [](const string& a, const string& b){ return a.size() < b.size(); })->size();
	transform(sv.begin(), sv.end(), sv.begin(), [maxwdth](const string& s) { return fmt::format("{:<{}}", s, maxwdth); });
}


/// __________________________________________________
/// Concatenate corresponding elements of two vector<string>, with separator; result in second vector<string>
inline void concat_vecstr_elmnts(const vector<string>& sv_a, vector<string>& sv_b, const string sep)
{
//	fmt::print("@concat_vecstr_elmnts(vector<string>&, const vector<string>&, sep = \" \")\n");
	transform(sv_a.begin(), sv_a.end(), sv_b.begin(), sv_b.begin(), [&sep](const string& str_a, const string& str_b) {
		return str_a + sep + str_b; }); 
}


/// __________________________________________________
/// Concatenate corresponding elements of vector<int> and vector<string>, with separator; result in vector<string>
inline void concat_vecstr_elmnts(const vector<int>& iv_a, vector<string>& sv_b, const string sep)
{
//	fmt::print("@concat_vecstr_elmnts(const vector<int>&, vector<string>&, sep = \" \")\n");
	transform(iv_a.begin(), iv_a.end(), sv_b.begin(), sv_b.begin(), [&sep](const int i, const string& str_b) {
		return (std::to_string(i)) + sep + str_b; }); 
}


/// __________________________________________________
/// Prefix vector<string> elements with elements of RObject
inline bool prefixwithnames(vector<string>& sv, RObject& namesobj)
{
//	fmt::print("@{}\n", "prefixwithnames(vector<string>&, RObject&)");
	if (is<CharacterVector>(namesobj)) {
		vector<string>&& names = as<vector<string>>(namesobj);
		stdlenstr(names);
		concat_vecstr_elmnts(names, sv);
	} else if(is<IntegerVector>(namesobj))
		concat_vecstr_elmnts(as<vector<int>>(namesobj), sv);
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
/// Coordlet class

/// __________________________________________________
/// Constructor of Coordlet
template<CoordType current_type>
Coordlet<current_type>::Coordlet(NumericVector nv) :
	nv{ nv },
	latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") }
{
//	fmt::print("§Coordlet<CoordType::{}>::Coordlet(NumericVector, bool); ", current_type);  _ctrsgn(typeid(*this));
}


/// __________________________________________________
/// Format Coordlet::nv as vector<string> of CoordType
template<CoordType current_type> template<CoordType required_type>
vector<string> Coordlet<current_type>::format() const
{
//	fmt::print("@Coordlet<CoordType::{}>::format<CoordType::{}>() const\n", current_type, required_type);
	auto out_sv = vector<string>(nv.size());

	if constexpr (isDecDeg_v<required_type>)
		transform(nv.begin(), nv.end(), out_sv.begin(), [this](double n){
				return fmt::format("{:>{}.{}f}\u00B0", ff.get_decdeg(n), 11, 6);
			});
	if constexpr (isDegMin_v<required_type>)
		transform(nv.begin(), nv.end(), out_sv.begin(), [this](double n){
				return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) +
					   fmt::format("{:0>{}.{}f}\u2032", fabs(ff.get_decmin(n)), 7, 4);
			});
	if constexpr (isDegMinSec_v<required_type>)
		transform(nv.begin(), nv.end(), out_sv.begin(), [this](double n){
				return fmt::format("{:>{}}\u00B0", abs(ff.get_deg(n)), 3) +
					   fmt::format("{:0>{}}\u2032", abs(ff.get_min(n)), 2) +
					   fmt::format("{:0>{}.{}f}\u2033", fabs(ff.get_sec(n)), 5, 2);
			});

	return out_sv;
}


/// __________________________________________________
/// Switch CoordType required for Coordlet<CoordType>::format()
template<CoordType current_type>
vector<string> Coordlet<current_type>::format_switch(CoordType required_type) const
{
//	fmt::print("@Coordlet<CoordType::{}>::format_switch(CoordType); required: {}\n", current_type, required_type);

	using enum CoordType;
	switch (required_type)
	{
		case decdeg:
			return format<decdeg>();

		case degmin:
			return format<degmin>();

		case degminsec:
			return format<degminsec>();

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

	if constexpr (isDecDeg_v<required_type>)
		transform(nv.begin(), nv.end(), nv.begin(), [this](double n){
				return ff.get_decdeg(n);
			});
	if constexpr (isDegMin_v<required_type>)
		transform(nv.begin(), nv.end(), nv.begin(), [this](double n){
				return ff.get_deg(n) * 1e2 + ff.get_decmin(n);
			});
	if constexpr (isDegMinSec_v<required_type>)
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

	using enum CoordType;
	switch (required_type)
	{
		case decdeg:
			convert0<decdeg>();
			break;

		case degmin:
			convert0<degmin>();
			break;

		case degminsec:
			convert0<degminsec>();
			break;

		default:
			stop("Coordlet<CoordType>::convert_switch(CoordType) my bad");
	}
}


/// __________________________________________________
/// Validate Coordlet::nv
template<CoordType current_type>
const vector<bool> Coordlet<current_type>::validate() const
{
//	fmt::print("@Coordlet<CoordType::{}>::validate(); latlon: {}\n", current_type, fmt::join(latlon, ", "));
	vector<bool>::const_iterator ll_it{ latlon.begin() };
	auto ll_size { latlon.size() };
	auto valid = vector<bool>{};
	valid.assign(nv.size(), {false});

	transform(nv.begin(), nv.end(), valid.begin(), [this, &ll_it, &ll_size](double n){
		return !((fabs(ff.get_decdeg(n)) > (ll_size && (ll_size > 1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(fabs(ff.get_decmin(n)) >= 60) ||
				(fabs(ff.get_sec(n)) >= 60));
	});

	if (all_of(valid.begin(), valid.end(), [](bool v) { return v;}))
		valid.assign({true});

	return valid;
}


/// __________________________________________________
/// __________________________________________________
/// CrdWptBase class

/// __________________________________________________
/// Destructor
CrdWptBase::~CrdWptBase() {}

/// __________________________________________________
/// Switch current CoordType to convert nv
void CrdWptBase::convert_switch_current(NumericVector nv, const CoordType required_type) const
{
//	fmt::print("@CrdWptBase::convert_switch_current(NumericVector, const CoordType); current_type: {}; required_type: {}\n", ct, required_type);
	using enum CoordType;
	switch (ct)
	{
		case decdeg:
			convert_dispatch<decdeg>(nv, required_type);
			break;

		case degmin:
			convert_dispatch<degmin>(nv, required_type);
			break;

		case degminsec:
			convert_dispatch<degminsec>(nv, required_type);
			break;

		default:
			stop("CrdWptBase::convert_switch_current(NumericVector, const CoordType) const my bad");
	}
}

/// __________________________________________________
/// Switch current CoordType to format nv
vector<string> CrdWptBase::format_switch_current(NumericVector nv, const CoordType required_type) const
{
//	fmt::print("@CrdWptBase::format_switch_current(NumericVector, const CoordType, bool); current: {}; required: {}\n", ct, required_type);
	using enum CoordType;
	switch (ct)
	{
		case decdeg:
			return format_dispatch<decdeg>(nv, required_type);

		case degmin:
			return format_dispatch<degmin>(nv, required_type);

		case degminsec:
			return format_dispatch<degminsec>(nv, required_type);

		default:
			stop("CrdWptBase::format_switch_current(NumericVector, const CoordType) const my bad");
	}
}


/// __________________________________________________
/// Switch current CoordType to validate nv
const vector<bool> CrdWptBase::validate_switch_current(const NumericVector nv) const
{
//	fmt::print("@CrdWptBase::validate_switch_current(const NumericVector); current_type: {}\n", ct);
	using enum CoordType;
	switch (ct)
	{
		case decdeg:
			return validate_dispatch<decdeg>(nv);

		case degmin:
			return validate_dispatch<degmin>(nv);

		case degminsec:
			return validate_dispatch<degminsec>(nv);

		default:
			stop("CrdWptBase::validate_switch_current(const NumericVector) const my bad");
	}
}

/// __________________________________________________
/// Dispatch nv to Coordlet<CoordType>::format_switch()
template<CoordType current_type> 
inline vector<string> CrdWptBase::format_dispatch(NumericVector nv, const CoordType required_type) const
{
//	fmt::print("@CrdWptBase::format_dispatch<CoordType::{}>(NumericVector, const CoordType); required: {}\n", current_type, required_type);
	return Coordlet<current_type>{ nv }.format_switch(required_type);
}


/// __________________________________________________
/// Dispatch nv to Coordlet<CoordType>::convert_switch()
template<CoordType current_type> 
inline void CrdWptBase::convert_dispatch(NumericVector nv, const CoordType new_type) const
{
//	fmt::print("@CrdWptBase::convert_dispatch<CoordType::{}>(NumericVector, const CoordType); new type: {}\n", current_type, new_type);
	Coordlet<current_type>{ nv }.convert_switch(new_type);
}


/// __________________________________________________
/// Dispatch nv to Coordlet<CoordType>::validate()
template<CoordType current_type> 
inline const vector<bool> CrdWptBase::validate_dispatch(const NumericVector nv) const
{
//	fmt::print("@CrdWptBase::validate_dispatch<CoordType::{}>(const NumericVector)\n", current_type);
	return Coordlet<current_type>{ nv }.validate();
}


/// __________________________________________________
/// __________________________________________________
/// Coords class

Coords::Coords(NumericVector nv) : CrdWptBase { get_coordtype(nv) }, nv{ nv }
{
//	fmt::print("§Coords::Coords(NumericVector); {} ", ct); _ctrsgn(typeid(*this));
}


void Coords::convert(CoordType newtype)
{
//	fmt::print("@Coords::convert(CoordType); current type: {}; new type: {}\n", ct, newtype);
	convert_switch_current(nv, newtype);
}


vector<string> Coords::format(CoordType required_type) const
{
//	fmt::print("@Coords::format(CoordType); current type: {}; required type: {}\n", ct, required_type);
	vector sv_out{ format_switch_current(nv, required_type) };
	format_suffix_switch(sv_out, required_type);
	return sv_out;
}


/// __________________________________________________
/// Switch current CoordType to format suffix
void Coords::format_suffix_switch(vector<string>& sv_out, const CoordType required_type) const
{
//	fmt::print("@Coords::format_suffix_switch(vector<string>&, const CoordType); current: {}; required: {}\n", ct, required_type);
	using enum CoordType;
	switch (required_type)
	{
		case decdeg:
			format_suffix<decdeg>(sv_out);
			break;

		case degmin:
			format_suffix<degmin>(sv_out);
			break;

		case degminsec:
			format_suffix<degminsec>(sv_out);
			break;

		default:
			stop("Coords::format_suffix_switch(vector<string>&, const CoordType) const my bad");
	}
}


template<CoordType required_type>
void Coords::format_suffix(vector<string>& sv_out) const
{
//	fmt::print("@Coords:: format_suffix(vector<string>& sv_out, CoordType) const; required type: {}\n", required_type);
	const auto latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") };
	vector<bool>::const_iterator ll_it { latlon.begin() };
	const auto ll_size { latlon.size() };

	if constexpr(isDecDeg_v<required_type>) {
		const auto lambda1 = [&ll_it](string& outstr, double n){ return outstr + (*ll_it++ ? " lat" : " lon"); };
		const auto lambda2 = [&ll_it](string& outstr, double n){ return outstr + (*ll_it ? " lat" : " lon"); };

		if (ll_size > 1)
			transform(sv_out.begin(), sv_out.end(), nv.begin(), sv_out.begin(), lambda1);
		else
			if (ll_size == 1)	// uniform coords
				transform(sv_out.begin(), sv_out.end(), nv.begin(), sv_out.begin(), lambda2);

	}

	if constexpr(isDegMin_v<required_type> || isDegMinSec_v<required_type>) {
		const auto lambda1 = [&ll_it](string& outstr, double n){ return outstr + cardpoint(n < 0, *ll_it++); };
		const auto lambda2 = [&ll_it](string& outstr, double n){ return outstr + cardpoint(n < 0, *ll_it); };
		const auto lambda3 = [](string& outstr, double n){ return outstr + cardi_b(n < 0); };

		if (ll_size > 1)
			transform(sv_out.begin(), sv_out.end(), nv.begin(), sv_out.begin(), lambda1);
		else
			if (ll_size == 1)	// uniform coords
				transform(sv_out.begin(), sv_out.end(), nv.begin(), sv_out.begin(), lambda2);
			else					// no latlon info
				transform(sv_out.begin(), sv_out.end(), nv.begin(), sv_out.begin(), lambda3);
	}
}


const bool Coords::validate() const
{
//	fmt::print("@Coords::validate(); current type: {}\n", ct);

	auto valid = validate_switch_current(nv);
	static_cast<NumericVector>(nv).attr("valid") = valid;
	return ( std::all_of(valid.begin(), valid.end(), [](bool i){ return i; } )
	);
}


/// __________________________________________________
/// __________________________________________________
/// Waypoints class

Waypoints::Waypoints(DataFrame df) :
	CrdWptBase { get_coordtype(df) }, df{ df },
	nvlat( df[get_vec_attr<DataFrame, int>(df, "llcols")[0] - 1] ), 
	nvlon( df[get_vec_attr<DataFrame, int>(df, "llcols")[1] - 1] )
{
//	fmt::print("§Waypoints::Waypoints(DataFrame); {} ", ct); _ctrsgn(typeid(*this));
	nvlat.attr("fmt") = coordtype_to_int(ct) + 1;
	nvlon.attr("fmt") = coordtype_to_int(ct) + 1;
	nvlat.attr("latlon") = true;
	nvlon.attr("latlon") = false;
}


Waypoints::~Waypoints()
{
//	fmt::print("§Waypoints::~Waypoints(); {} ", ct); _ctrsgn(typeid(*this), true);
	nvlat.attr("latlon") = R_NilValue;
	nvlon.attr("latlon") = R_NilValue;
	nvlat.attr("fmt") = R_NilValue;
	nvlon.attr("fmt") = R_NilValue;
}


void Waypoints::convert(CoordType newtype)
{
//	fmt::print("@ Waypoints::convert(CoordType); current type: {}; new type: {}\n", ct, newtype);
	convert_switch_current(nvlat, newtype);
	convert_switch_current(nvlon, newtype);
}


vector<string> Waypoints::format(CoordType required_type) const
{
//	fmt::print("@Waypoints::format(CoordType); current type: {}; required type: {}\n", ct, required_type);

	vector sv_lat{ format_switch_current(nvlat, required_type) };
	vector sv_lon{ format_switch_current(nvlon, required_type) };

	if (required_type != CoordType::decdeg) {
		format_suffix(sv_lat, true);
		format_suffix(sv_lon, false);
	}

	transform(sv_lat.begin(), sv_lat.end(), sv_lon.begin(), sv_lat.begin(), [](auto& latstr, auto& lonstr){return latstr + "  " + lonstr;});
	return sv_lat;
}


void Waypoints::format_suffix(vector<string>& out_sv, const bool latlon) const
{
//	fmt::print("@Waypoints::format_suffix(vector<string> out_sv) const; {}\n", latlon? "lat" : "lon");
	transform(out_sv.begin(), out_sv.end(), (latlon? nvlat : nvlon).begin(), out_sv.begin(), [latlon](string& outstr, double n){
		return outstr + cardpoint(n < 0, latlon); }
	);
}


const bool Waypoints::validate() const
{
//	fmt::print("@Waypoints::validate(); current type: {}\n", ct);

	auto validlat = validate_switch_current(nvlat);
	auto validlon = validate_switch_current(nvlon);

	static_cast<DataFrame>(df).attr("validlat") = validlat;
	static_cast<DataFrame>(df).attr("validlon") = validlon;

	return (
		std::all_of(validlat.begin(), validlat.end(), [](bool i){ return i; }) &&
		std::all_of(validlon.begin(), validlon.end(), [](bool i){ return i; })
	);
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
		return revalidate<NumericVector, Coords>(nv);
	return validated >> 1;
}


/// __________________________________________________
/// Check "lat_valid" and "lon_valid attributes of DataFrame are all true
bool check_valid(const DataFrame df)
{
//	fmt::print("@check_valid(const DataFrame)\n");

	int latvalidated = check_logical_attr(df, "validlat");
	int lonvalidated = check_logical_attr(df, "validlon");

	if (!(latvalidated & lonvalidated))
		return revalidate<DataFrame, Waypoints>(df);

	if (!(latvalidated >> 1))
		warning("Invalid latitude!");
	if (!(lonvalidated >> 1))
		warning("Invalid longitude!");
	return latvalidated >> 1 || lonvalidated >> 1;
}


/// __________________________________________________
/// Revalidate "coords" or "waypoints"
template<NumericVector_or_DataFrame T, Coords_or_Waypoints U>
bool revalidate(const T t)
{
//	fmt::print("@revalidate<NumericVector_or_DataFrame, Coords_or_Waypoints>(const T); T: {}; U: {}\n", demangle(typeid(t)), demangle(typeid(U)));
	warning("Revalidating %s…!", demangle(typeid(t)));
	if (!U{ t }.validate())
	    warning("Invalid coords or waypoints! [revalidate]");
	return check_valid(t);
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
/// Create coords - S3 method as_coords.default()
//' @rdname coords 
// [[Rcpp::export(name = "as_coords.default")]]
NumericVector as_coords(NumericVector object, int fmt = 1)
{
//	fmt::print("{1}@{0} fmt={2}\n", "as_coords(NumericVector, int)", exportstr, fmt);
	object.attr("fmt") = fmt;
	Coords{ object }.validate();
	object.attr("class") = "coords";
	return object;
}


/// __________________________________________________
/// Convert coords - S3 method convert.coords()
//' @rdname convert
// [[Rcpp::export(name = "convert.coords")]]
NumericVector convertcoords(NumericVector x, int fmt)
{
	checkinherits(x, "coords");
	CoordType type = get_coordtype(x);
	CoordType newtype = get_coordtype(fmt);
//	fmt::print("{}@convertcoords(NumericVector, int); from {} to {}\n", exportstr, type, newtype);
	if (!check_valid(x))
		stop("Invalid coords!");
	if (newtype != type) {
		Coords{ x }.convert(newtype);
		x.attr("fmt") = fmt;
	} else
		Rcout << "\t—— fmt out == fmt in! ——\n\n";
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
	Coords{ cd }.validate();
	return cd;
}


/// __________________________________________________
/// Validate coords - S3 method validate.coords()
//' @rdname validate
// [[Rcpp::export(name = "validate.coords")]]
NumericVector validatecoords(const NumericVector x, const bool force = true)
{
//	fmt::print("{1}@{0} force: {2}\n", "validatecoords(const NumericVector, const bool)", exportstr, force);
	checkinherits(x, "coords");
	if (force)									
		Coords{ x }.validate();
	if (!check_valid(x))
		warning("Invalid coords!");
	return x;
}


/// __________________________________________________
/// Format coords - S3 method format.coords()
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

	vector sv{ Coords{ x }.format(fmt ? get_coordtype(fmt) : get_coordtype(x)) };

	vector names{ get_vec_attr<NumericVector, string>(x, "names") };
	if (names.size() && usenames) {
		stdlenstr(names);
		concat_vecstr_elmnts(names, sv);
	}
	return wrap(sv);
}


/// __________________________________________________
/// Create waypoints - S3 method as_waypoints.default()
//' @rdname waypoints
// [[Rcpp::export(name = "as_waypoints.default")]]
DataFrame as_waypoints(DataFrame object, int fmt = 1)
{
//	fmt::print("{1}@{0} fmt={2}\n", "as_waypoints(DataFrame, int)", exportstr, fmt);
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
	Waypoints{ object }.validate();
	object.attr("class") = CharacterVector{"waypoints", "data.frame"};
	return object;
}


/// __________________________________________________
/// Convert waypoints type - S3 method convert.waypoints()
//' @rdname convert
// [[Rcpp::export(name = "convert.waypoints")]]
DataFrame convertwaypoints(DataFrame x, int fmt)
{
	checkinherits(x, "waypoints");
	CoordType type = get_coordtype(x);
	CoordType newtype = get_coordtype(fmt);
//	fmt::print("{}@convertwaypoints(DataFrame, int); from {} to {}\n", exportstr, type, newtype);
	if (!check_valid(x))
		stop("Invalid waypoints!");
	if(!valid_ll(x))
		stop("Invalid llcols attribute!");
	if (newtype != type) {
		Waypoints{ x }.convert(newtype);
		x.attr("fmt") = fmt;
	} else
		Rcout << "\t—— fmt out == fmt in! ——\n\n";
	return x;
}


/// __________________________________________________
/// Validate waypoints - S3 method validate.waypoints()
//' @rdname validate
// [[Rcpp::export(name = "validate.waypoints")]]
DataFrame validatewaypoints(DataFrame x, bool force = true)
{
//	fmt::print("{1}@{0} force: {2}\n", "validatewaypoints(DataFrame, bool)", exportstr, force);
	checkinherits(x, "waypoints");
	if(!valid_ll(x))
		stop("Invalid llcols attribute!");
	if (force)
		Waypoints{ x }.validate();
	if (!check_valid(x))
		warning("Invalid waypoints!");
	return x;
}


/// __________________________________________________
/// Format waypoints - S3 method format.waypoints()
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
	vector sv{ Waypoints{ x }.format(fmt ? get_coordtype(fmt) : get_coordtype(x)) };

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
