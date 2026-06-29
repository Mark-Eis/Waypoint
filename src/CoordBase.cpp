/// __________________________________________________
/// CoordBase.cpp
/// __________________________________________________

// [[Rcpp::plugins(cpp23)]]

#include <array>
#include <Rcpp.h>
#include <cxxabi.h>
#include <memory>
#include <string>
#include <utility>

using namespace Rcpp;

using std::array;
using std::vector;
using std::string;
using namespace std::literals;
using std::string_view;
using namespace std::string_view_literals;
using std::transform;
using std::unique_ptr;
using std::make_unique;

#include "CoordBase.h"

#define FMT_HEADER_ONLY
#include "fmt/format.h"		// …fmt/*.h copied to ~/Documents/R/Packages/Waypoint/src/fmt.
#include "fmt/ranges.h"		// …fmt/*.h copied to ~/Documents/R/Packages/Waypoint/src/fmt.


/// __________________________________________________
/// __________________________________________________
/// Development and Debugging functions

#if DEBUG > 0

/// __________________________________________________
/// Report object construction and destruction
void _ctrsgn(const std::type_info& obj, bool construct)
{ /*
*/	fmt::print("{}structing: ", construct ? "§§§Con" : "~§§De");
	std::fflush(nullptr);
	system(("c++filt -t " + string{ obj.name() }).data());
	std::fflush(nullptr);
}

#endif

/// __________________________________________________
/// Demangle object names
const string demangle(const std::type_info& obj)
{
	int status = 0;
	const auto objname { obj.name() };
//	char* p { abi::__cxa_demangle(obj.name(), NULL, NULL, &status) };
	char* p { abi::__cxa_demangle(objname, NULL, NULL, &status) };
//	string str { fmt::format("\"{}\" (status {})", p, std::to_string(status)) };
	string str { fmt::format("{}", p) };
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
template<NumVec_or_DataFrame T, typename U> 
inline vector<U> get_vec_attr(const T& t, const string attrname)
{
//	fmt::print("@{} attr=\"{}\" {}\n", "get_vec_attr<T, U>(const T&, const string)", attrname, t.hasAttribute(attrname) ? true : false);
	return t.hasAttribute(attrname) ? as<vector<U>>(t.attr(attrname)) : vector<U>{};
}

/// __________________________________________________
/// Return "fmt" attribute as int
inline int get_fmt_attribute(const NumVec_or_DataFrame auto& t)
{
//	fmt::print("@get_fmt_attribute<T>(const NumVec_or_DataFrame auto&); t: {}; fmt={}\n", demangle(typeid(t)), as<int>(t.attr("fmt")));
	return as<int>(t.attr("fmt"));
}

/// __________________________________________________
/// Check whether a NumericVector or DataFrame has a specified logical vector attribute and whether all true
int check_logical_attr(NumVec_or_DataFrame auto t, const string attrname)
{
//	fmt::print("@check_logical_attr(NumVec_or_DataFrame auto, const string); T: {}; attrname {}\n", demangle(typeid(t)), attrname);
	const vector vec_attr{ get_vec_attr<decltype(t), bool>(t, attrname) };
	if (vec_attr.size()) {
		return all_of(vec_attr.begin(), vec_attr.end(), [](bool v) { return v;}) ? 0b11 : 0b01;
	} else {
		return 0b00;
	}
}

/// __________________________________________________
/// Does object inherit given class?
inline void checkinherits(const NumVec_or_DataFrame auto& t, const string classname)
{
//	fmt::print("@checkinherits(const NumVec_or_DataFrame auto&, const string); t: {}; classname \"{}\"\n", demangle(typeid(t)), classname);
	if (!t.inherits(classname.c_str())) stop("Argument must be a \"%s\" object", classname.c_str());
}

/// __________________________________________________
/// Is item number present in data.frame? (Using C++ numbering)
inline bool is_item_in_df(const DataFrame df, int item_no)
{
//	fmt::print("@is_item_in_df(const DataFrame, int); item no. {}\n", item_no);
	if (NA_INTEGER == item_no)
		return false;
	else
		return !(item_no < 0) && item_no < df.size();
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
/// Find position of name within data.frame names
int name_pos_in_df(const DataFrame df, const string name)
{
//	fmt::print("@name_pos_in_df(const DataFrame, const string); name={}\n", name);
	vector names{ get_vec_attr<DataFrame, string>(df, "names"s) };
	if (!names.size())
		return -1;
	typedef decltype(names.size()) Tmp;
	Tmp i = 0;
	for (auto str : names ) {
		// fmt::print("@@name_pos_in_df(const DataFrame, const string); testing: {}\n", str);
		if (!str_tolower(str).compare(name)) {
			// fmt::print("@@@name_pos_in_df(const DataFrame, const string); found: {}\n", str);
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
	vector namescolvec{ get_vec_attr<DataFrame, int>(df, "namescol"s) };
	if (1 == namescolvec.size()) {
		int namescol = namescolvec[0] - 1;
		if (is_item_in_df(df, namescol))
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
	constexpr array names {"DecDeg"sv, "DegMin"sv, "DegMinSec"sv};
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
inline const CoordType get_coordtype(const NumVec_or_DataFrame auto& t)
{
//	fmt::print("@get_coordtype(const NumVec_or_DataFrame auto&); t: {}\n", demangle(typeid(t)));
	return get_coordtype(get_fmt_attribute(t));
}

/// __________________________________________________
/// Convert CoordType enum to int; + 1 for R
inline int coordtype_to_int(CoordType ct)
{
//	fmt::print("@{} ct={}\n", "coordtype_to_int(CoordType)", ct);
	return static_cast<char>(ct) + 1;
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
/// Coords class —— Constructor
template<DVecType T>
Coords<T>::Coords(T t, const vector<bool> latlon) :
	dv { std::move(t) },
	latlon { latlon } //,
{
#if DEBUG > 0
	_ctrsgn(typeid(*this)); fmt::print("\t(T, const vector<bool>)\n");
#endif
}

/// __________________________________________________
/// Format dv as a vectype object —— private
template<DVecType T> template<vectype U, functador V>
inline U Coords<T>::conform0() const
{
#if DEBUG > 0
	fmt::print("@Coords<T>::conform0<U, V>() const;\n\t\t\tT: {},\n\t\t\tU: {},\n\t\t\tV: {}\n",
		demangle(typeid(T)), demangle(typeid(U)), demangle(typeid(V)));
#endif
	U uv_out(dv.size());
	transform(dv.begin(), dv.end(), uv_out.begin(), V());	

	if constexpr (isDecDegVecString_v<U>)
		suffix_latlon(uv_out);
	else if constexpr (isDegMinVecString_v<U> || isDegMinSecVecString_v<U>)
		suffix_nesw(uv_out);
#if DEBUG > 0
	fmt::print("@ICoords<T>::conform0 <U, V>() const;\n\t{}\t &uv_out {}, &uv_out[0] {}, uv_out[0] {}, typeid: {}\n",
		padstr, address(uv_out), address(uv_out[0]), uv_out[0], demangle(typeid(uv_out)));
#endif
	return uv_out;
}

/// __________________________________________________
/// Add suffix of "N", "E", "S", "W"; or "(N/E)", "(S/W)"
template<DVecType T>
void Coords<T>::suffix_nesw(vector<string>& sv_out) const
{
#if DEBUG > 0
	fmt::print("@Coords<T>::suffix_nesw(vector<string>& sv_out) const\n");
#endif
	vector<bool>::const_iterator ll_it { latlon.begin() };
	const auto ll_size { latlon.size() };

	const auto lambda1 = [&ll_it](auto& outstr, auto n){ return outstr + cardpoint(n < 0, *ll_it++); };
	const auto lambda2 = [&ll_it](auto& outstr, auto n){ return outstr + cardpoint(n < 0, *ll_it); };
	const auto lambda3 = [](auto& outstr, auto n){ return outstr + cardi_b(n < 0); };

	if (ll_size > 1)
		transform(sv_out.begin(), sv_out.end(), dv.begin(), sv_out.begin(), lambda1);
	else
		if (ll_size == 1)	// uniform coords
			transform(sv_out.begin(), sv_out.end(), dv.begin(), sv_out.begin(), lambda2);
		else				// no latlon info
			transform(sv_out.begin(), sv_out.end(), dv.begin(), sv_out.begin(), lambda3);
}

/// __________________________________________________
/// Add suffix of "lat", "lon"
template<DVecType T>
void Coords<T>::suffix_latlon(vector<string>& sv_out) const
{
#if DEBUG > 0
	fmt::print("@Coords<T>::suffix_latlon(vector<string>& sv_out) const\n");
#endif
	vector<bool>::const_iterator ll_it { latlon.begin() };
	const auto ll_size { latlon.size() };

	const auto lambda1 = [&ll_it](auto& outstr, auto n){ return outstr + (*ll_it++ ? " lat" : " lon"); };
	const auto lambda2 = [&ll_it](auto& outstr, auto n){ return outstr + (*ll_it ? " lat" : " lon"); };

	if (ll_size > 1)
		transform(sv_out.begin(), sv_out.end(), dv.begin(), sv_out.begin(), lambda1);
	else
		if (ll_size == 1)	// uniform coords
			transform(sv_out.begin(), sv_out.end(), dv.begin(), sv_out.begin(), lambda2);
}

/// __________________________________________________
/// conform call entry point —— public
template<DVecType T> template<typename U, template <typename V> typename F>
vector<U> Coords<T>::conform(CoordType required) const
{
#if DEBUG > 0
	fmt::print("@Coords<T>::conform<U, F<V>>(CoordType) const;\n\t\t\tT: {},\n\t\t\tU: {},\n\t\t\tF: functador<T, vectype<U>>,\n\t\t\trequired: {}\n",
		demangle(typeid(T)), demangle(typeid(U)), required);
#endif
	using enum CoordType;
	switch (required)
	{
		case decdeg:
			return conform0<DecDegVec<U>, F<DecDegVec<U>>>();

		case degmin:
			return conform0<DegMinVec<U>, F<DegMinVec<U>>>();

		case degminsec:
			return conform0<DegMinSecVec<U>, F<DegMinSecVec<U>>>();

		default:
			stop("Coords<T>::conform<U, F<V>>(CoordType) const my bad");
	}
}

/// __________________________________________________
/// convert call entry point —— public												// ¡¡¡—— Deprecate ——!!!
template<DVecType T>
vector<double> Coords<T>::convert(CoordType newtype) const
{
#if DEBUG > 0
	fmt::print("@Coords<T>::convert(CoordType) const; T: {}, new type: {}\n", demangle(typeid(T)), newtype);
#endif
	using enum CoordType;
	switch (newtype)
	{
		case decdeg:
			return conform0<DecDegVecDouble, Convertidor<T, DecDegVecDouble>>();

		case degmin:
			return conform0<DegMinVecDouble, Convertidor<T, DegMinVecDouble>>();

		case degminsec:
			return conform0<DegMinSecVecDouble, Convertidor<T, DegMinSecVecDouble>>();

		default:
			stop("Coords<T>::convert(CoordType) const my bad");
	}
}

/// __________________________________________________
/// Format call entry point —— public												// ¡¡¡—— Deprecate ——!!!
template<DVecType T>
vector<string> Coords<T>::format(CoordType required_type) const
{
#if DEBUG > 0
	fmt::print("@Coords<T>::format(CoordType) const; T: {}, required type: {}\n", demangle(typeid(T)), required_type);
#endif
	using enum CoordType;
	switch (required_type)
	{
		case decdeg:
			return conform0<DecDegVecString, Formateador<T, DecDegVecString>>();

		case degmin:
			return conform0<DegMinVecString, Formateador<T, DegMinVecString>>();

		case degminsec:
			return conform0<DegMinSecVecString, Formateador<T, DegMinSecVecString>>();

		default:
			stop("Coords<T>::format(CoordType) const my bad");
	}
}

/// __________________________________________________
/// Validation call entry point —— public
template<DVecType T>
const vector<bool> Coords<T>::validate() const
{
#if DEBUG > 0
	fmt::print("@Coords<T>::validate(); latlon: {}\n", fmt::join(latlon, ", "));
	fmt::print("@ICoords<T>::validate();\n\t{}\t\t &dv {}, &dv[0] {}, dv[0] {}, typeid: {}\n",
		padstr, address(dv), address(dv[0]), dv[0], demangle(typeid(dv)));
#endif
	FamousFive<T> ff {};
	vector<bool>::const_iterator ll_it{ latlon.begin() };
	auto ll_size { latlon.size() };
	auto valid = vector<bool>{};
	valid.assign(dv.size(), {false});

	transform(dv.begin(), dv.end(), valid.begin(), [&ff, &ll_it, &ll_size](auto n){
		return !((fabs(ff.get_decdeg(n)) > (ll_size && (ll_size > 1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(fabs(ff.get_decmin(n)) >= 60) ||
				(fabs(ff.get_sec(n)) >= 60));
	});

	if (all_of(valid.begin(), valid.end(), [](auto v) { return v;}))
		valid.assign({true});

#if DEBUG > 0
	fmt::print("@IICoords<T>::validate();\n\t{}\t  &valid {}, typeid: {}, \n\t{}\n",
		padstr, address(valid), demangle(typeid(valid)), fmt::join(valid, ", "));
#endif
	return valid;
}

/// __________________________________________________
/// Instantiate Coords<T> object
template<DVecType T>
inline coords_t auto coordsmaker(NumericVector nv, vector<bool> latlon)
{
#if DEBUG > 0
	fmt::print("@coordsmaker(NumericVector, vector<bool>);\n\t{}\t\t &nv {}, &nv[0] {}, nv[0] {}, DVecType: {}\n",
		padstr, address(nv), address(nv[0]), nv[0], demangle(typeid(T)));
#endif
	if (!latlon.size())
		latlon = get_vec_attr<NumericVector, bool>(nv, "latlon"s);
	return Coords<T>(nv, latlon);
}

/// __________________________________________________
/// Convert "coords" NumericVector
vector<double> convert_switch(const NumericVector nv, CoordType newtype)
{
#if DEBUG > 0
	fmt::print("@convert_switch(const NumericVector, CoordType); current type: {}, new type: {}\n", get_coordtype(nv), newtype);
#endif
	using enum CoordType;
	switch (get_coordtype(nv))
	{
		case decdeg:
			return coordsmaker<DecDegVecDouble>(nv).conform<double, ConvertidorDecDegVec>(newtype);

		case degmin:
			return coordsmaker<DegMinVecDouble>(nv).conform<double, ConvertidorDegMinVec>(newtype);

		case degminsec:
			return coordsmaker<DegMinSecVecDouble>(nv).conform<double, ConvertidorDegMinSecVec>(newtype);

		default:
			stop("convert_switch(const NumericVector, CoordType) const my bad");
	}
}

/// __________________________________________________
/// Format "coords" NumericVector
vector<string> format_switch(const NumericVector nv, CoordType ct_required)
{
#if DEBUG > 0
	fmt::print("@format_switch(const NumericVector, CoordType); current type: {}, required type: {}\n", get_coordtype(nv), ct_required);
#endif
	using enum CoordType;
	switch (get_coordtype(nv))
	{
		case decdeg:
			return coordsmaker<DecDegVecDouble>(nv).conform<string, FormateadorDecDegVec>(ct_required);

		case degmin:
			return coordsmaker<DegMinVecDouble>(nv).conform<string, FormateadorDegMinVec>(ct_required);

		case degminsec:
			return coordsmaker<DegMinSecVecDouble>(nv).conform<string, FormateadorDegMinSecVec>(ct_required);

		default:
			stop("format_switch(const NumericVector, CoordType) const my bad");
	}
}

/// __________________________________________________
/// Validate "coords" NumericVector 
const vector<bool> validate_switch(const NumericVector nv)
{
#if DEBUG > 0
	fmt::print("@validate_switch(const NumericVector); current type: {}\n", get_coordtype(nv));
#endif
	using enum CoordType;
	switch (get_coordtype(nv))
	{
		case decdeg:
			return coordsmaker<DecDegVecDouble>(nv).validate();

		case degmin:
			return coordsmaker<DegMinVecDouble>(nv).validate();

		case degminsec:
			return coordsmaker<DegMinSecVecDouble>(nv).validate();

		default:
			stop("validate_switch(const NumericVector) const my bad");
	}
}


/// __________________________________________________
/// __________________________________________________
/// Waypoints class

/// __________________________________________________
/// Constructor
template<DVecType T>
Waypoints<T>::Waypoints(NumericVector nv_lat, NumericVector nv_lon) :
	crdlat{ coordsmaker<T>(nv_lat, vector<bool>{ true }) },
	crdlon{ coordsmaker<T>(nv_lon, vector<bool>{ false }) }
{
#if DEBUG > 0
	_ctrsgn(typeid(*this)); fmt::print("\t(vector<double>, vector<double>)\n");
#endif
}

/// __________________________________________________
/// convert call entry point -- public
template<DVecType T>
const bisvec<double> Waypoints<T>::convert(CoordType newtype) const
{
#if DEBUG > 0
	fmt::print("@Waypoints<T>::convert(CoordType) const; T: {}, new type: {}\n", demangle(typeid(T)), newtype);
#endif
	return { crdlat.convert(newtype),  crdlon.convert(newtype)};
}

/// __________________________________________________
/// format call entry point -- public
template<DVecType T>
const bisvec<string> Waypoints<T>::format(CoordType required_type) const
{
#if DEBUG > 0
	fmt::print("@Waypoints<T>::format(CoordType) const; T: {}, required type: {}\n", demangle(typeid(T)), required_type);
#endif
	return { crdlat.format(required_type),  crdlon.format(required_type)};
}

/// __________________________________________________
/// validate call entry point -- public
template<DVecType T>
const bisconstvec<bool> Waypoints<T>::validate() const
{
#if DEBUG > 0
	fmt::print("@Waypoints<T>::validate(CoordType) const; T: {}\n", demangle(typeid(T)));
#endif
	return { crdlat.validate(),  crdlon.validate()};
}


/// __________________________________________________
/// Instantiate Waypoints<T> object
template<DVecType T>
inline waypoints_t auto waypointsmaker(DataFrame df)
{
#if DEBUG > 0
	fmt::print("@waypointsmaker(DataFrame); DVecType: {}\n", demangle(typeid(T)));
#endif
	return Waypoints<T>(df[get_vec_attr<DataFrame, int>(df, "llcols")[0] - 1], df[get_vec_attr<DataFrame, int>(df, "llcols")[1] - 1]);
}

/// __________________________________________________
/// Convert "waypoints" DataFrame 
const bisvec<double> convert_switch(const DataFrame df, CoordType newtype)
{
#if DEBUG > 0
	fmt::print("@convert_switch(const DataFrame, CoordType); current type: {}, new type: {}\n", get_coordtype(df), newtype);
#endif
	using enum CoordType;
	switch (get_coordtype(df))
	{
		case decdeg:
			return waypointsmaker<DecDegVecDouble>(df).convert(newtype);

		case degmin:
			return waypointsmaker<DegMinVecDouble>(df).convert(newtype);

		case degminsec:
			return waypointsmaker<DegMinSecVecDouble>(df).convert(newtype);

		default:
			stop("convert_switch(const DataFrame) const my bad");
	}
}

/// __________________________________________________
/// Format "waypoints" DataFrame 
const bisvec<string> format_switch(const DataFrame df, CoordType ct_required)
{
#if DEBUG > 0
	fmt::print("@format_switch(const DataFrame, CoordType); current type: {}, required type: {}\n", get_coordtype(df), ct_required);
#endif
	using enum CoordType;
	switch (get_coordtype(df))
	{
		case decdeg:
			return waypointsmaker<DecDegVecDouble>(df).format(ct_required);

		case degmin:
			return waypointsmaker<DegMinVecDouble>(df).format(ct_required);

		case degminsec:
			return waypointsmaker<DegMinSecVecDouble>(df).format(ct_required);

		default:
			stop("format_switch(const DataFrame) const my bad");
	}
}

/// __________________________________________________
/// Validate "waypoints" DataFrame 
const bisconstvec <bool> validate_switch(const DataFrame df)
{
#if DEBUG > 0
	fmt::print("@validate_switch(const DataFrame); current type: {}\n", get_coordtype(df));
#endif
	using enum CoordType;
	switch (get_coordtype(df))
	{
		case decdeg:
			return waypointsmaker<DecDegVecDouble>(df).validate();

		case degmin:
			return waypointsmaker<DegMinVecDouble>(df).validate();

		case degminsec:
			return waypointsmaker<DegMinSecVecDouble>(df).validate();

		default:
			stop("validate_switch(const DataFrame) const my bad");
	}
}


/// __________________________________________________
/// __________________________________________________
/// Validation functions

/// __________________________________________________
/// Check "valid" attribute of NumericVector all true
bool check_valid(const NumericVector nv, bool newbie)
{
#if DEBUG > 0
	fmt::print("@check_valid(const NumericVector)\n");
#endif
	int validated = check_logical_attr(nv, "valid"s);
	if (!validated)
		return validate(nv, !newbie);
	return validated >> 1;
}

/// __________________________________________________
/// Check "lat_valid" and "lon_valid attributes of DataFrame are all true
bool check_valid(const DataFrame df, bool newbie)
{
#if DEBUG > 0
	fmt::print("@check_valid(const DataFrame)\n");
#endif
	int latvalidated = check_logical_attr(df, "validlat"s);
	int lonvalidated = check_logical_attr(df, "validlon"s);

	if (!(latvalidated & lonvalidated))
		return validate(df, !newbie);
	if (!(latvalidated >> 1))
		warning("Invalid latitude!");
	if (!(lonvalidated >> 1))
		warning("Invalid longitude!");
	return latvalidated >> 1 && lonvalidated >> 1;
}

/// __________________________________________________
/// Validate "coords" NumericVector or "waypoints" DataFrame
bool validate(const NumVec_or_DataFrame auto t, bool revalidate)
{
#if DEBUG > 0
	fmt::print("@validate(const NumVec_or_DataFrame auto, bool); revalidate: {}, t {}\n", revalidate, demangle(typeid(t)));
#endif
	using t_type = std::remove_const_t<decltype(t)>;
	bool iscoords {false};
	bool warn {false};
	auto valid { validate_switch(t) };
	if constexpr (isNumericVector_v<t_type>) {
#if DEBUG > 0
		fmt::print("@IIvalidate(const NumVec_or_DataFrame auto, bool);\n\t{}\t  &valid {}, typeid: {}, \n\t{}\n",
			padstr, address(valid), demangle(typeid(valid)), fmt::join(valid, ", "));
#endif
		iscoords = true;
		if (!std::all_of(valid.begin(), valid.end(), [](auto i){ return i; }))
			warn = true;
		static_cast<NumericVector>(t).attr("valid") = valid; 
	} else if constexpr (isDataFrame_v<t_type>) {
#if DEBUG > 0
		fmt::print("@IIIvalidate(const NumVec_or_DataFrame auto, bool); &valid {}, &valid[0] {}, &valid[1]{}, typeid: {}, \n\t{}, \n\t{}\n",
			address(valid), address(valid[0]), address(valid[1]), demangle(typeid(valid)), fmt::join(valid[0], ", "), fmt::join(valid[1], ", "));
#endif
		if (!std::all_of(valid[0].begin(), valid[0].end(), [](auto i){ return i; }) ||
			!std::all_of(valid[1].begin(), valid[1].end(), [](auto i){ return i; }))
			warn = true;
		static_cast<DataFrame>(t).attr("validlat") = valid[0];
		static_cast<DataFrame>(t).attr("validlon") = valid[1];
	} else
		stop("validate(const NumVec_or_DataFrame auto, bool revalidate) my bad!");
	if (warn)
		warning("%salidation detected invalid %s!", revalidate ? "Rev" : "V", iscoords ? "coords" : "waypoints");
	else if (revalidate)
		warning("%s revalidated!", iscoords ? "Coords" : "Waypoints");
	return check_valid(t);
}

/// __________________________________________________
/// Check df has valid "llcols" attribute
bool valid_ll(const DataFrame df)
{
#if DEBUG > 0
	fmt::print("@{}\n", "valid_ll(const DataFrame)");
#endif
	bool valid = false;
	vector llcols { get_vec_attr<DataFrame, int>(df, "llcols"s) };
	if (2 == llcols.size()) {
		transform(llcols.begin(), llcols.end(), llcols.begin(), [](auto x){ return --x; });
		if (is_item_in_df(df, llcols[0]) && is_item_in_df(df, llcols[1]) && llcols[0] != llcols[1])
			if (is<NumericVector>(df[llcols[0]]) && is<NumericVector>(df[llcols[1]]))
				valid = true;
	}
	return valid;
}


/// __________________________________________________
/// __________________________________________________
/// Exported functions

#if DEBUG < 3
/*
#endif

/// __________________________________________________
/// Dummy Function for Testing Only 	¡¡¡ ——— Temporary to Be Archived ——— !!!
//' @rdname cords
// [[Rcpp::export(name = "ctors")]]
NumericVector constructors(NumericVector object)
{ 
	fmt::print("{}@ctors(NumVec); fmt {}\n", exportstr, get_fmt_attribute(object));

	constexpr int x{ 7 };
	using addr = decltype(address(object));
	auto prevaddr = addr{ 0x00 }; 
	fmt::print("\t—initialise prevaddr {}\n", prevaddr);
	bool moved{false};

	constexpr auto mov { "moved"sv };
	constexpr auto cpy { "copied"sv };

// NumericVector "object" arg 
	fmt::print("{}I@ctors(NumVec); &object {}, object[x] {}, &object[x] {}\n", exportstr, address(object), object[x], address(object[x]));
	prevaddr = address(object[x]);
	
// Assign to vector<double> from NumericVector —— copies content of "object" to nv at a new address;
	auto nv { as<vector<double>>(object) };
	moved = address(nv[x]) == prevaddr; 
	fmt::print("{}II@ctors(NumVec); &nv {}, nv[x] {}, ({}) &nv[x] {}\n", exportstr, address(nv), nv[x], moved ? mov : cpy, address(nv[x]));
	prevaddr = address(nv[x]);

// Copy construct one vector<double> from another —— copies content of nv to nv2 at a new address
	auto nv2 { nv };
	moved = address(nv2[x]) == prevaddr; 
	fmt::print("{}III@ctors(NumVec); &nv2 {}, nv2[x] {}, ({}) &nv2[x] {}\n", exportstr, address(nv2), nv2[x], moved ? mov : cpy, address(nv2[x]));

// Move construct one vector<double> from another —— moves content of nv to nv3, content maintains its address
	auto nv3 { std::move(nv) };
	moved = address(nv3[x]) == prevaddr; 
	fmt::print("{}IV@ctors(NumVec); &nv3 {}, nv3[x] {}, ({}) &nv3[x] {}\n", exportstr, address(nv3), nv3[x], moved ? mov : cpy, address(nv3[x]));

// Check original "movee" nv —— content now undefined
	moved = address(nv[x]) != prevaddr; 
	fmt::print("{}V@ctors(NumVec); (post move) &nv {}, nv[x] {}, ({}) &nv[x] {}\n", exportstr, address(nv), "undefined!", moved ? mov : "same", address(nv[x]));

// Copy construct DecDegVecDouble from vector<double> —— copies content of nv3 to dv at a new address
	prevaddr = address(nv3[x]);
	DecDegVecDouble dv { nv3 }; 
	moved = address(dv[x]) == prevaddr; 
	fmt::print("{}VI@ctors(NumVec); &dv {}, dv[x] {}, ({}) &dv[x] {}\n", exportstr, address(dv), dv[x], moved ? mov : cpy, address(dv[x]));

// Check original "copyee" nv3 —— content unchanged
	moved = address(nv3[x]) != prevaddr; 
	fmt::print("{}VII@ctors(NumVec); (post copy) &nv3 {}, nv3[x] {}, ({}) &nv3[x] {}\n", exportstr, address(nv3), nv3[x], moved ? mov : "same", address(nv3[x]));

// Move construct DecDegVecDouble from vector<double> —— moves content of nv3 to dv2, content maintains its address
	DecDegVecDouble dv2 { std::move(nv3) }; 
	moved = address(dv2[x]) == prevaddr; 
	fmt::print("{}VIII@ctors(NumVec); &dv2 {}, dv2[x] {}, ({}) &dv2[x] {}\n", exportstr, address(dv2), dv2[x], moved ? mov : cpy, address(dv2[x]));

// Check original "movee" nv3 —— content now undefined 
	fmt::print("{}IX@ctors(NumVec); (post move) &nv3 {}, nv3[x] {}, ({}) &nv3[x] {}\n", exportstr, address(nv3), "undefined", moved ? mov : "same", address(nv3[x]));

// Copy construct one DecDegVecDouble from another —— copies content of dv2 to dv3 at a new address
	prevaddr = address(dv2[x]);
	DecDegVecDouble dv3 { dv2 }; 
	moved = address(dv3[x]) == prevaddr; 
	fmt::print("{}X@ctors(NumVec); &dv3 {}, dv3[x] {}, ({}) &dv3[x] {}\n", exportstr, address(dv3), dv3[x], moved ? mov : cpy, address(dv3[x]));

// Check original "copyee" dv2 —— content unchanged
	moved = address(dv2[x]) != prevaddr; 
	fmt::print("{}XI@ctors(NumVec); (post copy) &dv2 {}, dv2[x] {}, ({}) &dv2[x] {}\n", exportstr, address(dv2), dv2[x], moved ? mov : "same", address(dv2[x]));

// Move construct one DecDegVecDouble from another —— moves content of dv2 to dv4, content maintains its address (was in dv2, nv3, nv)
	DecDegVecDouble dv4 { std::move(dv2) }; 
	moved = address(dv4[x]) == prevaddr; 
	fmt::print("{}XII@ctors(NumVec); &dv4 {}, dv4[x] {}, ({}) &dv4[x] {}\n", exportstr, address(dv4), dv4[x], moved ? mov : cpy, address(dv4[x]));

// Check original "movee" dv2 —— content now undefined
	moved = address(dv2[x]) != prevaddr; 
	fmt::print("{}XIII@ctors(NumVec); (post move) &dv2 {}, dv2[x] {}, ({}) &dv2[x] {}\n", exportstr, address(dv2), "undefined", moved ? mov : "same", address(dv2[x]));
	
// Copy assign to one DecDegVecDouble from another —— using signature dv3[x] = 123.456789; copies content of dv3 to dv4 at existing address
	dv3[x] = 123.456789;
	prevaddr = address(dv3[x]); 
	fmt::print("{}XIV@ctors(NumVec); (post change dv3[x]) &dv3 {}, dv3[x] {}, &dv3[x] {}\n", exportstr, address(dv3), dv3[x], address(dv3[x]));
	dv4 = dv3; 
	fmt::print("{}XV@ctors(NumVec); (post copy assign) &dv4 {}, dv4[x] {}, &dv4[x] {}\n", exportstr, address(dv4), dv4[x], address(dv4[x]));

// Check original "copyee" dv3 —— content unchanged
	moved = address(dv3[x]) != prevaddr; 
	fmt::print("{}XVI@ctors(NumVec); (post copy assign) &dv3 {}, dv3[x] {}, ({}) &dv3[x] {}\n", exportstr, address(dv3), dv3[x], moved ? mov : "same", address(dv3[x]));

// Move assign to one DecDegVecDouble from another —— moves content of dv4 to dv3, content maintains its address (was in dv2, nv3, nv)
	prevaddr = address(dv4[x]);
	dv3 = (std::move(dv4)); 
	moved = address(dv3[x]) == prevaddr; 
	fmt::print("{}XVII@ctors(NumVec); (post move assign) &dv3 {}, dv3[x] {}, ({}) &dv3[x] {}\n", exportstr, address(dv3), dv3[x], moved ? mov : cpy, address(dv3[x]));

// Check original "movee" dv4 —— content now undefined
	moved = address(dv4[x]) != prevaddr; 
	fmt::print("{}XVIII@ctors(NumVec); (post move assign) &dv4 {}, dv4[x] {}, ({}) &dv4[x] {}\n", exportstr, address(dv4), "undefined", moved ? mov : "same", address(dv4[x]));

	return object;
}

/// __________________________________________________
/// Dummy Function for Testing Only 	¡¡¡ ——— Temporary to Be Archived ——— !!!
//' @rdname cords
// [[Rcpp::export(name = "cordelia")]]
NumericVector movit(NumericVector object)
{ 
	fmt::print("{}@movit(NumericVector); fmt {}\n", exportstr, get_fmt_attribute(object));
	using enum CoordType;

// NumericVector "object" arg 
	fmt::print("{}I@movit(NumVec); &object {}, object[0] {}, &object[0] {}\n", exportstr, address(object), object[0], address(object[0]));

// Assign to vector<double> nv from NumericVector —— copies content of "object" to a new address
	auto nv { as<vector<double>>(object) }; 
	fmt::print("{}II@movit(NumVec); &nv {}, nv[0] {}, &nv[0] {}\n", exportstr, address(nv), nv[0], address(nv[0]));

// Move construct DecDegVecDouble from vector<double> —— moves content from nv to dv, maintains address
	DecDegVecDouble dv { std::move(nv) }; 
	fmt::print("{}III@movit(NumVec); &dv {}, dv[0] {}, &dv[0] {}\n", exportstr, address(dv), dv[0], address(dv[0])); 
	fmt::print("{}IV@movit(NumVec); &nv {}, nv[0] {}, &nv[0] {}\n", exportstr, address(nv), "undefined", address(nv[0]));

// Move construct Coordlet<DecDegVecDouble> from DecDegVecDouble —— moves content  from dv to clt, maintains address 
	fmt::print("{}V@movit(NumericVector); Constructing Coordlet<DecDegVecDouble> by moving dv\n", exportstr);
	auto clt = Coordlet<DecDegVecDouble>{ std::move(dv), vector<bool>{} };
	clt.report(); 
	fmt::print("{}VI@movit(NumVec); &dv {}, dv[0] {}, &dv[0] {}\n", exportstr, address(dv), "undefined", address(dv[0]));

// Re-assign to vector<double> nv from NumericVector —— copies "object" content to new address
	nv = as<vector<double>>(object); 
	fmt::print("{}VII@movit(NumVec); &nv {}, nv[0] {}, &nv[0] {}\n", exportstr, address(nv), nv[0], address(nv[0]));

// Construct Coordlet<DecDegVecDouble> indirectly from vector<double>; using *, new and delete to observe destruction 
	fmt::print("{}VIII@movit(NumericVector); Making ptr1\n", exportstr);
	const auto* ptr1 = new Coordlet<DecDegVecDouble>{ std::move(nv), vector<bool>{} };	// std::move() needed here—no copy elision.
	ptr1->report(); 
	fmt::print("{}IX@movit(NumVec); &nv {}, nv[0] {}, &nv[0] {}\n", exportstr, address(nv), "undefined", address(nv[0]));
	delete ptr1;

// Re-assign to vector<double> vn from NumericVector —— copies "object" content to new address (may re-use &dv[0] from deleted Coordlet)
	nv = as<vector<double>>(object); 
	fmt::print("{}X@movit(NumVec); &nv {}, nv[0] {}, &nv[0] {}\n", exportstr, address(nv), nv[0], address(nv[0]));

// Construct Coordlet<DecDegVecDouble> directly from DecDegVecDouble; using *, new and delete to observe destruction 
	fmt::print("{}XI@movit(NumericVector); Making ptr2\n", exportstr);
	const auto* ptr2 = new Coordlet<DecDegVecDouble>{ DecDegVecDouble{ std::move(nv) }, vector<bool>{} };
	ptr2->report(); 
	fmt::print("{}XI@movit(NumVec); &nv {}, nv[0] {}, &nv[0] {}\n", exportstr, address(nv), "undefined", address(nv[0]));
	delete ptr2; 
	fmt::print("{}XI@movit(NumericVector); ptr2 deleted, clt to be deleted on exit\n", exportstr);
	return object;
}
#if DEBUG < 2
*/
#endif

/// __________________________________________________
/// Create coords - S3 method as_coords.default()
//' @rdname coords 
// [[Rcpp::export(name = "as_coords.default")]]
NumericVector as_coords(NumericVector object, int fmt = 1)
{
#if DEBUG > 0
	fmt::print("{}@as_coords(NumericVector, int); fmt={}\n", exportstr, fmt);
#endif
	object.attr("fmt") = fmt;
	if (!check_valid(object, true))
		warning("[Use review() to show invalid elements]");
	object.attr("class") = "coords";
	return object;
}

/// __________________________________________________
/// Convert coords - S3 method convert.coords()
//' @rdname convert
// [[Rcpp::export(name = "convert.coords")]]
NumericVector convertcoords(const NumericVector x, int fmt)
{
	checkinherits(x, "coords"s);
	CoordType ct_current = get_coordtype(x);
	CoordType newtype = get_coordtype(fmt);
#if DEBUG > 0
	fmt::print("{}@convertcoords(NumericVector, int); from {} to {}\n", exportstr, ct_current, newtype);
	fmt::print("{}@Iconvertcoords(NumericVector, int);\n\t{}\t\t  &x {}, &x[0] {}, x[0] {}, typeid: {}\n",
		exportstr, padstr, address(x), address(x[0]), x[0], demangle(typeid(x)));
#endif
	if (!check_valid(x))
		stop("Invalid coords! Conversion aborted.\n [Use review() to show invalid elements]");
	if (newtype != ct_current) {
		auto vd_out { convert_switch(x, newtype) };
#if DEBUG > 0
		fmt::print("{}@IIconvertcoords(NumericVector, int);\n\t{}\t &vd_out {}, &vd_out[0] {}, vd_out[0] {}, typeid: {}\n",
			exportstr, padstr, address(vd_out), address(vd_out[0]), vd_out[0], demangle(typeid(vd_out)));
#endif
		NumericVector nv_out { wrap(vd_out) };									// Copies output string
#if DEBUG > 0
		fmt::print("{}@IIconvertcoords(NumericVector, int);\n\t{}\t &nv_out {}, &nv_out[0] {}, nv_out[0] {}, typeid: {}\n",
			exportstr, padstr, address(nv_out), address(nv_out[0]), nv_out[0], demangle(typeid(nv_out)));
#endif
		nv_out.attr("class") = "coords";
		nv_out.attr("fmt") = fmt;
		nv_out.attr("valid") = x.attr("valid");
		nv_out.attr("latlon") = x.attr("latlon");
		nv_out.names() = x.names();
		return nv_out;
	} else { 
		Rcout << "\t—— fmt out == fmt in: returning original x!——\n\n";
		return x;														//	¡¡¡—— D'you think that's wise, sir? ——!!!
	}
}

/// __________________________________________________
/// Set latlon attribute on "coords" NumericVector and revalidate
//' @rdname coords
// [[Rcpp::export(name = "`latlon<-`")]]
NumericVector latlon(NumericVector cd, LogicalVector value)
{
#if DEBUG > 0
	fmt::print("{}@latlon(NumericVector, LogicalVector)\n", exportstr);
#endif
	checkinherits(cd, "coords"s);
	if (value.size() != cd.size() && value.size() != 1)
		stop("value must be either length 1 or length(cd)");
	else
		cd.attr("latlon") = value;
	if (!validate(cd))
		warning("[Use review() to show invalid elements]");	
	return cd;
}

/// __________________________________________________
/// Format coords - S3 method format.coords()
//' @rdname format
// [[Rcpp::export(name = "format.coords")]]
CharacterVector formatcoords(const NumericVector x, bool usenames = true, bool validate = true, int fmt = 0)
{
#if DEBUG > 0
	fmt::print("{}@formatcoords(const NumericVector, bool, bool, int); usenames: {}, validate: {}, fmt: {}\n", exportstr, usenames, validate, fmt);
	fmt::print("{}@Iformatcoords(const NumericVector, bool, bool, int);\n\t{}\t\t  &x {}, &x[0] {}, x[0] {}, typeid: {}\n",
		exportstr, padstr, address(x), address(x[0]), x[0], demangle(typeid(x)));
#endif
	checkinherits(x, "coords"s);
	if(!x.size())
		stop("x has 0 length!");
	if (validate)
		if (!check_valid(x))
			warning("Formatting invalid coords!\n [Use review() to show invalid elements]");
	CoordType ct_current { get_coordtype(x) };
	CoordType ct_required { fmt ? get_coordtype(fmt) : ct_current };
	auto sv_out { format_switch(x, ct_required) };
#if DEBUG > 0
	fmt::print("{}@IIformatcoords(const NumericVector, bool, bool, int);\n\t{}\t\t&sv_out {}, &sv_out[0] {}, sv_out[0] {}, typeid: {}\n",
		exportstr, padstr, address(sv_out), address(sv_out[0]), sv_out[0], demangle(typeid(sv_out).name()));
#endif
	vector names{ get_vec_attr<NumericVector, string>(x, "names"s) };
	if (names.size() && usenames) {
		stdlenstr(names);
		concat_vecstr_elmnts(names, sv_out);
	}
	return wrap(sv_out);
}

/// __________________________________________________
/// Validate coords - S3 method validate.coords()
//' @rdname validate
// [[Rcpp::export(name = "validate.coords")]]
NumericVector validatecoords(const NumericVector x, const bool force = true)
{
#if DEBUG > 0
	fmt::print("{}@validatecoords(const NumericVector, const bool); force: {}\n", exportstr, force);
	fmt::print("{}@Ivalidatecoords(const NumericVector, const bool);\n\t{}\t\t  &x {}, &x[0] {}, x[0] {}, typeid: {}\n",
		exportstr, padstr, address(x), address(x[0]), x[0], demangle(typeid(x)));
#endif
	checkinherits(x, "coords"s);
	bool warn { false };
	if (force)	{			
		if (!validate(x))
			warn = true;
	} else if (!check_valid(x))
		warn = true;
	if (warn)
		warning("[Use review() to show invalid elements]");
	return x;
}

/// __________________________________________________
/// Create waypoints - S3 method as_waypoints.default()
//' @rdname waypoints
// [[Rcpp::export(name = "as_waypoints.default")]]
DataFrame as_waypoints(DataFrame object, int fmt = 1)
{
#if DEBUG > 0
#endif
#if DEBUG > 0
	fmt::print("{}@as_waypoints(DataFrame, int); fmt={}\n", exportstr, fmt);
#endif
	object.attr("fmt") = fmt;
	int namescol = 0;
	if (!object.hasAttribute("namescol")) {
		namescol = name_pos_in_df(object, "name"s);
		if (++namescol)
			object.attr("namescol") = namescol;
	}
	if (!object.hasAttribute("llcols")) {
		const vector llcols{ namescol + 1, namescol + 2 };
		object.attr("llcols") = llcols;
	}
	if(!valid_ll(object))
		stop("Invalid llcols attribute!");
	if (!check_valid(object, true))
		warning("[Use review() to show invalid elements]");
	object.attr("class") = CharacterVector{"waypoints", "data.frame"};
	return object;
}


/// __________________________________________________
/// Convert waypoints type - S3 method convert.waypoints()
//' @rdname convert
// [[Rcpp::export(name = "convert.waypoints")]]
DataFrame convertwaypoints(DataFrame x, int fmt)
{
	checkinherits(x, "waypoints"s);
	CoordType ct_current = get_coordtype(x);
	CoordType newtype = get_coordtype(fmt);
#if DEBUG > 0
	fmt::print("{}@convertwaypoints(DataFrame, int); from {} to {}\n", exportstr, ct_current, newtype);
#endif
	if (!check_valid(x))
		stop("Invalid waypoints! Conversion aborted.\n [Use review() to show invalid elements]");
	if(!valid_ll(x))
		stop("Invalid llcols attribute!");
	if (newtype != ct_current) {
		auto vds_out { convert_switch(x, newtype) };
#if DEBUG > 0
		fmt::print("{}@Iconvertwaypoints(DataFrame, int);  vds_out[0][0] {}, &vds_out[0] {}, &vds_out[0][0] {}, typeid: {}\n",
			exportstr, vds_out[0][0], address(vds_out[0]), address(vds_out[0][0]), demangle(typeid(vds_out[0])));
		fmt::print("{}@IIconvertwaypoints(DataFrame, int);  vds_out[1][0] {}, &vds_out[1] {}, &vds_out[1][0] {}, typeid: {}\n",
			exportstr, vds_out[1][0], address(vds_out[1]), address(vds_out[1][0]), demangle(typeid(vds_out[1])));
#endif
		auto llcols { get_vec_attr<decltype(x), int>(x, "llcols") };
		auto namescol { get_vec_attr<decltype(x), int>(x, "namescol") };
		auto names { get_vec_attr<decltype(x), string>(x, "names") };
		auto row_names { get_vec_attr<decltype(x), int>(x, "row.names") };
		auto validlat { get_vec_attr<decltype(x), bool>(x, "validlat") };
		auto validlon { get_vec_attr<decltype(x), bool>(x, "validlon") };

		auto llcol_it { x.erase(llcols[0] - 1) };
		x.insert(llcol_it, vds_out[0]);

		llcol_it = x.erase(llcols[1] - 1);
		x.insert(llcol_it, vds_out[1]);

		x.attr("names") = names;
		x.attr("class") = vector{"waypoints", "data.frame"};
		x.attr("row.names") = row_names;
		x.attr("fmt") = fmt;
		x.attr("namescol") = namescol;
		x.attr("llcols") = llcols;
		x.attr("validlat") = validlat;
		x.attr("validlon") = validlon;
	} else
		Rcout << "\t—— fmt out == fmt in! ——\n\n";
	return x;
}


/// __________________________________________________
/// Format waypoints - S3 method format.waypoints()
//' @rdname format
// [[Rcpp::export(name = "format.waypoints")]]
CharacterVector formatwaypoints(DataFrame x, bool usenames = true, bool validate = true, int fmt = 0)
{
#if DEBUG > 0
	fmt::print("{}@formatwaypoints(DataFrame, bool, bool, int); usenames: {}, validate: {}, fmt: {}\n", exportstr, usenames, validate, fmt);
#endif
	checkinherits(x, "waypoints"s);
	if(!x.nrows())
		stop("x has 0 rows!");
	if(!valid_ll(x))
		stop("Invalid llcols attribute!");
	if (validate)
		if (!check_valid(x))
			warning("Formatting invalid waypoints!");
	auto svs_out { format_switch(x, fmt ? get_coordtype(fmt) : get_coordtype(x)) };
	auto& sv_out { svs_out[0] };	// AKA "sv_lat"
	auto& sv_lon { svs_out[1] };
    transform(sv_out.begin(), sv_out.end(), sv_lon.begin(), sv_out.begin(), [](auto& latstr, auto& lonstr){ return latstr + "  " + lonstr; });
#if DEBUG > 0
	fmt::print("{}@Iformatwaypoints(DataFrame, bool, bool, int); sv_out[0] {}, &sv_out {}, &sv_out[0] {}, typeid: {}\n",
		exportstr, sv_out[0], address(sv_out), address(sv_out[0]), demangle(typeid(sv_out)));
#endif
	if (usenames) {
		RObject names = getnames(x);
		if (!prefixwithnames(sv_out, names))
			stop("Invalid \"namescol\" attribute!");
	}
#if DEBUG > 0
	fmt::print("{}@IIformatwaypoints(DataFrame, bool, bool, int); sv_out[0] {}, &sv_out {}, &sv_out[0] {}, typeid: {}\n",
		exportstr, sv_out[0], address(sv_out), address(sv_out[0]), demangle(typeid(sv_out)));
#endif
	return wrap(sv_out);
}

/// __________________________________________________
/// Validate waypoints - S3 method validate.waypoints()
//' @rdname validate
// [[Rcpp::export(name = "validate.waypoints")]]
DataFrame validatewaypoints(DataFrame x, bool force = true)
{
#if DEBUG > 0
	fmt::print("{}@validatewaypoints(DataFrame, bool); force: {}\n", exportstr, force);
#endif
	checkinherits(x, "waypoints"s);
	if(!valid_ll(x))
		stop("Invalid llcols attribute!");
	bool warn { false };
	if (force)	{			
		if (!validate(x))
			warn = true;
	} else if (!check_valid(x))
		warn = true;
	if (warn)
		warning("[Use review() to show invalid elements]");
	return x;
}

/// __________________________________________________
/// Latitude and longitude headers for S3 print.waypoint()
//' @rdname format
// [[Rcpp::export]]
CharacterVector ll_headers(int width, int fmt)
{
#if DEBUG > 0
	fmt::print("{}@ll_headers(int, int); width={}, fmt={}\n", exportstr, width, fmt);
#endif
	--fmt;  //	  to C++ array numbering
	constexpr int spacing[][3] { {15,  17,  18}, {11, 13, 14} };
	return wrap(vector {
		fmt::format("{:>{}}{:>{}}", "Latitude", width - spacing[0][fmt], "Longitude", spacing[0][fmt] - 1), // --fmt —> C++ array numbering
		fmt::format("{:>{}}", string(spacing[1][fmt], '_') + string(2, ' ') + string(spacing[1][fmt] + 1, '_'), width),
	});
}
/*
/// __________________________________________________
/// Clone coords object from waypoints vector
//' @rdname coords
// [[Rcpp::export(name = "as_coords.waypoints")]]
NumericVector as_coordswaypoints(DataFrame object, bool which)
{
#if DEBUG > 0
	fmt::print("{}@as_coord(DataFrame); which: {}\n", exportstr, which ? "lat" : "lon");
#endif
	checkinherits(object, "waypoints"s);
	NumericVector nv = object[get_vec_attr<DataFrame, int>(object, "llcols"s)[which ? 0 : 1] - 1];
	nv = clone(nv);
	nv.attr("class") = "coords";
	nv.attr("fmt") = object.attr("fmt");
	nv.attr("valid") = object.attr(which ? "validlat" : "validlon");
	nv.attr("latlon") = which;
	nv.attr("names") = getnames(object);
	return nv;
}
*/
/// __________________________________________________
/// __________________________________________________

#if DEBUG > 1
.
#endif
