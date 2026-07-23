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
#include <iostream>     // std::cout
#include <sstream>      // std::ostringstream

using namespace Rcpp;

using std::array;
using std::vector;
using std::string;
using namespace std::literals;
using std::string_view;
using namespace std::string_view_literals;
using std::transform;
using std::ostringstream;
using std::fixed;
using std::left;
using std::setw;
using std::setfill;
using std::setprecision;

#include "CoordBase.h"

#define FMT_HEADER_ONLY
#include "fmt/format.h"		// …fmt/*.h copied to ~/Documents/R/Packages/Waypoint/src/fmt.

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
inline double round2(double x, int n)
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
template<typename U> 
inline vector<U> get_vec_attr(const NumVec_or_DataFrame auto& t, const string attrname)
{
	return t.hasAttribute(attrname) ? as<vector<U>>(t.attr(attrname)) : vector<U>{};
}

/// __________________________________________________
/// Return "fmt" attribute as int
inline int get_fmt_attribute(const NumVec_or_DataFrame auto& t)
{
	return as<int>(t.attr("fmt"));
}

/// __________________________________________________
/// Check whether a NumericVector or DataFrame has a specified logical vector attribute and whether all true
int check_logical_attr(NumVec_or_DataFrame auto t, const string attrname)
{
	const vector vec_attr{ get_vec_attr<bool>(t, attrname) };
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
	if (!t.inherits(classname.c_str())) stop("Argument must be a \"%s\" object", classname.c_str());
}

/// __________________________________________________
/// Is item number present in data.frame? (Using C++ numbering)
inline bool is_item_in_df(const DataFrame df, int item_no)
{
	if (NA_INTEGER == item_no)
		return false;
	else
		return !(item_no < 0) && item_no < df.size();
}

/// __________________________________________________
/// Standarise width of strings in vector to that of the longest
inline void stdlenstr(vector<string>& sv)
{
	auto maxwdth = max_element(sv.begin(), sv.end(), [](const string& a, const string& b){ return a.size() < b.size(); })->size();
	transform(sv.begin(), sv.end(), sv.begin(), [maxwdth](const string& s) { return fmt::format("{:<{}}", s, maxwdth); });
}

/// __________________________________________________
/// Concatenate corresponding elements of two vector<string>, with separator; result in second vector<string>
inline void concat_vecstr_elmnts(const vector<string>& sv_a, vector<string>& sv_b, const string sep)
{
	transform(sv_a.begin(), sv_a.end(), sv_b.begin(), sv_b.begin(), [&sep](const string& str_a, const string& str_b) {
		return str_a + sep + str_b; }); 
}

/// __________________________________________________
/// Concatenate corresponding elements of vector<int> and vector<string>, with separator; result in vector<string>
inline void concat_vecstr_elmnts(const vector<int>& iv_a, vector<string>& sv_b, const string sep)
{
	transform(iv_a.begin(), iv_a.end(), sv_b.begin(), sv_b.begin(), [&sep](const int i, const string& str_b) {
		return (std::to_string(i)) + sep + str_b; }); 
}

/// __________________________________________________
/// Prefix vector<string> elements with elements of RObject
inline bool prefixwithnames(vector<string>& sv, RObject& namesobj)
{
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
	transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return tolower(c); });
	return s;
}

/// __________________________________________________
/// Find position of name within data.frame names
int name_pos_in_df(const DataFrame df, const string name)
{
	vector names{ get_vec_attr<string>(df, "names"s) };
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
	vector namescolvec{ get_vec_attr<int>(df, "namescol"s) };
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
	return get_coordtype(get_fmt_attribute(t));
}

/// __________________________________________________
/// Convert CoordType enum to int; + 1 for R
inline int coordtype_to_int(CoordType ct)
{
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
template<DVecType T, typename S>
Coords<T, S>::Coords(NumericVector nv) :
	dv { std::move(as<vector<double>>(nv)) },
	latlon { get_vec_attr<bool>(nv, "latlon"s) }
{
	static_assert(sufijo<S> && std::derived_from<S, Coords>);
}

/// __________________________________________________
/// Format dv as a vectype object —— private
template<DVecType T, typename S> template<vectype U, functador V>
inline U Coords<T, S>::conform0() const
{
	U uv_out(dv.size());
	transform(dv.begin(), dv.end(), uv_out.begin(), V());	
	if constexpr (SVecType<U>)
		add_suffix(uv_out);
	return uv_out;
}

/// __________________________________________________
/// conform call entry point —— public
template<DVecType T, typename S> template<typename U, template <typename V> typename F>
vector<U> Coords<T, S>::conform(CoordType required) const
{
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
/// Validation call entry point —— public
template<DVecType T, typename S>
const vector<bool> Coords<T, S>::validate() const
{
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

	return valid;
}


/// __________________________________________________
/// Suffix call entry point
template<DVecType T, typename S>
void Coords<T, S>::add_suffix(vectype auto& uv_out) const
{
	static_cast<const S *>(this)->suffix(uv_out);
}


/// __________________________________________________
/// __________________________________________________
/// SufijoCoords class —— add suffix
template<DVecType T>
void SufijoCoords<T>::suffix(vectype auto& uv_out) const
{
	using uv_out_type = std::remove_cvref_t<decltype(uv_out)>;
	vector<bool>::const_iterator ll_it { latlon.begin() };
	const auto ll_size { latlon.size() };

	if constexpr (isDecDegVecString_v<uv_out_type>) {
		const auto lambda1 = [&ll_it](auto& outstr, auto n){ return outstr + (*ll_it++ ? " lat" : " lon"); };
		const auto lambda2 = [&ll_it](auto& outstr, auto n){ return outstr + (*ll_it ? " lat" : " lon"); };
	
		if (ll_size > 1)
			transform(uv_out.begin(), uv_out.end(), dv.begin(), uv_out.begin(), lambda1);
		else
			if (ll_size == 1)   // uniform coords
				transform(uv_out.begin(), uv_out.end(), dv.begin(), uv_out.begin(), lambda2);
	} else if constexpr (isDegMinVecString_v<uv_out_type> || isDegMinSecVecString_v<uv_out_type>) {
		const auto lambda1 = [&ll_it](auto& outstr, auto n){ return outstr + cardpoint(n < 0, *ll_it++); };
		const auto lambda2 = [&ll_it](auto& outstr, auto n){ return outstr + cardpoint(n < 0, *ll_it); };
		const auto lambda3 = [](auto& outstr, auto n){ return outstr + cardi_b(n < 0); };
	
		if (ll_size > 1)
			transform(uv_out.begin(), uv_out.end(), dv.begin(), uv_out.begin(), lambda1);
		else
			if (ll_size == 1)   // uniform coords
				transform(uv_out.begin(), uv_out.end(), dv.begin(), uv_out.begin(), lambda2);
			else				// no latlon info
				transform(uv_out.begin(), uv_out.end(), dv.begin(), uv_out.begin(), lambda3);
	}
}


/// __________________________________________________
/// __________________________________________________
/// SufijoWaypoints class —— add suffix
template<DVecType T>
void SufijoWaypoints<T>::suffix(vectype auto& uv_out) const
{
	using uv_out_type = std::remove_cvref_t<decltype(uv_out)>;
	if constexpr (!isDecDegVecString_v<uv_out_type>)
		transform(uv_out.begin(), uv_out.end(), dv.begin(), uv_out.begin(), [this](auto& outstr, auto n){
		   return outstr + cardpoint(n < 0, latlon[0]); }
		);
}


/// __________________________________________________
/// Convert "coords" NumericVector
vector<double> convert_switch(const NumericVector nv, CoordType newtype)
{
	using enum CoordType;
	switch (get_coordtype(nv))
	{
		case decdeg:
			return SufijoCoords<DecDegVecDouble>(nv).template conform<double, ConvertidorDecDegVec>(newtype);

		case degmin:
			return SufijoCoords<DegMinVecDouble>(nv).template conform<double, ConvertidorDegMinVec>(newtype);

		case degminsec:
			return SufijoCoords<DegMinSecVecDouble>(nv).template conform<double, ConvertidorDegMinSecVec>(newtype);

		default:
			stop("convert_switch<sufijo>(const NumericVector, CoordType) const my bad");
	}
}

/// __________________________________________________
/// Format "coords" NumericVector with coords suffixes
vector<string> format_switch_c(const NumericVector nv, CoordType ct_required)
{
	using enum CoordType;
	switch (get_coordtype(nv))
	{
		case decdeg:
			return SufijoCoords<DecDegVecDouble>(nv).template conform<string, FormateadorDecDegVec>(ct_required);

		case degmin:
			return SufijoCoords<DegMinVecDouble>(nv).template conform<string, FormateadorDegMinVec>(ct_required);

		case degminsec:
			return SufijoCoords<DegMinSecVecDouble>(nv).template conform<string, FormateadorDegMinSecVec>(ct_required);

		default:
			stop("format_switch_c(const NumericVector, CoordType) const my bad");
	}
}

/// __________________________________________________
/// Format "coords" NumericVector with waypoints suffixes
vector<string> format_switch_w(const NumericVector nv, CoordType ct_required)
{
	using enum CoordType;
	switch (get_coordtype(nv))
	{
		case decdeg:
			return SufijoWaypoints<DecDegVecDouble>(nv).template conform<string, FormateadorDecDegVec>(ct_required);

		case degmin:
			return SufijoWaypoints<DegMinVecDouble>(nv).template conform<string, FormateadorDegMinVec>(ct_required);

		case degminsec:
			return SufijoWaypoints<DegMinSecVecDouble>(nv).template conform<string, FormateadorDegMinSecVec>(ct_required);

		default:
			stop("format_switch_w(const NumericVector, CoordType) const my bad");
	}
}


/// __________________________________________________
/// Validate "coords" NumericVector 
const vector<bool> validate_switch(const NumericVector nv)
{
	using enum CoordType;
	switch (get_coordtype(nv))
	{
		case decdeg:
			return SufijoCoords<DecDegVecDouble>(nv).validate();

		case degmin:
			return SufijoCoords<DegMinVecDouble>(nv).validate();

		case degminsec:
			return SufijoCoords<DegMinSecVecDouble>(nv).validate();

		default:
			stop("validate_switch(const NumericVector) const my bad");
	}
}


/// __________________________________________________
/// __________________________________________________
/// Waypoints class

/// __________________________________________________
/// Constructor
Waypoints::Waypoints(const DataFrame& df) :
	nv_lat( df[get_vec_attr<int>(df, "llcols")[0] - 1] ), 
	nv_lon( df[get_vec_attr<int>(df, "llcols")[1] - 1] )
{
	nv_lat.attr("fmt") = get_vec_attr<int>(df, "fmt");
	nv_lon.attr("fmt") = get_vec_attr<int>(df, "fmt");
	nv_lat.attr("latlon") = true;
	nv_lon.attr("latlon") = false;
}

/// __________________________________________________
/// Destructor
Waypoints::~Waypoints()
{
	nv_lat.attr("latlon") = R_NilValue;
	nv_lon.attr("latlon") = R_NilValue;
	nv_lat.attr("fmt") = R_NilValue;
	nv_lon.attr("fmt") = R_NilValue;
}

/// __________________________________________________
/// Convert nv_lat, nv_lon
vector<double> Waypoints::convert(CoordType newtype, bool latlon) const
{
	return convert_switch(latlon ? nv_lat : nv_lon, newtype);
}

/// __________________________________________________
/// Format nv_lat, nv_lon
vector<string> Waypoints::format(CoordType required_type, bool latlon) const
{
	using enum CoordType;
	auto sv_out { format_switch_w(latlon ? nv_lat : nv_lon, required_type) };
	return sv_out;
}

/// __________________________________________________
/// Convert nv_lat, nv_lon
const vector<bool> Waypoints::validate(bool latlon) const
{
	return validate_switch(latlon ? nv_lat : nv_lon);
}


/// __________________________________________________
/// Validate "waypoints" DataFrame 
inline const bisconstvec <bool> validate_switch(const DataFrame df)
{
	Waypoints wp{ df };
	return { wp.validate(true), wp.validate(false) };
}


/// __________________________________________________
/// __________________________________________________
/// Validation functions

/// __________________________________________________
/// Check "valid" attribute of NumericVector all true
bool check_valid(const NumericVector nv, bool newbie)
{
	int validated = check_logical_attr(nv, "valid"s);
	if (!validated)
		return validate(nv, !newbie);
	return validated >> 1;
}

/// __________________________________________________
/// Check "lat_valid" and "lon_valid attributes of DataFrame are all true
bool check_valid(const DataFrame df, bool newbie)
{
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
	using t_type = std::remove_const_t<decltype(t)>;
	bool iscoords {false};
	bool warn {false};
	auto valid { validate_switch(t) };
	if constexpr (isNumericVector_v<t_type>) {
		iscoords = true;
		if (!std::all_of(valid.begin(), valid.end(), [](auto i){ return i; }))
			warn = true;
		static_cast<NumericVector>(t).attr("valid") = valid; 
	} else if constexpr (Is_DataFrame<t_type>) {
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
	bool valid = false;
	vector llcols { get_vec_attr<int>(df, "llcols"s) };
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

/// __________________________________________________
/// Create coords - S3 method as_coords.default()
//' @rdname coords 
// [[Rcpp::export(name = "as_coords.default")]]
NumericVector as_coords(NumericVector object, int fmt = 1)
{
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
	if (!check_valid(x))
		stop("Invalid coords! Conversion aborted.\n [Use review() to show invalid elements]");
	if (newtype != ct_current) {
		auto vd_out { convert_switch(x, newtype) };
		NumericVector nv_out { wrap(vd_out) };									// Copies output string
		nv_out.attr("class") = "coords";
		nv_out.attr("fmt") = fmt;
		nv_out.attr("valid") = x.attr("valid");
		nv_out.attr("latlon") = x.attr("latlon");
		nv_out.names() = x.names();
		return nv_out;
	} else { 
		warning("Returning x (fmt = %i, unchanged)", fmt);
		return x;
	}
}

/// __________________________________________________
/// Set latlon attribute on "coords" NumericVector and revalidate
//' @rdname coords
// [[Rcpp::export(name = "`latlon<-`")]]
NumericVector latlon(NumericVector cd, LogicalVector value)
{
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
	using enum CoordType;
	checkinherits(x, "coords"s);
	if(!x.size())
		stop("x has 0 length!");
	if (validate)
		if (!check_valid(x))
			warning("Formatting invalid coords!\n [Use review() to show invalid elements]");
	CoordType ct_current { get_coordtype(x) };
	CoordType ct_required { fmt ? get_coordtype(fmt) : ct_current };
	auto sv_out { format_switch_c(x, ct_required) };
	vector names{ get_vec_attr<string>(x, "names"s) };
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
	if (!check_valid(x))
		stop("Invalid waypoints! Conversion aborted.\n [Use review() to show invalid elements]");
	if(!valid_ll(x))
		stop("Invalid llcols attribute!");
	if (newtype != ct_current) {
		Waypoints wp{ x };
		auto vd_lat { wp.convert(newtype, true) };
		auto vd_lon { wp.convert(newtype, false) };
		auto llcols{ get_vec_attr<int>(x, "llcols") };
		for (auto& llcol : llcols)							// llcols to C++ zero-based indexing
			--llcol;
		auto namescol{ get_vec_attr<int>(x, "namescol") };
		auto names{ get_vec_attr<string>(x, "names") };
		auto row_names{ get_vec_attr<int>(x, "row.names") };
		auto validlat{ get_vec_attr<bool>(x, "validlat") };
		auto validlon{ get_vec_attr<bool>(x, "validlon") };

		auto llcol_it{ x.erase(llcols[0]) };
		x.insert(llcol_it, vd_lat);
		llcol_it = x.erase(llcols[1]);
		x.insert(llcol_it, vd_lon);

		x.attr("names") = names;
		x.attr("class") = vector{"waypoints", "data.frame"};
		x.attr("row.names") = row_names;
		x.attr("fmt") = fmt;
		x.attr("namescol") = namescol;
		for (auto& llcol : llcols)	// llcols to R one-based indexing
			++llcol;
		x.attr("llcols") = llcols;
		x.attr("validlat") = validlat;
		x.attr("validlon") = validlon;
	} else
		warning("Returning x (fmt = %i, unchanged)", fmt);
	return x;
}


/// __________________________________________________
/// Format waypoints - S3 method format.waypoints()
//' @rdname format
// [[Rcpp::export(name = "format.waypoints")]]
CharacterVector formatwaypoints(DataFrame x, bool usenames = true, bool validate = true, int fmt = 0)
{
	checkinherits(x, "waypoints"s);
	if(!x.nrows())
		stop("x has 0 rows!");
	if(!valid_ll(x))
		stop("Invalid llcols attribute!");
	if (validate)
		if (!check_valid(x))
			warning("Formatting invalid waypoints!");
	Waypoints wp{ x };
	auto required { fmt ? get_coordtype(fmt) : get_coordtype(x) };
	auto vs_lat { wp.format(required, true) };
	auto vs_lon { wp.format(required, false) };
	transform(vs_lat.begin(), vs_lat.end(), vs_lon.begin(), vs_lat.begin(), [](auto& latstr, auto& lonstr){ return latstr + "  " + lonstr; });
	if (usenames) {
		RObject names = getnames(x);
		if (!prefixwithnames(vs_lat, names))
			stop("Invalid \"namescol\" attribute!");
	}
	return wrap(vs_lat);
}

/// __________________________________________________
/// Validate waypoints - S3 method validate.waypoints()
//' @rdname validate
// [[Rcpp::export(name = "validate.waypoints")]]
DataFrame validatewaypoints(DataFrame x, bool force = true)
{
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
	--fmt;														// -> C++ array numbering
	constexpr auto spacing{ array{ 0, 2, 3 } };
	const auto llstring{ "Latitude"s + string(5 + spacing[fmt], ' ') + "Longitude"s };
	const auto u_string{ string(11 + spacing[fmt], '_') + "  "s + string(12 + spacing[fmt], '_') };
	
	return wrap(vector{
		string(width - llstring.length() - 1, ' ') + llstring,
		string(width - u_string.length(), ' ') + u_string
	});
}

/// __________________________________________________
/// Clone coords object from waypoints vector
//' @rdname coords
// [[Rcpp::export(name = "as_coords.waypoints")]]
NumericVector as_coordswaypoints(DataFrame object, bool which)
{
	checkinherits(object, "waypoints"s);
	NumericVector nv = object[get_vec_attr<int>(object, "llcols"s)[which ? 0 : 1] - 1];
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
