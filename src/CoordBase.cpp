/// __________________________________________________
/// CoordBase.cpp
/// __________________________________________________

// [[Rcpp::plugins(cpp23)]]

#include <Rcpp.h>
#include <cxxabi.h>
#include <memory>
#include <array>
#include <utility>

using namespace Rcpp;

using std::array;
using std::vector;
using std::string;
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
template<NumericVector_or_DataFrame T, typename U> 
inline vector<U> get_vec_attr(const T& t, const char* attrname)
{
//	fmt::print("@{} attr=\"{}\" {}\n", "get_vec_attr<T, U>(const T&, const char*)", attrname, t.hasAttribute(attrname) ? true : false);
	return t.hasAttribute(attrname) ? as<vector<U>>(t.attr(attrname)) : vector<U>{};
}

/// __________________________________________________
/// Return "fmt" attribute as int
inline int get_fmt_attribute(const NumericVector_or_DataFrame auto& t)
{
//	fmt::print("@get_fmt_attribute<T>(const NumericVector_or_DataFrame auto&); t: {}; fmt={}\n", demangle(typeid(t)), as<int>(t.attr("fmt")));
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
inline void checkinherits(const NumericVector_or_DataFrame auto& t, const char* classname)
{
//	fmt::print("@checkinherits(const NumericVector_or_DataFrame auto&, const char*); t: {}; classname \"{}\"\n", demangle(typeid(t)), classname);
	if (!t.inherits(classname)) stop("Argument must be a \"%s\" object", classname);
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
int name_pos_in_df(const DataFrame df, const char* name)
{
//	fmt::print("@name_pos_in_df(const DataFrame, const char*); name={}\n", name);
	vector names{ get_vec_attr<DataFrame, string>(df, "names") };
	if (!names.size())
		return -1;
	typedef decltype(names.size()) Tmp;
	Tmp i = 0;
	for (auto str : names ) {
		// fmt::print("@@name_pos_in_df(const DataFrame, const char*); testing: {}\n", str);
		if (!str_tolower(str).compare(name)) {
			// fmt::print("@@@name_pos_in_df(const DataFrame, const char*); found: {}\n", str);
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
inline const CoordType get_coordtype(const NumericVector_or_DataFrame auto& t)
{
//	fmt::print("@get_coordtype(const NumericVector_or_DataFrame auto&); t: {}\n", demangle(typeid(t)));
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
/// New Business  !!!!!!!!!!!!!!
/// __________________________________________________
/// FamousFiveNew0 class

/// __________________________________________________
/// Destructor
FamousFiveNew0::~FamousFiveNew0()
{
	_ctrsgn(typeid(*this), false);
}

/// __________________________________________________
/// __________________________________________________
/// FamousFive0 class

/// __________________________________________________
/// Destructor
FamousFive0::~FamousFive0()
{
	fmt::print("§~FamousFive0() "); _ctrsgn(typeid(*this), false);
}

/// __________________________________________________
/// New Business —— !!!!!!!!!!!!!!
/// __________________________________________________
/// CoordletNew class

/// __________________________________________________
/// Constructor of CoordletNew
template<DVecType T>
CoordletNew<T>::CoordletNew(T&& _dv, const vector<bool> _latlon) :
	ff { make_unique<FamousFiveNew<T>>() },
	dv { static_cast<T&&>(_dv) },
	latlon{ _latlon }
{
	 _ctrsgn(typeid(*this)); fmt::print("\t(T&&, const vector<bool>); &_dv[0]: {}, &dv[0]: {}\n", address(_dv[0]), address(dv[0]));
}

/// __________________________________________________
/// Format dv as SVecType for printing
template<DVecType T> template<SVecType U>
U CoordletNew<T>::format() const
{
	fmt::print("@CoordletNew<T>::format<U>() const; T: {}, U: {}\n", demangle(typeid(T)), demangle(typeid(U)));
	U sv_out{ std::move(vector<string>(dv.size())) };

	if constexpr (isDecDegVecString_v<U>)
		transform(dv.begin(), dv.end(), sv_out.begin(), [this](auto n)
			{
				return fmt::format("{:>{}.{}f}\u00B0", ff->get_decdeg(n), 11, 6);
			});	
	else if constexpr (isDegMinVecString_v<U>)
		transform(dv.begin(), dv.end(), sv_out.begin(), [this](auto n)
			{
				return fmt::format("{:>{}}\u00B0", abs(ff->get_deg(n)), 3) + 
					   fmt::format("{:0>{}.{}f}\u2032", fabs(ff->get_decmin(n)), 7, 4);
			});	
	else if constexpr (isDegMinSecVecString_v<U>)
		transform(dv.begin(), dv.end(), sv_out.begin(), [this](auto n)
			{
				return fmt::format("{:>{}}\u00B0", abs(ff->get_deg(n)), 3) +
					   fmt::format("{:0>{}}\u2032", abs(ff->get_min(n)), 2) +
					   fmt::format("{:0>{}.{}f}\u2033", fabs(ff->get_sec(n)), 5, 2);
			});	

	fmt::print("@ICoordletNew<T>::format<U>() const; latlon: {}\n", fmt::join(latlon, ", "));
	vector<bool>::const_iterator ll_it { latlon.begin() };
	const auto ll_size { latlon.size() };

	if constexpr (isDecDegVecString_v<U>) {
		fmt::print("@IICoordletNew<T>::format<U>() const; T: {}, if constexpr (isDecDegVecString_v<U>)\n", demangle(typeid(T)));
		const auto lambda1 = [&ll_it](auto& outstr, auto n){ return outstr + (*ll_it++ ? " lat" : " lon"); };
		const auto lambda2 = [&ll_it](auto& outstr, auto n){ return outstr + (*ll_it ? " lat" : " lon"); };

		if (ll_size > 1)
			transform(sv_out.begin(), sv_out.end(), dv.begin(), sv_out.begin(), lambda1);
		else
			if (ll_size == 1)	// uniform coords
				transform(sv_out.begin(), sv_out.end(), dv.begin(), sv_out.begin(), lambda2);

	} else if constexpr (isDegMinVecString_v<U> || isDegMinSecVecString_v<U>) {
		fmt::print("@IIICoordletNew<T>::format<U>() const; T: {}, if constexpr (isDegMinVecString_v<U> || isDegMinSecVecString_v<U>)\n", demangle(typeid(T)));
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

	return sv_out;
}

/// __________________________________________________
/// Validate CoordletNew::dv
template<DVecType T>
const vector<bool> CoordletNew<T>::validate() const
{
	fmt::print("@CoordletNew<T>::validate(); latlon: {}\n", fmt::join(latlon, ", "));
	vector<bool>::const_iterator ll_it{ latlon.begin() };
	auto ll_size { latlon.size() };
	auto valid = vector<bool>{};
	valid.assign(dv.size(), {false});

	transform(dv.begin(), dv.end(), valid.begin(), [this, &ll_it, &ll_size](auto n){
		return !((fabs(ff->get_decdeg(n)) > (ll_size && (ll_size > 1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(fabs(ff->get_decmin(n)) >= 60) ||
				(fabs(ff->get_sec(n)) >= 60));
	});

	if (all_of(valid.begin(), valid.end(), [](auto v) { return v;}))
		valid.assign({true});

	return valid;
}

/// __________________________________________________
/// Temporary —— to be deleted
template<DVecType T>
void CoordletNew<T>::report() const
{
	fmt::print("@CoordletNew<T>::report() const; &dv: {}, dv[0]: {}, &dv[0]: {}\n", address(dv), dv[0], address(dv[0]));
}


/// __________________________________________________
/// New Business  !!!!!!!!!!!!!!
/// __________________________________________________
/// CrdWptBaseNew class

/// __________________________________________________
/// Constructor
CrdWptBaseNew::CrdWptBaseNew()
{
	_ctrsgn(typeid(*this)); fmt::print("\t()\n");
}

/// __________________________________________________
/// Destructor
CrdWptBaseNew::~CrdWptBaseNew()
{
	_ctrsgn(typeid(*this), false);
}

/// __________________________________________________
/// New Business  !!!!!!!!!!!!!!
/// __________________________________________________
/// CoordsNew class
template<DVecType T>
CoordsNew<T>::CoordsNew(vector<double> nv, const vector<bool> latlon) :
	cdlt { CoordletNew<T>{ std::move(nv), latlon }}
{
	_ctrsgn(typeid(*this)); fmt::print("\t(vector<double>, const vector<bool>)\n");
}

/// __________________________________________________
/// Convert call entry point -- public
template<DVecType T>
void CoordsNew<T>::convert(CoordType newtype)
{
	fmt::print("@CoordsNew<T>::convert(CoordType); new type: {}\n", newtype);
//	Coordlet{ nv }.convert(newtype);
//	nv.attr("fmt") = coordtype_to_int(newtype);
}

/// __________________________________________________
/// Format call entry point -- public
template<DVecType T>
vector<string> CoordsNew<T>::format(CoordType required_type) const
{
	fmt::print("@CoordsNew<T>::format(CoordType); required type: {}\n", required_type);
	using enum CoordType;

	switch (required_type)
	{
		case decdeg:
			return cdlt.template format<DecDegVecString>();

		case degmin:
			return cdlt.template format<DegMinVecString>();

		case degminsec:
			return cdlt.template format<DegMinSecVecString>();

		default:
			stop("CoordsNew<T>::format(CoordType) my bad");
	}
}

/// __________________________________________________
/// Validation call entry point -- public
template<DVecType T>
const vector<bool> CoordsNew<T>::validate() const							//	¡¡¡—— NB return type -> const vector<bool> ——!!!
{
	fmt::print("@CoordsNew<T>::validate()\n");
	return cdlt.validate();
}

/// __________________________________________________
/// Temporary —— to be deleted
template<DVecType T>
void CoordsNew<T>::report() const
{
	fmt::print("@CoordsNew<T>::report() const\n");
	cdlt.report();
}


/// __________________________________________________
/// Make Coords<DVecType>
unique_ptr<CrdWptBaseNew> coordsmaker(NumericVector nv)
{
	fmt::print("@coordsmaker(NumericVector); {}, &nv {}, &nv[0] {}\n", get_coordtype(nv), address(nv), address(nv[0]));
	using enum CoordType;
	const auto latlon { get_vec_attr<NumericVector, bool>(nv, "latlon") };
	switch (get_coordtype(nv))
	{
		case decdeg:
			return make_unique<CoordsNew<DecDegVecDouble>>( DecDegVecDouble{ nv }, latlon);

		case degmin:
			return make_unique<CoordsNew<DegMinVecDouble>>( DegMinVecDouble{ nv }, latlon );

		case degminsec:
			return make_unique<CoordsNew<DegMinSecVecDouble>>( DegMinSecVecDouble{ nv }, latlon );

		default:
			stop("coordsmaker(NumericVector) my bad");
	}
}


/// __________________________________________________
/// Old Hat  !!!!!!!!!!!!!!
/// __________________________________________________
/// Coordlet class

/// __________________________________________________
/// Constructor of Coordlet
Coordlet::Coordlet(NumericVector _nv) :				// Needs to know CoordType in another way to intialise ff
	ff { switch_ff(_nv) },
	nv{ _nv },
	latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") }
{
//	fmt::print("§Coordlet(NumericVector, bool) ");  _ctrsgn(typeid(*this));
}

/// __________________________________________________
/// Switch CoordType for ff in Coordlet constructor
unique_ptr<FamousFive0> Coordlet::switch_ff(NumericVector nv)
{
//	fmt::print("@§Coordlet::switch_ff(NumericVector)\n");
	using enum CoordType;

	switch (get_coordtype(nv))
	{
		case decdeg:
			return make_unique<FamousFive<decdeg>>();

		case degmin:
			return make_unique<FamousFive<degmin>>();

		case degminsec:
			return make_unique<FamousFive<degminsec>>();

		default:
			stop("Coordlet::switch_ff(NumericVector) my bad");
	}
}

/// __________________________________________________
/// Switch CoordType to convert format
void Coordlet::convert(CoordType required_type)
{
//	fmt::print("@Coordlet::convert(CoordType); required_type: {}\n", required_type);
	using enum CoordType;

	const auto lambdd = [this](auto n){ return ff->get_decdeg(n); };
	const auto lambdm = [this](auto n){ return ff->get_deg(n) * 1e2 + ff->get_decmin(n); };
	const auto lambdms = [this](auto n){ return ff->get_deg(n) * 1e4 + ff->get_min(n) * 1e2 + ff->get_sec(n); };

	switch (required_type)
	{
		case decdeg:
			transform(nv.begin(), nv.end(), nv.begin(), lambdd);
			break;

		case degmin:
			transform(nv.begin(), nv.end(), nv.begin(), lambdm);
			break;

		case degminsec:
			transform(nv.begin(), nv.end(), nv.begin(), lambdms);
			break;

		default:
			stop("Coordlet::convert(CoordType) my bad");
	}
}

/// __________________________________________________
/// Switch CoordType to format nv
vector<string> Coordlet::format(CoordType required_type) const
{
//	fmt::print("@Coordlet::format(CoordType) const; required: {}\n", required_type);
	auto sv_out = vector<string>(nv.size());

	const auto lambdd = [this](auto n){ return fmt::format("{:>{}.{}f}\u00B0", ff->get_decdeg(n), 11, 6); };
	const auto lambdm = [this](auto n){ return fmt::format("{:>{}}\u00B0", abs(ff->get_deg(n)), 3) +
											   fmt::format("{:0>{}.{}f}\u2032", fabs(ff->get_decmin(n)), 7, 4);
											};
	const auto lambdms = [this](auto n){ return fmt::format("{:>{}}\u00B0", abs(ff->get_deg(n)), 3) +
												fmt::format("{:0>{}}\u2032", abs(ff->get_min(n)), 2) +
												fmt::format("{:0>{}.{}f}\u2033", fabs(ff->get_sec(n)), 5, 2);
											};
	using enum CoordType;
	switch (required_type)
	{
		case decdeg:
			transform(nv.begin(), nv.end(), sv_out.begin(), lambdd);
			break;

		case degmin:
			transform(nv.begin(), nv.end(), sv_out.begin(), lambdm);
			break;

		case degminsec:
			transform(nv.begin(), nv.end(), sv_out.begin(), lambdms);
			break;

		default:
			stop("Coordlet<CoordType>::format(CoordType) const my bad");
	}
	return sv_out;
}

/// __________________________________________________
/// Validate Coordlet::nv
const vector<bool> Coordlet::validate() const
{
//	fmt::print("@Coordlet::validate(); latlon: {}\n", fmt::join(latlon, ", "));
	vector<bool>::const_iterator ll_it{ latlon.begin() };
	auto ll_size { latlon.size() };
	auto valid = vector<bool>{};
	valid.assign(nv.size(), {false});

	transform(nv.begin(), nv.end(), valid.begin(), [this, &ll_it, &ll_size](auto n){
		return !((fabs(ff->get_decdeg(n)) > (ll_size && (ll_size > 1 ? *ll_it++ : *ll_it) ? 90 : 180)) ||
				(fabs(ff->get_decmin(n)) >= 60) ||
				(fabs(ff->get_sec(n)) >= 60));
	});

	if (all_of(valid.begin(), valid.end(), [](auto v) { return v;}))
		valid.assign({true});

	return valid;
}


/// __________________________________________________
/// Old Hat  !!!!!!!!!!!!!!
/// __________________________________________________
/// CrdWptBase class

/// __________________________________________________
/// Constructor
CrdWptBase::CrdWptBase(CoordType _ct) : ct { _ct }
{
//	fmt::print("§CrdWptBase(CoordType); {}", ct); _ctrsgn(typeid(*this));
}

/// __________________________________________________
/// Destructor
CrdWptBase::~CrdWptBase()
{
//	fmt::print("§~CrdWptBase(); {}", ct); _ctrsgn(typeid(*this), false);
}

/// __________________________________________________
/// Old Hat  !!!!!!!!!!!!!!
/// __________________________________________________
/// Coords class

/// __________________________________________________
/// Constructor
Coords::Coords(NumericVector nv) : CrdWptBase { get_coordtype(nv) }, nv{ nv }
{
//	fmt::print("§Coords(NumericVector); {}", ct); _ctrsgn(typeid(*this));
}

/// __________________________________________________
/// Convert call entry point -- public
void Coords::convert(CoordType newtype)
{
//	fmt::print("@Coords::convert(CoordType); current type: {}; new type: {}\n", ct, newtype);
	Coordlet{ nv }.convert(newtype);
	nv.attr("fmt") = coordtype_to_int(newtype);
}

/// __________________________________________________
/// Format call entry point -- public
vector<string> Coords::format(CoordType required_type) const
{
//	fmt::print("@Coords::format(CoordType); current type: {}; required type: {}\n", ct, required_type);
	using enum CoordType;
	vector sv_out{ Coordlet{ nv }.format(required_type) };
	if (decdeg == required_type)
		suffix_latlon(sv_out);
	else
		suffix_nesw(sv_out);
	return sv_out;
}

/// __________________________________________________
/// Add suffix of "N", "E", "S", "W"; or "(N/E)", "(S/W)"
void Coords::suffix_nesw(vector<string>& sv_out) const
{
//	fmt::print("@Coords::suffix_nesw(vector<string>& sv_out) const\n");
	const auto latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") };
	vector<bool>::const_iterator ll_it { latlon.begin() };
	const auto ll_size { latlon.size() };

	const auto lambda1 = [&ll_it](auto& outstr, auto n){ return outstr + cardpoint(n < 0, *ll_it++); };
	const auto lambda2 = [&ll_it](auto& outstr, auto n){ return outstr + cardpoint(n < 0, *ll_it); };
	const auto lambda3 = [](auto& outstr, auto n){ return outstr + cardi_b(n < 0); };

	if (ll_size > 1)
		transform(sv_out.begin(), sv_out.end(), nv.begin(), sv_out.begin(), lambda1);
	else
		if (ll_size == 1)	// uniform coords
			transform(sv_out.begin(), sv_out.end(), nv.begin(), sv_out.begin(), lambda2);
		else				// no latlon info
			transform(sv_out.begin(), sv_out.end(), nv.begin(), sv_out.begin(), lambda3);
}

/// __________________________________________________
/// Add suffix of "lat", "lon"
void Coords::suffix_latlon(vector<string>& sv_out) const
{
//	fmt::print("@Coords::suffix_latlon(vector<string>& sv_out) const\n");
	const auto latlon{ get_vec_attr<NumericVector, bool>(nv, "latlon") };
	vector<bool>::const_iterator ll_it { latlon.begin() };
	const auto ll_size { latlon.size() };

	const auto lambda1 = [&ll_it](auto& outstr, auto n){ return outstr + (*ll_it++ ? " lat" : " lon"); };
	const auto lambda2 = [&ll_it](auto& outstr, auto n){ return outstr + (*ll_it ? " lat" : " lon"); };

	if (ll_size > 1)
		transform(sv_out.begin(), sv_out.end(), nv.begin(), sv_out.begin(), lambda1);
	else
		if (ll_size == 1)	// uniform coords
			transform(sv_out.begin(), sv_out.end(), nv.begin(), sv_out.begin(), lambda2);
}

/// __________________________________________________
/// Validation call entry point -- public
const bool Coords::validate() const
{
//	fmt::print("@Coords::validate(); current type: {}\n", ct);
	auto valid = Coordlet{ nv }.validate();

	static_cast<NumericVector>(nv).attr("valid") = valid;
	return ( std::all_of(valid.begin(), valid.end(), [](auto i){ return i; } )
	);
}


/// __________________________________________________
/// __________________________________________________
/// TBC  !!!!!!!!!!!!!!
/// __________________________________________________
/// Waypoints class

/// __________________________________________________
/// Constructor
Waypoints::Waypoints(DataFrame df) :
	CrdWptBase { get_coordtype(df) }, df{ df },
	nvlat( df[get_vec_attr<DataFrame, int>(df, "llcols")[0] - 1] ), 
	nvlon( df[get_vec_attr<DataFrame, int>(df, "llcols")[1] - 1] )
{
//	fmt::print("§Waypoints(DataFrame); {}", ct); _ctrsgn(typeid(*this));
	nvlat.attr("fmt") = coordtype_to_int(ct);
	nvlon.attr("fmt") = coordtype_to_int(ct);
	nvlat.attr("latlon") = true;
	nvlon.attr("latlon") = false;
}

/// __________________________________________________
/// Destructor
Waypoints::~Waypoints()
{
//	fmt::print("§~Waypoints(); {}", ct); _ctrsgn(typeid(*this), false);
	nvlat.attr("latlon") = R_NilValue;
	nvlon.attr("latlon") = R_NilValue;
	nvlat.attr("fmt") = R_NilValue;
	nvlon.attr("fmt") = R_NilValue;
}

/// __________________________________________________
/// Convert call entry point -- public
void Waypoints::convert(CoordType newtype)
{
//	fmt::print("@ Waypoints::convert(CoordType); current type: {}; new type: {}\n", ct, newtype);
	Coordlet{ nvlat }.convert(newtype);
	Coordlet{ nvlon }.convert(newtype);
	df.attr("fmt") = coordtype_to_int(newtype);
}

/// __________________________________________________
/// Format call entry point -- public
vector<string> Waypoints::format(CoordType required_type) const
{
//	fmt::print("@Waypoints::format(CoordType); current type: {}; required type: {}\n", ct, required_type);
	using enum CoordType;

	vector sv_lat{ Coordlet{ nvlat }.format(required_type) };
	vector sv_lon{ Coordlet{ nvlon }.format(required_type) };
	if (decdeg != required_type) {
		suffix_nesw(sv_lat, true);
		suffix_nesw(sv_lon, false);
	}
	transform(sv_lat.begin(), sv_lat.end(), sv_lon.begin(), sv_lat.begin(), [](auto& latstr, auto& lonstr){return latstr + "  " + lonstr;});
	return sv_lat;
}

/// __________________________________________________
/// Add suffix of  "N", "S", "E", "W" if CoordType::degmin or CoordType::degminsec
void Waypoints::suffix_nesw(vector<string>& sv_out, bool latlon) const
{
//	fmt::print("@Waypoints::suffix_nesw(vector<string> sv_out) const; {}\n", latlon ? "lat" : "lon");
	transform(sv_out.begin(), sv_out.end(), (latlon ? nvlat : nvlon).begin(), sv_out.begin(), [latlon](auto& outstr, auto n){
		return outstr + cardpoint(n < 0, latlon); }
	);
}

/// __________________________________________________
/// Validation call entry point -- public
const bool Waypoints::validate() const
{
//	fmt::print("@Waypoints::validate(); current type: {}\n", ct);
	auto validlat = Coordlet{ nvlat }.validate();
	auto validlon = Coordlet{ nvlon }.validate();

	static_cast<DataFrame>(df).attr("validlat") = validlat;
	static_cast<DataFrame>(df).attr("validlon") = validlon;

	return (
		std::all_of(validlat.begin(), validlat.end(), [](auto i){ return i; }) &&
		std::all_of(validlon.begin(), validlon.end(), [](auto i){ return i; })
	);
}

/// __________________________________________________
/// Revised  !!!!!!!!!!!!!!
/// __________________________________________________
/// Validation functions

/// __________________________________________________
/// Check "valid" attribute of NumericVector all true
bool check_valid(const NumericVector nv)
{
	fmt::print("@check_valid(const NumericVector)\n");
	int validated = check_logical_attr(nv, "valid");
	if (!validated)
		return revalidate<NumericVector>(nv);
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
		return revalidate<DataFrame>(df);

	if (!(latvalidated >> 1))
		warning("Invalid latitude!");
	if (!(lonvalidated >> 1))
		warning("Invalid longitude!");
	return latvalidated >> 1 || lonvalidated >> 1;
}

/// __________________________________________________
/// Revalidate "coords" or "waypoints"
template<NumericVector_or_DataFrame T>
bool revalidate(const T t)
{
	fmt::print("@revalidate<NumericVector_or_DataFrame>(const T); T: {}\n", demangle(typeid(t)));

	if constexpr (std::is_same_v<NumericVector, T>) { 
		auto valid { coordsmaker(t)->validate() };
		static_cast<NumericVector>(t).attr("valid") = valid; 
		if (!std::all_of(valid.begin(), valid.end(), [](auto i){ return i; }))
			warning("Revalidation found invalid coords!");
		else
			warning("Coords revalidated!");
	}

	if constexpr (std::is_same_v<DataFrame, T>) { 
		if (!Waypoints{ t }.validate())
			warning("Revalidation found invalid waypoints!");
		else
			warning("Waypoints revalidated!");
	}

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
		transform(llcols.begin(), llcols.end(), llcols.begin(), [](auto x){ return --x; });
		if (is_item_in_df(df, llcols[0]) && is_item_in_df(df, llcols[1]) && llcols[0] != llcols[1])
			if (is<NumericVector>(df[llcols[0]]) && is<NumericVector>(df[llcols[1]]))
				valid = true;
	}
	return valid;
}

/// __________________________________________________
/// In Progress  !!!!!!!!!!!!!!
/// __________________________________________________
/// Exported functions

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
	constexpr const char* mov { "moved" };
	constexpr const char* cpy { "copied" };

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

// Move construct CoordletNew<DecDegVecDouble> from DecDegVecDouble —— moves content  from dv to clt, maintains address
	fmt::print("{}V@movit(NumericVector); Constructing CoordletNew<DecDegVecDouble> by moving dv\n", exportstr);
	auto clt = CoordletNew<DecDegVecDouble>{ std::move(dv), vector<bool>{} };
	clt.report();
	fmt::print("{}VI@movit(NumVec); &dv {}, dv[0] {}, &dv[0] {}\n", exportstr, address(dv), "undefined", address(dv[0]));

// Re-assign to vector<double> nv from NumericVector —— copies "object" content to new address
	nv = as<vector<double>>(object);
	fmt::print("{}VII@movit(NumVec); &nv {}, nv[0] {}, &nv[0] {}\n", exportstr, address(nv), nv[0], address(nv[0]));

// Construct CoordletNew<DecDegVecDouble> indirectly from vector<double>; using *, new and delete to observe destruction
	fmt::print("{}VIII@movit(NumericVector); Making ptr1\n", exportstr);
	const auto* ptr1 = new CoordletNew<DecDegVecDouble>{ std::move(nv), vector<bool>{} };	// std::move() needed here—no copy elision.
	ptr1->report();
	fmt::print("{}IX@movit(NumVec); &nv {}, nv[0] {}, &nv[0] {}\n", exportstr, address(nv), "undefined", address(nv[0]));
	delete ptr1;

// Re-assign to vector<double> vn from NumericVector —— copies "object" content to new address (may re-use &dv[0] from deleted CoordletNew)
	nv = as<vector<double>>(object);
	fmt::print("{}X@movit(NumVec); &nv {}, nv[0] {}, &nv[0] {}\n", exportstr, address(nv), nv[0], address(nv[0]));

// Construct CoordletNew<DecDegVecDouble> directly from DecDegVecDouble; using *, new and delete to observe destruction
	fmt::print("{}XI@movit(NumericVector); Making ptr2\n", exportstr);
	const auto* ptr2 = new CoordletNew<DecDegVecDouble>{ DecDegVecDouble{ std::move(nv) }, vector<bool>{} };
	ptr2->report();
	fmt::print("{}XI@movit(NumVec); &nv {}, nv[0] {}, &nv[0] {}\n", exportstr, address(nv), "undefined", address(nv[0]));
	delete ptr2;
	fmt::print("{}XI@movit(NumericVector); ptr2 deleted, clt to be deleted on exit\n", exportstr);
	return object;
}

/// __________________________________________________
/// Dummy Function for Testing Only 	¡¡¡ ——— Temporary to Be Archived ——— !!!
//' @rdname cords
// [[Rcpp::export(name = "coordle")]]
NumericVector CoordsNewTest(NumericVector object)
{
	CoordType type = get_coordtype(object);
	fmt::print("{}@CoordsNewTest(NumericVector); fmt {}, CoordType {}\n", exportstr, get_fmt_attribute(object), type);

	auto Coords_ptr { unique_ptr<CrdWptBaseNew>{ coordsmaker(object) } };
	Coords_ptr->report();
	Coords_ptr->convert(CoordType::decdeg);
	Coords_ptr->validate();
	
	return object;
}

/// __________________________________________________
/// Create coords - S3 method as_coords.default()
//' @rdname coords 
// [[Rcpp::export(name = "as_coords.default")]]
NumericVector as_coords(NumericVector object, int fmt = 1)
{
//	fmt::print("{}@as_coords(NumericVector, int); fmt={}\n", exportstr, fmt);
	object.attr("fmt") = fmt;
	auto valid = coordsmaker(object)->validate();
	object.attr("valid") = valid;
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
	if (newtype != type)
		Coords{ x }.convert(newtype);
	else
		Rcout << "\t—— fmt out == fmt in! ——\n\n";
	return x;
}

/// __________________________________________________
/// Set latlon attribute on "coords" NumericVector and revalidate
//' @rdname coords
// [[Rcpp::export(name = "`latlon<-`")]]
NumericVector latlon(NumericVector cd, LogicalVector value)
{
	fmt::print("{}@latlon(NumericVector, LogicalVector)\n", exportstr);
	checkinherits(cd, "coords");
	if (value.size() != cd.size() && value.size() != 1)
		stop("value must be either length 1 or length(cd)");
	else
		cd.attr("latlon") = value;
	auto valid { coordsmaker(cd)->validate() };
	cd.attr("valid") = valid; 
	return cd;
}

/// __________________________________________________
/// Format coords - S3 method format.coords()
//' @rdname format
// [[Rcpp::export(name = "format.coords")]]
CharacterVector formatcoords(NumericVector x, bool usenames = true, bool validate = true, int fmt = 0)
{
	fmt::print("{}@formatcoords(NumericVector, bool, bool, int); usenames: {}, validate: {}, fmt: {}\n", exportstr, usenames, validate, fmt);
	fmt::print("{}@Iformatcoords(NumericVector, bool, bool, int); &x {}, &x[0] {}\n", exportstr, address(x), address(x[0]));

	checkinherits(x, "coords");
	if(!x.size())
		stop("x has 0 length!");
	if (validate)
		if (!check_valid(x))
			warning("Formatting invalid coords!");

	CoordType ct_current { get_coordtype(x) };
	CoordType ct_required { fmt ? get_coordtype(fmt) : ct_current };
	vector sv_out{ coordsmaker(x)->format(ct_required) };
	vector names{ get_vec_attr<NumericVector, string>(x, "names") };
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
	fmt::print("{}@validatecoords(const NumericVector, const bool); force: {}\n", exportstr, force);
	checkinherits(x, "coords");
	if (force)	{			
		auto valid { coordsmaker(x)->validate() };
		if (!std::all_of(valid.begin(), valid.end(), [](auto i){ return i; })) {
			warning("Validation of coords failed in Mimiland!");
			static_cast<NumericVector>(x).attr("valid") = valid;
		}
	}
	if (!check_valid(x))
		warning("Invalid coords!");
	return x;
}

/// __________________________________________________
/// Create waypoints - S3 method as_waypoints.default()
//' @rdname waypoints
// [[Rcpp::export(name = "as_waypoints.default")]]
DataFrame as_waypoints(DataFrame object, int fmt = 1)
{
//	fmt::print("{}@as_waypoints(DataFrame, int); fmt={}\n", exportstr, fmt);
	object.attr("fmt") = fmt;
	int namescol = 0;
	if (!object.hasAttribute("namescol")) {
		namescol = name_pos_in_df(object, "name");
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
	if (newtype != type)
		Waypoints{ x }.convert(newtype);
	else
		Rcout << "\t—— fmt out == fmt in! ——\n\n";
	return x;
}

/// __________________________________________________
/// Format waypoints - S3 method format.waypoints()
//' @rdname format
// [[Rcpp::export(name = "format.waypoints")]]
CharacterVector formatwaypoints(DataFrame x, bool usenames = true, bool validate = true, int fmt = 0)
{
//	fmt::print("{}@formatwaypoints(DataFrame, bool, bool, int); usenames: {}, validate: {}, fmt: {}\n", exportstr, usenames, validate, fmt);
	checkinherits(x, "waypoints");
	if(!x.nrows())
		stop("x has 0 rows!");
	if(!valid_ll(x))
		stop("Invalid llcols attribute!");
	if (validate)
		if (!check_valid(x))
			warning("Formatting invalid waypoints!");
	vector sv_out{ Waypoints{ x }.format(fmt ? get_coordtype(fmt) : get_coordtype(x)) };

	if (usenames) {
		RObject names = getnames(x);
		if (!prefixwithnames(sv_out, names))
			stop("Invalid \"namescol\" attribute!");
	}
	return wrap(sv_out);
}

/// __________________________________________________
/// Validate waypoints - S3 method validate.waypoints()
//' @rdname validate
// [[Rcpp::export(name = "validate.waypoints")]]
DataFrame validatewaypoints(DataFrame x, bool force = true)
{
//	fmt::print("{}@validatewaypoints(DataFrame, bool); force: {}\n", exportstr, force);
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
/// Latitude and longitude headers for S3 print.waypoint()
//' @rdname format
// [[Rcpp::export]]
CharacterVector ll_headers(int width, int fmt)
{
//	fmt::print("{}@ll_headers(int, int); width={}, fmt={}\n", exportstr, width, fmt);
	--fmt;  //	  to C++ array numbering
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
//	fmt::print("{}@as_coord(DataFrame); which: {}\n", exportstr, which ? "lat" : "lon");
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

#if DEBUG > 2
.
#endif
