/// __________________________________________________
/// CoordBase.h
/// __________________________________________________

#ifndef COORDBASE_H_
#define COORDBASE_H_

/*
#define FMT_HEADER_ONLY
#include </opt/homebrew/Cellar/fmt/11.1.4/include/fmt/base.h>		// verbose path needs sorting!

/// Either use format_as() function or specialise the formatter struct template.
// #define USE_FORMAT_AS

int fmtdemo();
int fmtdemo2(int);
std::string fmtdemo3(int);

/// __________________________________________________
/// CoordType enum class

namespace CoordType_enum {

	enum class CoordType : char { decdeg, degmin, degminsec };
	
	inline const CoordType get_coordtype(int i);


#ifdef USE_FORMAT_AS
// Use format_as()

	auto format_as(CoordType ct) -> std::string;
#endif  // ifdef USE_FORMAT_AS

}

#ifndef USE_FORMAT_AS
// Use formatter struct template specialisation (can't be within namespace CoordType_enum)

template <> struct fmt::formatter<CoordType_enum::CoordType>: formatter<string_view> {
	// parse is inherited from formatter<string_view>.

	auto format(CoordType_enum::CoordType, format_context&) const
		-> format_context::iterator;
};
#endif  // ifndef USE_FORMAT_AS
*/

#endif  // COORDBASE_H_
