#ifndef TOMATL_DSP_UTILITY
#define TOMATL_DSP_UTILITY

#ifndef NULL
#define NULL 0
#endif

#ifndef TOMATL_DELETE
	#define TOMATL_DELETE(x) if (x != NULL) { delete x; x = NULL;}
#endif

#ifndef TOMATL_BRACE_DELETE
	#define TOMATL_BRACE_DELETE(x) if (x != NULL) { delete[] x; x = NULL;}
#endif

#ifndef TOMATL_BOUND_VALUE
	#define TOMATL_BOUND_VALUE(x, minval, maxval) std::min((maxval), std::max((x), (minval)))
#endif

#ifndef TOMATL_IS_IN_BOUNDS_INCLUSIVE
	#define TOMATL_IS_IN_BOUNDS_INCLUSIVE(x, minval, maxval) (x >= minval && x <= maxval)
#endif

#ifndef TOMATL_TO_DB
	#define TOMATL_TO_DB(x) (20. * std::log10(x))
#endif

#ifndef TOMATL_FROM_DB
	#define TOMATL_FROM_DB(x) (std::pow(2., x / 6.))
#endif

#ifndef TOMATL_FAST_DIVIDE_BY_255
	#define TOMATL_FAST_DIVIDE_BY_255(x) (((x) + 1 + ((x) >> 8)) >> 8);
#endif

// TODO: define for all platforms
#ifndef forcedinline
	#define forcedinline __forceinline
#endif

#include "spsc_queue.h"
#include "Coord.h"
#include "Scaling.h"
#include "EnvelopeWalker.h"
#include "GonioCalculator.h"
#include "SimpleWindowedDft.h"
#include "SpectroCalculator.h"
#include "FrequencyDomainGrid.h"

#endif
