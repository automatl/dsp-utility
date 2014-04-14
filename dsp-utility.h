#ifndef TOMATL_DSP_UTILITY
#define TOMATL_DSP_UTILITY

#ifndef NULL
#define NULL 0
#endif

#ifndef TOMATL_DELETE(x)
	#define TOMATL_DELETE(x) if (x != NULL) { delete x; x = NULL;}
#endif

#ifndef TOMATL_BRACE_DELETE(x)
	#define TOMATL_BRACE_DELETE(x) if (x != NULL) { delete[] x; x = NULL;}
#endif

#ifndef TOMATL_BOUND_VALUE(x, minval, maxval)
	#define TOMATL_BOUND_VALUE(x, minval, maxval) std::min((maxval), std::max((x), (minval)))
#endif

#ifndef TOMATL_TO_DB(x)
	#define TOMATL_TO_DB(x) 20 * std::log10(x)
#endif

#ifndef TOMATL_FROM_DB(x)
	#define TOMATL_FROM_DB(x) std::pow(2, x / 6)
#endif

#include "spsc_queue.h"
#include "Coord.h"
#include "Scaling.h"
#include "EnvelopeWalker.h"
#include "GonioCalculator.h"
#include "SimpleWindowedDft.h"
#include "SpectroCalculator.h"

#endif
