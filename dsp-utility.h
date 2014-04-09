#ifndef TOMATL_DSP_UTILITY
#define TOMATL_DSP_UTILITY

#ifndef NULL
#define NULL 0
#endif

#ifndef TOMATL_DELETE
#define TOMATL_DELETE(x) if (x != NULL) { delete x; x = NULL;}
#endif

#include "spsc_queue.h"
#include "Coord.h"
#include "Scaling.h"
#include "EnvelopeWalker.h"
#include "GonioCalculator.h"
#include "SimpleWindowedDft.h"
#include "SpectroCalculator.h"

#endif
