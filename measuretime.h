#ifndef MEASURETIME_H
#define MEASURETIME_H

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

double gettime(); // Fix for Mac OS X

#endif // MEASURETIME_H
