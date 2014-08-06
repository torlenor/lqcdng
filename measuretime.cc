#include <ctime>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

double gettime(){ // Fix for Mac OS X
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
	clock_serv_t cclock;
	mach_timespec_t spec;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &spec);
	mach_port_deallocate(mach_task_self(), cclock);

#else
	struct timespec spec;
	clockid_t types[] = { CLOCK_REALTIME, CLOCK_MONOTONIC, CLOCK_PROCESS_CPUTIME_ID,
			CLOCK_THREAD_CPUTIME_ID, (clockid_t) - 1 };
	clock_gettime (types[0], &spec);
	clock_gettime(CLOCK_REALTIME, &spec);
#endif

	return spec.tv_nsec/(double)1E9 + spec.tv_sec;
}
