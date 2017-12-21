#include "StopWatch.h"
#include <windows.h>

//LARGE_INTEGER
//getFILETIMEoffset()
//{
//	SYSTEMTIME s;
//	FILETIME f;
//	LARGE_INTEGER t;
//
//	s.wYear = 1970;
//	s.wMonth = 1;
//	s.wDay = 1;
//	s.wHour = 0;
//	s.wMinute = 0;
//	s.wSecond = 0;
//	s.wMilliseconds = 0;
//	SystemTimeToFileTime(&s, &f);
//	t.QuadPart = f.dwHighDateTime;
//	t.QuadPart <<= 32;
//	t.QuadPart |= f.dwLowDateTime;
//	return (t);
//}
//
//int clock_gettime(int X, struct timeval *tv)
//{
//	LARGE_INTEGER           t;
//	FILETIME            f;
//	double                  microseconds;
//	static LARGE_INTEGER    offset;
//	static double           frequencyToMicroseconds;
//	static int              initialized = 0;
//	static BOOL             usePerformanceCounter = 0;
//
//	if (!initialized) {
//		LARGE_INTEGER performanceFrequency;
//		initialized = 1;
//		usePerformanceCounter = QueryPerformanceFrequency(&performanceFrequency);
//		if (usePerformanceCounter) {
//			QueryPerformanceCounter(&offset);
//			frequencyToMicroseconds = (double)performanceFrequency.QuadPart / 1000000.;
//		}
//		else {
//			offset = getFILETIMEoffset();
//			frequencyToMicroseconds = 10.;
//		}
//	}
//	if (usePerformanceCounter) QueryPerformanceCounter(&t);
//	else {
//		GetSystemTimeAsFileTime(&f);
//		t.QuadPart = f.dwHighDateTime;
//		t.QuadPart <<= 32;
//		t.QuadPart |= f.dwLowDateTime;
//	}
//
//	t.QuadPart -= offset.QuadPart;
//	microseconds = (double)t.QuadPart / frequencyToMicroseconds;
//	t.QuadPart = microseconds;
//	tv->tv_sec = t.QuadPart / 1000000;
//	tv->tv_usec = t.QuadPart % 1000000;
//	return (0);
//}

int clock_gettime(int, timespec *spec)      //C-file part
{
	__int64 wintime; GetSystemTimeAsFileTime((FILETIME*)&wintime);
	wintime -= 116444736000000000i64;  //1jan1601 to 1jan1970
	spec->tv_sec = wintime / 10000000i64;           //seconds
	spec->tv_nsec = wintime % 10000000i64 * 100;      //nano-seconds
	return 0;
}

CStopWatch::CStopWatch()
{
}

void CStopWatch::Start()
{
	clock_gettime(0, &ts_begin);
	mIsStarted = true;
}

void CStopWatch::Stop()
{
	if (!mIsStarted)
	{
		return;
	}

	clock_gettime(0, &ts_end);
	double time = (ts_end.tv_sec - ts_begin.tv_sec) +
		(ts_end.tv_nsec - ts_begin.tv_nsec) / 1e9;
	mTimeRunning += time;
}

void CStopWatch::Reset()
{
	mIsStarted = false;
	mTimeRunning = 0.0;
}

double CStopWatch::GetTime()
{
	return mTimeRunning;
}