#pragma once
#include <time.h>

class CStopWatch
{
public:
	CStopWatch();
	void Start();
	void Stop();
	void Reset();
	double GetTime();

private:
	bool mIsStarted = false;
	timespec ts_begin, ts_end;
	double mTimeRunning;
};
