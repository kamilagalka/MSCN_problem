#pragma once
#include <Windows.h>

class CTimer
{
public:
	CTimer();
	~CTimer();

	void start();
	void stop();
	double getSecs();

private:
	double d_secs;
	LARGE_INTEGER li_start;
	LARGE_INTEGER li_end;
	LARGE_INTEGER li_freq;
};

inline CTimer::CTimer()
{
}

inline CTimer::~CTimer()
{
}

inline void CTimer::start()
{
	QueryPerformanceFrequency(&li_freq);
	QueryPerformanceCounter(&li_start);
}

inline void CTimer::stop()
{
	QueryPerformanceCounter(&li_end);
}

inline double CTimer::getSecs()
{
	d_secs = li_end.QuadPart - li_start.QuadPart;
	d_secs = d_secs / li_freq.QuadPart;
	return d_secs;
}
