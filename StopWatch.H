#pragma once

#include <ctime>

class StopWatch
{
	private:
		std::clock_t startTime;
		double totalTime;
		double intervallTime;
	public:
		StopWatch()
		{
			totalTime = 0;
		}
		void start()
		{
			startTime = std::clock();
		}
		void intervall()
		{
			intervallTime = ( std::clock() - startTime) / (double) CLOCKS_PER_SEC;
		}
		void stop()
		{
			totalTime += ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
		}
		double getIntervallTime()
		{
			return intervallTime;
		}
		double getTotalTime()
		{
			return totalTime;
		}
};
