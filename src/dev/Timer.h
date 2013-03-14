/*
 * Timer.h
 *
 *  Created on: Jan 2, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <vector>
#include <cutil.h>
#include "../host/log.h"

#ifdef __MACH__		/* For OSX */
#define timespec timeval
#endif


class Timer {
protected:
	timespec _host_start_time, _host_end_time;
	cudaEvent_t _device_start_event, _device_end_event;
	float _host_elapsed_time, _device_elapsed_time;
	bool _host_running, _device_running;
	std::vector< std::pair<double, const char*> > _timer_splits;
public:
	Timer();
	virtual ~Timer();
	__host__ void startHostTimer();
	__host__ void stopHostTimer();
	__host__ void restartHostTimer();
	__host__ void resetHostTimer();
	__host__ void recordHostSplit(const char* msg);
	__host__ double getHostTime();
	__host__ void startDeviceTimer();
	__host__ void stopDeviceTimer();
	__host__ void restartDeviceTimer();
	__host__ void resetDeviceTimer();
	__host__ void recordDeviceSplit(const char* msg);
	__host__ double getDeviceTime();
	double diff(timespec start, timespec end);
	void printSplits();
};

#endif /* TIMER_H_ */
