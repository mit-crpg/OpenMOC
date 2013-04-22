/*
 * Timer.cpp
 *
 *  Created on: Jan 2, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 *  The timer class is for profiling code. It outputs running time in
 *  seconds but has resolution of microseconds on OSX and nanoseconds
 *  on Linux machines.
 */


#include "Timer.h"


/**
 * Timer class constructor
 */
Timer::Timer() {
	_host_running = false;
	_host_elapsed_time = 0;
}


/**
 * Default Timer destructor
 */
Timer::~Timer() { }


/**
 * Starts the host timer - similar to starting a stopwatch for
 * code running on the host
 */
__host__ void Timer::startHostTimer() {

	if (!_host_running) {
		#ifdef __MACH__ 	/* For OSX */
			gettimeofday(&_host_start_time, NULL);
		#else				/* For linux */
			clock_gettime(CLOCK_MONOTONIC, &_host_start_time);
		#endif
		_host_running = true;
	}
	return;
}


/**
 * Stops the host timer - similar to stopping a stopwatch for code
 * running on the host
 */
__host__ void Timer::stopHostTimer() {
	if (_host_running) {
		#ifdef __MACH__ /* For OSX */
			gettimeofday(_host_end_time, NULL);
		#else			/* For linux */
		  clock_gettime(CLOCK_MONOTONIC, &_host_end_time);
		#endif
		_host_running = false;
		_host_elapsed_time += diff(_host_start_time, _host_end_time);
	}
	return;
}


/**
 * Resets the host timer - similar to resetting a stopwatch for code
 * running on the host
 */
__host__ void Timer::resetHostTimer() {
	_host_elapsed_time = 0;
	_host_running = false;
}


/**
 * Restarts the host timer. The elapsed time will accumulate along with the
 * previous time(s) the timer was running. If the timer was already running
 * this function does nothing
 */
__host__ void Timer::restartHostTimer() {
	if (!_host_running) {
		_host_elapsed_time += diff(_host_start_time, _host_end_time);
		startHostTimer();
	}
}


/**
 * Records a message corresponding to a given time recorded by the host
 * timer. When this method is called it assumes that the timer has been
 * stopped and has the current time for the process corresponding to the
 * message
 * @param msg a msg corresponding to this time split
 */
__host__ void Timer::recordHostSplit(const char* msg) {
	_timer_splits.push_back(std::pair<double, const char*>(getHostTime(), msg));
}



/**
 * Returns the amount of time elapsed from startHostTimer to stopHostTimer
 * of the timer. If the timer is currently running, returns the time
 * relative from the host timer to the present time.
 * @return the elapsed time
 */
__host__ double Timer::getHostTime() {
	/* If the timer is not running */
	if (!_host_running) {
		#ifdef __MACH__		/* For OSX */
			return _host_elapsed_time * 1.0E-6;
		#else				/* For Linux */
			return _host_elapsed_time * 1.0E-9;
		#endif
	}

	/* If the timer is currently running */
	else {
		timespec temp;
		#ifdef __MACH__ 	/* For OSX */
			gettimeofday(&temp, NULL);
		#else				/* For Linux */
		  clock_gettime(CLOCK_MONOTONIC, &temp);
		#endif

		_host_elapsed_time += diff(_host_start_time, temp);

		#ifdef __MACH__		/* For OSX */
			return _host_elapsed_time * 1.0E-6;
		#else				/* For Linux */
			return _host_elapsed_time * 1.0E-9;
		#endif
	}
}


/**
 * Starts the device timer - similar to starting a stopwatch for code
 * running on the device
 */
__host__ void Timer::startDeviceTimer() {

	if (!_device_running) {
		CUDA_SAFE_CALL(cudaEventCreate(&_device_start_event));
		CUDA_SAFE_CALL(cudaEventCreate(&_device_end_event));
		CUDA_SAFE_CALL(cudaEventRecord(_device_start_event, 0));
		_device_running = true;
	}
	return;
}


/**
 * Stops the device timer - similar to stopping a stopwatch for code
 * running on the device
 */
__host__ void Timer::stopDeviceTimer() {
	if (this->_device_running) {
		CUDA_SAFE_CALL(cudaEventRecord(_device_end_event, 0));
		CUDA_SAFE_CALL(cudaEventSynchronize(_device_end_event));
		_device_running = false;
		CUDA_SAFE_CALL(cudaEventElapsedTime(&_device_elapsed_time,
 	 	 	 	 	 	 _device_start_event, _device_end_event));
	}
	return;
}


/**
 * Resets the device timer - similar to resetting a stopwatch for code
 * running on the device
 */
__host__ void Timer::resetDeviceTimer() {
	_device_elapsed_time = 0;
	_device_running = false;
}


/**
 * Restarts the device timer. The elapsed time will accumulate along with the
 * previous time(s) the timer was running. If the timer was already running
 * this function does nothing
 */
__host__ void Timer::restartDeviceTimer() {
	if (!_device_running) {
		_device_elapsed_time += _device_elapsed_time;
		startDeviceTimer();
	}
}


/**
 * Records a message corresponding to a given time recorded by the timer.
 * When this method is called it assumes that the timer has been stopped
 * and has the current time for the process corresponding to the message
 * @param msg a msg corresponding to this time split
 */
__host__ void Timer::recordDeviceSplit(const char* msg) {
	_timer_splits.push_back(std::pair<double, const char*>(getDeviceTime(),
																	msg));
}



/**
 * Returns the amount of time elapsed from startHostTimer to stopHostTimer
 * of the timer. If the timer is currently runnning, returns the time from
 * the timer startHostTimer to the present time.
 * @return the elapsed time
 */
__host__ double Timer::getDeviceTime() {
	/* If the timer is not running */
	if (!_device_running)
		return _device_elapsed_time * 1.0E-3;

	/* If the timer is currently running */
	else {
		log_printf(WARNING, "Unable to return the time for a currently"
					" running device process");
		return 0.0;
	}
}


/**
 * Helper function which computes the time between the values of
 * two timespec structs representing the start and end times for code
 * running on the host
 * @param start timespec representing the host start time
 * @param end timespec representing the host end time
 */
__host__ double Timer::diff(timespec start, timespec end) {
	double sec, delta;
	#ifdef __MACH__
		double usec;
		delta = end.tv_usec - startHostTimer.tv_usec;
	#else
		double nsec;
		delta = end.tv_nsec - start.tv_nsec;
	#endif

	if (delta < 0) {
		sec = end.tv_sec - start.tv_sec;
		#ifdef __MACH__
			usec = 1.0E6 + delta;
		#else
			nsec = 1.0E9 + delta;
		#endif

	} else {
		sec = end.tv_sec - start.tv_sec;
		#ifdef __MACH__
			usec = delta;
		#else
			nsec = delta;
		#endif
	}

	#ifdef __MACH__
		return (sec * 1.0E6 + usec);
	#else
		return(sec*1.0E9 + nsec);
	#endif
}


/**
 * This method will loop through all of the Timer's splits and print a
 * formatted message string (80 characters in length) to the console
 * with the message and the time corresponding to that message
 */
__host__ void Timer::printSplits() {

	const char* curr_msg;
	double curr_split;
	int msg_length, num_whitespaces;

	for (int i=0; i < (int)_timer_splits.size(); i++) {
		std::stringstream formatted_msg;

		curr_split = _timer_splits.at(i).first;
		curr_msg = _timer_splits.at(i).second;
		msg_length = strlen(curr_msg);

		/* Num whitespaces for formatting is:
		 * 80 max char - 13 char for logger - 13 for time - msg length */
		num_whitespaces = 80 - 13 - 11 - msg_length -3;

		formatted_msg << curr_msg;

		/* Add periods to format message to 80 characters length */
		for (int i=0; i < num_whitespaces; i++)
			formatted_msg << ".";

		log_printf(RESULT, "%s%.7f sec", formatted_msg.str().c_str(),
															curr_split);
	}
}
