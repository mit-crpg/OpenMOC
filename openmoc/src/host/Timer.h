/**
 * @file Timer.h
 * @brief The Timer class.
 * @date January 2, 2012
 * @author  William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef TIMER_H_
#define TIMER_H_

#ifdef __cplusplus
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <map>
#include <string>
#include "log.h"
#endif

#ifdef __MACH__         /* For OSX */
#define timespec timeval
#endif


/**
 * @class Timer Timer.h "openmoc/src/host/Timer.cpp"
 * @brief The timer class is for timing and profiling regions of code.
 * @details The timer outputs running time in seconds but has resolution
 *          of microseconds on OSX and nanoseconds on Linux machines.
 */
class Timer {

private:
    /** A timespec struct with the current split's start time */
    timespec _start_time;
    /** A timespec struct with the current split's end time */
    timespec _end_time;
    /** The time elapsed (seconds) for the current split */
    float _elapsed_time;
    /** Whether or not the timer is running for the current split */
    bool _running;
    /** A vector of the times and messages for each split */
    static std::map<std::string, double> _timer_splits;

    /**
     * @brief Assignment operator for static referencing of the Timer.
     * @param & the Timer static class object
     * @return a pointer to the Timer static class object
     */
    Timer &operator=(const Timer &) { return *this; }

    /**
     * @brief Timer constructor.
     * @param & The Timer static reference pointer.
     */
    Timer(const Timer &) { }

    double diff(timespec start, timespec end);

public:
    /**
     * @brief Constructor sets the current split elapsed time to zero.
     */
    Timer() {
        _running = false;
        _elapsed_time = 0;
    }


    /**
     * @brief Destructor
     */
    virtual ~Timer() { }


    /**
     * @brief Returns a static instance of the Timer class.
     * @return a pointer to the static Timer class
     */
    static Timer *Get() {
        static Timer instance;
	return &instance;
    }

    void startTimer();
    void stopTimer();
    void resetTimer();
    void restartTimer();
    void recordSplit(std::string msg);
    void recordSplit(const char* msg);
    double getTime();
    double getSplit(std::string msg);
    double getSplit(const char* msg);
    void printSplits();
    void clearSplits();
};

#endif /* TIMER_H_ */
