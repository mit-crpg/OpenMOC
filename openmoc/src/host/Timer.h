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
#include <vector>
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
    std::vector< std::pair<double, const char*> > _timer_splits;

    double diff(timespec start, timespec end);

public:
    Timer();
    virtual ~Timer();
    void startTimer();
    void stopTimer();
    void resetTimer();
    void restartTimer();
    void recordSplit(const char* msg);
    double getTime();
    void printSplits();
};

#endif /* TIMER_H_ */
