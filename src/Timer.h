/**
 * @file Timer.h
 * @brief The Timer class.
 * @date January 2, 2012
 * @author  William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef TIMER_H_
#define TIMER_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "log.h"
#include "constants.h"
#include <time.h>
#include <omp.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <map>
#include <vector>
#include <string>
#include <ios>
#include <unistd.h>
#endif

#ifdef MPIx
#include <mpi.h>
#endif


/**
 * @class Timer Timer.h "src/Timer.cpp"
 * @brief The Timer class is for timing and profiling regions of code.
 */
class Timer {

private:

  /** A vector of floating point start times at each inclusive level
   *  at which we are timing */
  static std::vector<double> _start_times;

  /** The time elapsed (seconds) for the current split */
  float _elapsed_time;

  /** Whether or not the Timer is running for the current split */
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

public:
  /**
   * @brief Constructor sets the current split elapsed time to zero.
   */
  Timer() {
    _running = false;
    _elapsed_time = 0;
  }

  /**
   * @brief Destructor.
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
  void recordSplit(const char* msg);
  double getTime();
  double getSplit(const char* msg);
  void printSplits();
  void clearSplit(const char* msg);
  void clearSplits();
  void processMemUsage(double& vm_usage, double& resident_set);
#ifdef MPIx
  void reduceTimer(MPI_Comm comm);
#endif
};

#endif /* TIMER_H_ */
