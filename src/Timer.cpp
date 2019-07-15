#include "Timer.h"


std::map<std::string, double> Timer::_timer_splits;
std::vector<double> Timer::_start_times;


/**
 * @brief Starts the Timer.
 * @details This method is similar to starting a stopwatch.
 */
void Timer::startTimer() {

  double start_time = omp_get_wtime();
  _start_times.push_back(start_time);
  _running = true;

  return;
}


/**
 * @brief Stops the Timer.
 * @details This method is similar to stopping a stopwatch.
 */
void Timer::stopTimer() {

  if (_running) {

    double end_time = omp_get_wtime();
    double start_time = _start_times.back();

    _elapsed_time = end_time - start_time;

    if (_start_times.empty())
      _running = false;

    _start_times.pop_back();
  }

  return;
}


/**
 * @brief Records a message corresponding to a time for the current split.
 * @details When this method is called it assumes that the Timer has been
 *          stopped and has the current time for the process corresponding
 *          to the message.
 * @param msg a msg corresponding to this time split
 */
void Timer::recordSplit(const char* msg) {

  double time = getTime();
  std::string msg_string = std::string(msg);

  if (_timer_splits.find(msg_string) != _timer_splits.end())
    _timer_splits.at(msg_string) += time;
  else
    _timer_splits.insert(std::pair<std::string, double>(msg_string, time));
}


/**
 * @brief Returns the time elapsed from startTimer() to stopTimer().
 * @return the elapsed time in seconds
 */
double Timer::getTime() {
  return _elapsed_time;
}


/**
 * @brief Returns the time associated with a particular split.
 * @details If the split does not exist, returns 0.
 * @param msg the message tag for the split
 * @return the time recorded for the split (seconds)
 */
double Timer::getSplit(const char* msg) {

  std::string msg_string = std::string(msg);

  if (_timer_splits.find(msg_string) == _timer_splits.end())
    return 0.0;
  else
    return _timer_splits.at(msg_string);
}


/**
 * @brief Prints the times and messages for each split to the console.
 * @details This method will loop through all of the Timer's splits and print a
 *          formatted message string (80 characters in length) to the console
 *          with the message and the time corresponding to that message.
 */
void Timer::printSplits() {

  std::string curr_msg;
  double curr_split;
  std::map<std::string, double>::iterator iter;

  for (iter = _timer_splits.begin(); iter != _timer_splits.end(); ++iter) {

    std::stringstream formatted_msg;

    curr_msg = (*iter).first;
    curr_split = (*iter).second;

    curr_msg.resize(53, '.');
    formatted_msg << curr_msg;

    log_printf(RESULT, "%s%1.4E sec", formatted_msg.str().c_str(), curr_split);
  }
}


/**
 * @brief Clears the time split for this message and deletes the message's
 *        entry in the Timer's splits log.
 * @param msg the message tag for the split
 */
void Timer::clearSplit(const char* msg) {

  std::string msg_string = std::string(msg);

  if (_timer_splits.find(msg_string) == _timer_splits.end())
    return;
  else
    _timer_splits.erase(msg_string);
}


/**
 * @brief Clears all times split messages from the Timer.
 */
void Timer::clearSplits() {
  _timer_splits.clear();
}


//TODO This function is not really about time, it could be moved elsewhere.
/**
 * @brief Read memory usage file (on a HPC installation), and process it to make
 *        it more readable. Used for profiling.
 * @param vm_usage total use of virtual memory
 * @param resident_set total use of resident memory
 */
void Timer::processMemUsage(double& vm_usage, double& resident_set) {

   vm_usage     = 0.0;
   resident_set = 0.0;

   /* Open file containing memory info */
   std::ifstream stat_stream("/proc/self/stat", std::ios_base::in);

   /* Read in dummy data */
   std::string tmp;
   for (int i=0; i < 22; i++)
     stat_stream >> tmp;

   /* Read in virtual and resident memory */
   unsigned long vsize;
   stat_stream >> vsize;
   long rss;
   stat_stream >> rss;
   stat_stream.close();

   /* Calculate memory usage */
   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
   vm_usage = (double) vsize / 1024.0 / 1024.0;
   resident_set = rss * page_size_kb / 1024.0;
}


/**
 * @brief Transfer timer data across all domains.
 * @param comm a MPI communicator to transfer data
 */
#ifdef MPIx
void Timer::reduceTimer(MPI_Comm comm) {

  std::map<std::string, double>::iterator iter;
  for (iter = _timer_splits.begin(); iter != _timer_splits.end(); ++iter) {

    /* Collapse timing results down to one value for each category */
    double curr_split = (*iter).second;
    double total_split = 0;
    MPI_Reduce(&curr_split, &total_split, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    /* On the main node, average over the number of ranks, update result */
    if (fabs(total_split) > FLT_EPSILON) {
      int num_ranks;
      MPI_Comm_size(comm, &num_ranks);
      total_split /= num_ranks;
      (*iter).second = total_split;
    }
  }
}
#endif
