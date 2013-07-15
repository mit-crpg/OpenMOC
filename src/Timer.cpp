#include "Timer.h"


std::map<std::string, double> Timer::_timer_splits;
std::vector<timespec*> Timer::_start_times;


/**
 * @brief Starts the timer.
 * @details This method is similar to starting a stopwatch.
 */
void Timer::startTimer() {

    timespec* start_time = new timespec;
    _start_times.push_back(start_time);	

    #ifdef __MACH__         /* For OSX */
    gettimeofday(start_time, NULL);
    #else                           /* For linux */
    clock_gettime(CLOCK_MONOTONIC, start_time);
    #endif
    _running = true;

    return;
}


/**
 * @brief Stops the timer.
 * @details This method is similar to stopping a stopwatch.
 */
void Timer::stopTimer() {

    if (_running) {

        timespec* end_time = new timespec;
	timespec* start_time = _start_times.back();
	
        #ifdef __MACH__                  /* For OSX */
        gettimeofday(end_time, NULL);
        #else                            /* For linux */
        clock_gettime(CLOCK_MONOTONIC, end_time);
        #endif

        _elapsed_time = diff(*start_time, *end_time);

	if (_start_times.empty())
	    _running = false;

	_start_times.pop_back();
	delete start_time;
	delete end_time;
    }

    return;
}


/**
 * @brief Records a message corresponding to a time for the current split.
 * @details When this method is called it assumes that the timer has been
 *          stopped and has the current time for the process corresponding 
 *          to the message.
 * @param msg a msg corresponding to this time split
 */
void Timer::recordSplit(const char* msg) {

    double time = getTime();
    std::string msg_string = std::string(msg);

    if (_timer_splits.find(msg_string) != _timer_splits.end())
        _timer_splits[msg_string] += time;
    else
        _timer_splits.insert(std::pair<std::string, double>(msg_string, time));
}


/**
 * @brief Returns the amount of time elapsed from startTimer to stopTimer. 
 * @return the elapsed time in seconds
 */
double Timer::getTime() {
    #ifdef __MACH__                  /* For OSX */
    return _elapsed_time * 1.0E-6;
    #else                           /* For Linux */
    return _elapsed_time * 1.0E-9;
    #endif
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
        return _timer_splits[msg_string];
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
    int msg_length, num_whitespaces;
    std::map<std::string, double>::iterator iter;

    for (iter = _timer_splits.begin(); iter != _timer_splits.end(); ++iter) {

        std::stringstream formatted_msg;

        curr_msg = (*iter).first;
        curr_split = (*iter).second;
        msg_length = curr_msg.length();

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


/**
 * @brief Clears the time split for this message and deletes 
 *        the message's entry in the Timer's splits log
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


/**
* @brief Helper function which computes the time between the values of
*        two timespec structs representing the start and end times for code.
* @param start timespec representing the start time
* @param end timespec representing the end time
*/
double Timer::diff(timespec start, timespec end) {

    double sec, delta;
    #ifdef __MACH__
    double usec;
    delta = end.tv_usec - startTimer.tv_usec;
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

    } 
    else {
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
