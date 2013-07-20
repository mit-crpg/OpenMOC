/**
 * @file log.h
 * @brief Utility functions for writing log messages to the screen
 * @details Applies level-based logging to print formatted messages
 *          to the screen and to a log file. 
 * @author William Boyd (wboyd@mit.edu)
 * @date January 22, 2012
 *
 */

#ifndef LOG_H_
#define LOG_H_

#ifdef __cplusplus
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdexcept>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif


/**
 * @enum logLevels
 * @brief Logging levels characterize an ordered set of message types
 *        which may be printed to the screen.
 */


/**
 * @var logLevel
 * @brief Logging levels characterize an ordered set of message types
 *        which may be printed to the screen.
 */
typedef enum logLevels {
    /** A debugging message */
    DEBUG,

    /** An informational but verbose message */
    INFO,

    /** A brief progress update on run progress */
    NORMAL,

    /** A message of a single line of characters */
    SEPARATOR,

    /** A message centered within a line of characters */
    HEADER,

    /** A message sandwiched between two lines of characters */
    TITLE,

    /** A message for to warn the user */
    WARNING,

    /** A message to warn of critical program conditions */
    CRITICAL,

    /** A message containing program results */
    RESULT,

    /** A messsage for unit testing */
    UNITTEST,

    /** A message reporting error conditions */
    ERROR
} logLevel;


/**
 * @brief A function stub used to convert C++ exceptions into Python exceptions
 *        through SWIG.
 * @param msg a character array for the exception message
 */
extern void set_err(const char *msg);

void setOutputDirectory(char* directory);
const char* getOutputDirectory();
void setLogfileName(char* filename);
const char* getLogfileName();

void setSeparatorCharacter(char c);
char getSeparatorCharacter();
void setHeaderCharacter(char c);
char getHeaderCharacter();
void setTitleCharacter(char c);
char getTitleCharacter();
void setLineLength(int length);

void setLogLevel(const char* newlevel);
int getLogLevel();

void log_printf(logLevel level, const char *format, ...);
std::string createMultilineMsg(std::string level, std::string message);


#endif /* LOG_H_ */
