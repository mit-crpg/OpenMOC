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
#ifdef SWIG
#include "Python.h"
#endif
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
#include <omp.h>
#endif

#ifdef MPIx
#include <mpi.h>
#endif

#ifdef SWIG
#define printf PySys_WriteStdout
#endif

/**
 * @enum logLevels
 * @brief Logging levels characterize an ordered set of message types
 *        which may be printed to the screen.
 */
typedef enum logLevels{
  /** A debugging message */
  DEBUG,

  /** An informational but verbose message */
  INFO,

  /** An informational verbose message - printed by rank 0 process only */
  INFO_ONCE,

  /** A brief progress update on run progress */
  NORMAL,

  /** A brief progress update by node on run progress */
  NODAL,

  /** A message of a single line of characters */
  SEPARATOR,

  /** A message centered within a line of characters */
  HEADER,

  /** A message sandwiched between two lines of characters */
  TITLE,

  /** A message to warn the user */
  WARNING,

  /** A message to warn the user - to be printed by rank 0 process only */
  WARNING_ONCE,

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
 * @details This method is not defined in the C++ source. It is defined in the
 *          SWIG inteface files (i.e., openmoc/openmoc.i)
 * @param msg a character array for the exception message
 */
extern void set_err(const char *msg);

void initialize_logger();
void set_output_directory(char* directory);
const char* get_output_directory();
void set_log_filename(char* filename);
const char* get_log_filename();

void set_separator_character(char c);
char get_separator_character();
void set_header_character(char c);
char get_header_character();
void set_title_character(char c);
char get_title_character();
void set_line_length(int length);
void set_log_level(const char* new_level);
void set_log_level(int new_level);
int get_log_level();

void log_printf(logLevel level, const char *format, ...);
std::string create_multiline_msg(std::string level, std::string message);
#ifdef MPIx
void log_set_ranks(MPI_Comm comm);
#endif

#endif /* LOG_H_ */
