#ifdef __cplusplus
#include "log.h"
#endif

#define LOG_C


/**
 * @var log_level
 * @brief Minimum level of logging messages printed to the screen and log file.
 * @details The default logging level is NORMAL.
 */
static logLevel log_level = NORMAL;


/**
 * @var log_filename
 * @brief The name of the output log file. By default this is an empty character
 *        array and must be set for messages to be redirected to a log file.
 */
static std::string log_filename = "";


/**
 * @var output_directory
 * @brief The directory in which a "log" folder will be created for logfiles.
 * @details By default this is the same directory as that where the logger
 *          is invoked by an input file.
 */
static std::string output_directory = ".";


/**
 * @var logging
 * @brief A switch which is set to true once the first message is logged.
 * @details The logging switch is needed to indicate whether the output
 *          log file has been created or not.
 */
static bool logging = false;


/**
 * @var separator_char
 * @brief The character to use for SEPARATOR log messages. The default is "-".
 */
static char separator_char = '*';


/**
 * @var header_char
 * @brief The character to use for HEADER log messages. The default is "*".
 */
static char header_char = '*';


/**
 * @var title_char
 * @brief The character to use for TITLE log messages. The default is "*".
 */
static char title_char = '*';


/**
 * @var line_length
 * @brief The maximum line length for a log message.
 * @details The default is 67 which contributes to the standard 80 characters
 *          when accounting for the log level prepended to each log message.
 */
static int line_length = 67;


/**
 * @var log_error_lock
 * @brief OpenMP mutex lock for ERROR messages which throw exceptions
 */
static omp_lock_t log_error_lock;


/**
 * @brief Initializes the logger for use.
 * @details This should be immediately called when the logger is imported
 *          into Python and before any of its other routines are called. The
 *          routine initializes an OpenMP mutual exclusion lock which is used
 *          to preclude race conditions from occurring when an ERROR message
 *          is reported and program execution is terminated.
 */
void initialize_logger() {
  /* Initialize OpenMP mutex lock for ERROR messages with exceptions */
  omp_init_lock(&log_error_lock);
}


/**
 * @brief Sets the output directory for log files.
 * @details If the directory does not exist, it creates it for the user.
 * @param directory a character array for the log file directory
 */
void set_output_directory(char* directory) {

  output_directory = std::string(directory);
  std::string log_directory;

  /* Check to see if directory exists - if not, create it */
  struct stat st;
  if ((stat(directory, &st)) == 0) {
    log_directory = std::string("") + directory + std::string("/log");
    mkdir(log_directory.c_str(), S_IRWXU);
  }
}


/**
 * @brief Returns the output directory for log files.
 * @return a character array for the log file directory
 */
const char* get_output_directory() {
  return output_directory.c_str();
}


/**
 * @brief Sets the name for the log file.
 * @param filename a character array for log filename
 */
void set_log_filename(char* filename) {
  log_filename = std::string(filename);
}


/**
 * @brief Returns the log filename.
 * @return a character array for the log filename
 */
const char* get_log_filename() {
  return log_filename.c_str();
}


/**
 * @brief Sets the character to be used when printing SEPARATOR log messages.
 * @param c the character for SEPARATOR log messages
 */
void set_separator_character(char c) {
  separator_char = c;
}


/**
 * @brief Returns the character used to format SEPARATOR log messages.
 * @return the character used for SEPARATOR log messages
 */
char get_separator_character() {
  return separator_char;
}

/**
 * @brief Sets the character to be used when printing HEADER log messages.
 * @param c the character for HEADER log messages
 */
void set_header_character(char c) {
  header_char = c;
}


/**
 * @brief Returns the character used to format HEADER type log messages.
 * @return the character used for HEADER type log messages
 */
char get_header_character() {
  return header_char;
}


/**
 * @brief Sets the character to be used when printing TITLE log messages.
 * @param c the character for TITLE log messages
 */
void set_title_character(char c) {
  title_char = c;
}


/**
 * @brief Returns the character used to format TITLE log messages.
 * @return the character used for TITLE log messages
 */
char get_title_character() {
  return title_char;
}


/**
 * @brief Sets the maximum line length for log messages.
 * @details Messages longer than this amount will be broken up into
            multiline messages.
 * @param length the maximum log message line length in characters
 */
void set_line_length(int length) {
  line_length = length;
}


/**
 * @brief Sets the minimum log message level which will be printed to the
 *        console and to the log file.
 * @param new_level the minimum logging level as a character array
 */
void set_log_level(const char* new_level) {

  if (strcmp("DEBUG", new_level) == 0) {
    log_level = DEBUG;
    log_printf(INFO, "Logging level set to DEBUG");
  }
  else if (strcmp("INFO", new_level) == 0) {
    log_level = INFO;
    log_printf(INFO, "Logging level set to INFO");
  }
  else if (strcmp("NORMAL", new_level) == 0) {
    log_level = NORMAL;
    log_printf(INFO, "Logging level set to NORMAL");
  }
  else if (strcmp("SEPARATOR", new_level) == 0) {
    log_level = SEPARATOR;
    log_printf(INFO, "Logging level set to SEPARATOR");
  }
  else if (strcmp("HEADER", new_level) == 0) {
    log_level = HEADER;
    log_printf(INFO, "Logging level set to HEADER");
  }
  else if (strcmp("TITLE", new_level) == 0) {
    log_level = TITLE;
    log_printf(INFO, "Logging level set to TITLE");
  }
  else if (strcmp("WARNING", new_level) == 0) {
    log_level = WARNING;
    log_printf(INFO, "Logging level set to WARNING");
  }
  else if (strcmp("CRITICAL", new_level) == 0) {
    log_level = CRITICAL;
    log_printf(INFO, "Logging level set to CRITICAL");
  }
  else if (strcmp("RESULT", new_level) == 0) {
    log_level = RESULT;
    log_printf(INFO, "Logging level set to RESULT");
  }
  else if (strcmp("ERROR", new_level) == 0) {
      log_level = ERROR;
      log_printf(INFO, "Logging level set to ERROR");
  }
}


/**
 * @brief Return the minimum level for log messages printed to the screen.
 * @return the minimum level for log messages
 */
const char* get_log_level() {

  std::string level;

  switch (log_level) {
  case (DEBUG):
    level = "DEBUG";
    break;
  case (INFO):
    level = "INFO";
    break;
  case (NORMAL):
    level = "NORMAL";
    break;
  case (SEPARATOR):
    level = "SEPARATOR";
    break;
  case (HEADER):
    level = "HEADER";
    break;
  case (TITLE):
    level = "TITLE";
    break;
  case (WARNING):
    level = "WARNING";
    break;
  case (CRITICAL):
    level = "CRITICAL";
    break;
  case (RESULT):
    level = "RESULT";
    break;
  case (ERROR):
    level = "ERROR";
    break;
  }

  return level.c_str();
}


/**
 * @brief Print a formatted message to the console.
 * @details If the logging level is ERROR, this function will throw a
 *          runtime exception
 * @param level the logging level for this message
 * @param format variable list of C++ formatted arguments
 */
void log_printf(logLevel level, const char* format, ...) {

  if (level >= log_level) {
    va_list args;
    va_start(args, format);

    int size = 256;
    char* buffer = new char[size];
    std::string msg_string;
    int n = vsnprintf(buffer, size, format, args);

    /* Resize buffer as needed for strings longer than 256 characters */
    if (size <= n) {
      delete [] buffer;
      buffer = new char[n+1];
      n = vsnprintf(buffer, n+1, format, args);
    }

    std::string msg = std::string(buffer);
    va_end(args);

    /* Append the log level to the message */
    switch (level) {
    case (DEBUG):
      {
        std::string level_prefix = "[  DEBUG  ]  ";

        /* If message is too long for a line, split into many lines */
        if (int(msg.length()) > line_length)
          msg_string = create_multiline_msg(level_prefix, msg);

        /* Puts message on single line */
        else
          msg_string = level_prefix + msg + "\n";

        break;
      }
    case (INFO):
      {
        std::string level_prefix = "[  INFO   ]  ";

        /* If message is too long for a line, split into many lines */
        if (int(msg.length()) > line_length)
          msg_string = create_multiline_msg(level_prefix, msg);

        /* Puts message on single line */
        else
          msg_string = level_prefix + msg + "\n";

        break;
      }
    case (NORMAL):
      {
        std::string level_prefix = "[  NORMAL ]  ";

        /* If message is too long for a line, split into many lines */
        if (int(msg.length()) > line_length)
          msg_string = create_multiline_msg(level_prefix, msg);

        /* Puts message on single line */
        else
          msg_string = level_prefix + msg + "\n";

        break;
      }
    case (SEPARATOR):
      {
        std::string pad = std::string(line_length, separator_char);
        std::string prefix = std::string("[SEPARATOR]  ");
        std::stringstream ss;
        ss << prefix << pad << "\n";
        msg_string = ss.str();
        break;
      }
    case (HEADER):
      {
        int size = strlen(buffer);
        int halfpad = (line_length - 4 - size) / 2;
        std::string pad1 = std::string(halfpad, header_char);
        std::string pad2 = std::string(halfpad +
                           (line_length - 4 - size) % 2, header_char);
        std::string prefix = std::string("[  HEADER ]  ");
        std::stringstream ss;
        ss << prefix << pad1 << "  " << buffer << "  " << pad2 << "\n";
        msg_string = ss.str();
        break;
      }
    case (TITLE):
      {
        int size = strlen(buffer);
        int halfpad = (line_length - size) / 2;
        std::string pad = std::string(halfpad, ' ');
        std::string prefix = std::string("[  TITLE  ]  ");
        std::stringstream ss;
        ss << prefix << std::string(line_length, title_char) << "\n";
        ss << prefix << pad << buffer << pad << "\n";
        ss << prefix << std::string(line_length, title_char) << "\n";
        msg_string = ss.str();
        break;
      }
    case (WARNING):
      {
        std::string level_prefix = "[ WARNING ]  ";

        /* If message is too long for a line, split into many lines */
        if (int(msg.length()) > line_length)
          msg_string = create_multiline_msg(level_prefix, msg);

        /* Puts message on single line */
        else
          msg_string = level_prefix + msg + "\n";

        break;
      }
    case (CRITICAL):
      {
        std::string level_prefix = "[ CRITICAL]  ";

        /* If message is too long for a line, split into many lines */
        if (int(msg.length()) > line_length)
          msg_string = create_multiline_msg(level_prefix, msg);

        /* Puts message on single line */
        else
          msg_string = level_prefix + msg + "\n";

        break;
      }
    case (RESULT):
      {
        std::string level_prefix = "[  RESULT ]  ";

        /* If message is too long for a line, split into many lines */
        if (int(msg.length()) > line_length)
          msg_string = create_multiline_msg(level_prefix, msg);

        /* Puts message on single line */
        else
          msg_string = level_prefix + msg + "\n";

        break;
      }
    case (ERROR):
      {
        /* Create message based on runtime error stack */
        std::string level_prefix = "";

        /* If message is too long for a line, split into many lines */
        if (int(msg.length()) > line_length)
          msg_string = create_multiline_msg(level_prefix, msg);

        /* Puts message on single line */
        else
          msg_string = level_prefix + msg + "\n";
      }
    }

    /* If this is our first time logging, add a header with date, time */
    if (!logging) {

      /* If output directory was not defined by user, then log file is
       * written to a "log" subdirectory. Create it if it doesn't exist */
      if (output_directory.compare(".") == 0)
        set_output_directory((char*)".");

      /* Write the message to the output file */
      std::ofstream log_file;
      log_file.open((output_directory + "/log/" + log_filename).c_str(),
                   std::ios::app);

      /* Append date, time to the top of log output file */
      time_t rawtime;
      struct tm * timeinfo;
      time (&rawtime);
      timeinfo = localtime (&rawtime);
      log_file << "Current local time and date: " << asctime(timeinfo);
      logging = true;

      log_file.close();
    }

    /* Write the log message to the log_file */
    std::ofstream log_file;
    log_file.open((output_directory + "/log/" + log_filename).c_str(),
                  std::ios::app);
    log_file << msg_string;
    log_file.close();

    /* Write the log message to the shell */
    if (level == ERROR) {
      omp_set_lock(&log_error_lock);
      {
        throw std::logic_error(msg_string.c_str());
      }
      omp_unset_lock(&log_error_lock);
    }
    else
      printf("%s", msg_string.c_str());
  }
}


/**
 * @brief Breaks up a message which is too long for a single line into a
 *        multiline message.
 * @details This is an internal function which is called by log_printf and
 *          should not be called directly by the user.
 * @param level a string containing log level prefix
 * @param message a string containing the log message
 * @return a string with a formatted multiline message
 */
std::string create_multiline_msg(std::string level, std::string message) {

  int size = message.length();

  std::string substring;
  int start = 0;
  int end = line_length;

  std::string msg_string;

  /* Loop over msg creating substrings for each line */
  while (end <= size + line_length) {

    /* Append log level to the beginning of each line */
    msg_string += level;

    /* Begin multiline messages with ellipsis */
    if (start != 0)
      msg_string += "... ";

    /* Find the current full length substring for line*/
    substring = message.substr(start, line_length);

    /* Truncate substring to last complete word */
    if (end <= size) {
      int endspace = substring.find_last_of(" ");
      if (message.at(endspace+1) != ' ' &&
          endspace != int(std::string::npos)) {
        end -= line_length - endspace;
        substring = message.substr(start, end-start);
      }
    }

    /* concatenate substring to output message */
    msg_string += substring + "\n";

    /* Reduce line length to account for ellipsis prefix */
    if (start == 0)
      line_length -= 4;

    /* Update substring indices */
    start = end;
    end += line_length + 1;
  }

  /* Reset line length */
  line_length += 4;

  return msg_string;
}
