#ifdef __cplusplus
#include "log.h"
#endif

#define LOG_C


/**
 * @var log_level
 * @brief Minimum level of logging messages printed to the screen and log file.
 *        Default logging level is NORMAL.
 */
static logLevel log_level = NORMAL;

/**
 * @var logfile_name
 * @brief The name of the output logfile. By default this is an empty.
 *        and must be set for messages to be redirected to a logfile.
 */
static std::string logfile_name = "";


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
 *          logfile has been created or not.
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
 * @details The default is 67 which adds to the standard 80 characters when 
 *          accounting for the log level prepended to each log message.
 */
static int line_length = 67;


/**
 * @brief Sets the output directory for log files. 
 * @details If the directory does not exist, it creates it for the user.
 * @param directory a character array for the log file directory
 */
void setOutputDirectory(char* directory) {

    output_directory = std::string(directory);

    /* Check to see if directory exists - if not, create it */
    struct stat st;
    if (!stat(directory, &st) == 0) {
		mkdir(directory, S_IRWXU);
        mkdir((output_directory+"/log").c_str(), S_IRWXU);
    }

    return;
}


/**
 * @brief Returns the output directory for log files.
 * @return a character array for the log file directory
 */
const char* getOutputDirectory() {
    return output_directory.c_str();
}


/**
 * @brief Sets the name for the log file.
 * @param filename a character array for log filename
 */
void setLogfileName(char* filename) {
    logfile_name = std::string(filename);
}


/**
 * @brief Returns the log filename.
 * @return a character array for the log filename
 */
const char* getLogfileName() {
    return logfile_name.c_str();
}


/**
 * @brief Sets the character to be used when printing SEPARATOR type log messages.
 * @param c the character for SEPARATOR type log messages
 */
void setSeparatorCharacter(char c) {
    separator_char = c;
}


/**
 * @brief Returns the character used to format SEPARATOR type log messages.
 * @return the character used for SEPARATOR type log messages
 */
char getSeparatorCharacter() {
    return separator_char;
}

/**
 * @brief Sets the character to be used when printing HEADER type log messages.
 * @param c the character for HEADER type log messages
 */
void setHeaderCharacter(char c) {
    header_char = c;
}


/**
 * @brief Returns the character used to format HEADER type log messages.
 * @return the character used for HEADER type log messages
 */
char getHeaderCharacter() {
    return header_char;
}


/**
 * @brief Sets the character to be used when printing TITLE type log messages.
 * @param c the character for TITLE type log messages
 */
void setTitleCharacter(char c) {
    title_char = c;
}


/**
 * @brief Returns the character used to format TITLE type log messages.
 * @return the character used for TITLE type log messages
 */
char getTitleCharacter() {
    return title_char;
}


/**
 * @brief Sets the maximum line length for log messages. 
 * @details Messages longer than this amount will be broken up into 
            multiline messages.
 * @param length the maximum log message line length in characters
 */
void setLineLength(int length) {
    line_length = length;
}


/**
 * @brief Sets the minimum log message level which will be printed to the 
 *        console and to the log file.
 * @param newlevel the minimum logging level
 */
void setLogLevel(const char* newlevel) {

    if (strcmp("DEBUG", newlevel) == 0) {
        log_level = DEBUG;
        log_printf(INFO, "Logging level set to DEBUG");
    }
    else if (strcmp("INFO", newlevel) == 0) {
        log_level = INFO;
        log_printf(INFO, "Logging level set to INFO");
    }
    else if (strcmp("NORMAL", newlevel) == 0) {
        log_level = NORMAL;
        log_printf(INFO, "Logging level set to NORMAL");
    }
    else if (strcmp("SEPARATOR", newlevel) == 0) {
        log_level = HEADER;
        log_printf(INFO, "Logging level set to SEPARATOR");
    }
    else if (strcmp("HEADER", newlevel) == 0) {
        log_level = HEADER;
        log_printf(INFO, "Logging level set to HEADER");
    }
    else if (strcmp("TITLE", newlevel) == 0) {
        log_level = TITLE;
        log_printf(INFO, "Logging level set to TITLE");
    }
    else if (strcmp("WARNING", newlevel) == 0) {
        log_level = WARNING;
        log_printf(INFO, "Logging level set to WARNING");
    }
    else if (strcmp("CRITICAL", newlevel) == 0) {
        log_level = CRITICAL;
        log_printf(INFO, "Logging level set to CRITICAL");
    }
    else if (strcmp("RESULT", newlevel) == 0) {
        log_level = RESULT;
        log_printf(INFO, "Logging level set to RESULT");
    }
    else if (strcmp("UNITTEST", newlevel) == 0) {
        log_level = UNITTEST;
        log_printf(INFO, "Logging level set to UNITTEST");
    }
    else if (strcmp("ERROR", newlevel) == 0) {
        log_level = ERROR;
        log_printf(INFO, "Logging level set to ERROR");
    }

    return;
}


/**
 * @brief Return the minimum level for log messages printed to the screen.
 * @return the minimum level for log messages
 */
int getLogLevel(){
     return log_level;
}


/**
 * @brief Print a formatted message to the console. 
 * @details If logging level is ERROR, this function will throw a
 *          runtime exception
 * @param level the logging level for this message
 * @param *format variable list of C++ formatted arguments
 */
void log_printf(logLevel level, const char *format, ...) {

    char message[512];
    std::string msg_string;
    if (level >= log_level) {
    	va_list args;

        va_start(args, format);
        vsprintf(message, format, args);
        va_end(args);

    	/* Append the log level to the message */
    	switch (level) {
            case (DEBUG):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[  DEBUG  ]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

                /* Puts message on single line */
                else
                    msg_string = level_prefix + msg + "\n";

	        break;
            }
	    case (INFO):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[  INFO   ]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

                /* Puts message on single line */
                else
                    msg_string = level_prefix + msg + "\n";

	        break;
            }
	    case (NORMAL):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[  NORMAL ]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

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
                int size = strlen(message);
                int halfpad = (line_length - 4 - size) / 2;
                std::string pad1 = std::string(halfpad, header_char);
                std::string pad2 = std::string(halfpad + 
                                    (line_length - 4 - size) % 2, header_char);
                std::string prefix = std::string("[  HEADER ]  ");
                std::stringstream ss;
                ss << prefix << pad1 << "  " << message << "  " << pad2 << "\n";
                msg_string = ss.str();
	        break;
            }
	    case (TITLE):
            {
                int size = strlen(message);
                int halfpad = (line_length - size) / 2;
                std::string pad = std::string(halfpad, ' ');
                std::string prefix = std::string("[  TITLE  ]  ");
                std::stringstream ss;
                ss << prefix << std::string(line_length, title_char) << "\n";
                ss << prefix << pad << message << pad << "\n";
                ss << prefix << std::string(line_length, title_char) << "\n";
                msg_string = ss.str();
	        break;
            }
	    case (WARNING):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[ WARNING ]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

                /* Puts message on single line */
                else
                    msg_string = level_prefix + msg + "\n";

	        break;
            }
	    case (CRITICAL):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[ CRITICAL]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

                /* Puts message on single line */
                else
                    msg_string = level_prefix + msg + "\n";

	        break;
            }
	    case (RESULT):
                msg_string = std::string("[  RESULT ]  ") + message + "\n";
                break;
	    case (UNITTEST):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[ UNITTEST]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

                /* Puts message on single line */
                else
                    msg_string = level_prefix + msg + "\n";

	            break;
            }
            case (ERROR):
	    {
                va_start(args, format);
                vsprintf(message, format, args);
                va_end(args);
                set_err(message);
                throw std::runtime_error(message);
                break;
	    }
          }


        /* If this is our first time logging, add a header with date, time */
        if (!logging) {

            /* If output directory was not defined by user, then log file is
             * written to a "log" subdirectory. Create it if it doesn't exist */
            if (output_directory.compare(".") == 0) {
                struct stat st;
                if (!stat("log", &st) == 0)
                    mkdir("log", S_IRWXU);
            }

            /* Write the message to the output file */
            std::ofstream logfile;
            logfile.open ((output_directory + "/" + logfile_name).c_str(), 
			  std::ios::app); 

            /* Append date, time to the top of log output file */
            time_t rawtime;
            struct tm * timeinfo;
            time (&rawtime);
            timeinfo = localtime (&rawtime);
            logfile << "Current local time and date: " << asctime(timeinfo);
            logging = true;

            logfile.close();
        }

        /* Write the log message to the logfile */
        std::ofstream logfile;
        logfile.open((output_directory + "/" + logfile_name).c_str(), 
		     std::ios::app);
        logfile << msg_string;
        logfile.close();

        /* Write the log message to the shell */
        std::cout << msg_string;
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
std::string createMultilineMsg(std::string level, std::string message) {

    int size = message.length();

    std::string substring;
    int start = 0;
    int end = line_length;

    std::string msg_string;

    /* Loop over msg creating substrings for each line */
    while (end < size + line_length) {

        /* Append log level to the beginning of each line */
        msg_string += level;

        /* Begin multiline messages with ellipsis */
        if (start != 0)
            msg_string += "... ";

        /* Find the current full length substring for line*/
        substring = message.substr(start, line_length);

        /* Truncate substring to last complete word */
        if (end < size-1) {
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
        start = end + 1;
        end += line_length + 1;
    }

    /* Reset line length */
    line_length += 4;

    return msg_string;
}
