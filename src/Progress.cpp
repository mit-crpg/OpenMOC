#include "Progress.h"

Progress::Progress(int num_iterations, std::string name, double interval) {

  _num_iterations = num_iterations;
  _name = name;
  _counter = 0;
  _last_interval = 0.0;
  _interval = interval;
}


Progress::~Progress() {
}


void Progress::incrementCounter() {

  #pragma omp critical
  {
    /* Output first result */
    if (_counter == 0) {
      std::string msg = "Progress " + _name + ": %4.2f %%";
      log_printf(NORMAL, msg.c_str(), 0.0);
    }

    _counter++;
    double progress = double(_counter) / _num_iterations;

    /* If next interval reached, output progress */
    if (progress >= _last_interval + _interval) {
      _last_interval += _interval;
      std::string msg = "Progress " + _name + ": %4.2f %%";
      log_printf(NORMAL, msg.c_str(), _last_interval*100.0);
    }
  }
}


void Progress::reset() {
  _counter = 0;
  _last_interval = 0.0;
}
