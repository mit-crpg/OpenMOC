#include "Progress.h"
#include "Geometry.h"

Progress::Progress(int num_iterations, std::string name, double interval,
                   Geometry* geometry, bool mpi_comm) {

  if (geometry != NULL) {
    if (!geometry->isDomainDecomposed()) {
      geometry = NULL;
      mpi_comm = false;
    }
  }
  _mpi_comm = mpi_comm;
  _geometry = geometry;
  _num_iterations = num_iterations;
  _name = name;
  _counter = 0;
  _curr_interval = 0;
  int num_intervals = 1. / interval + 1;
  _intervals.resize(num_intervals);
  int interval_stride = interval * num_iterations;
  if (interval_stride == 0)
    interval_stride = 1;
  for (int i=0; i < num_intervals; i++)
    _intervals.at(i) = std::min(i * interval_stride, num_iterations-1);
  _intervals.at(num_intervals-1) = num_iterations-1;

}


Progress::~Progress() {
}


void Progress::incrementCounter() {

  int curr_count;
  bool found_interval = false;
  #pragma omp critical
  {
    curr_count = _counter++;
    while (_curr_interval < _intervals.size()) {
      if (curr_count != _intervals.at(_curr_interval))
        break;
      double count = curr_count;
      if (count != 0)
        count += 1;
      double num_iters = _num_iterations;
      double percent = count / num_iters * 100.0;
/*
#ifdef MPIx
      if (_mpi_comm)
        MPI_Barrier(_geometry->getMPICart());
#endif
*/
      std::string msg = "Progress " + _name + ": %4.2f %%";
      log_printf(NORMAL, msg.c_str(), percent);
      _curr_interval++;
      if (_curr_interval >= _intervals.size())
        break;
    }
#pragma omp flush (_counter, _curr_interval)
  }
}


void Progress::reset() {
  _counter = 0;
  _curr_interval = 0;
}
