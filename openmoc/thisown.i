/** Rules for Python/C++ transferrable memory ownership */

%module thisown

/* A Cell owns the memory for each Surface it contains */
%pythonappend Surface::addSurface %{
  args[1].thisown = 0
%}

/* Python must free memory for each Surface that is not in a Cell */
%pythonappend Cell::removeSurface %{
  args[0].thisown = 1
%}

/* A Universe owns the memory for each Cell it contains */
%pythonappend Universe::addCell %{
  args[0].thisown = 0
%}

/* Python must free memory for each Cell that is not in a Universe */
%pythonappend Universe::removeCell %{
  args[0].thisown = 1
%}

/* A Lattice owns the memory for each Universe it contains */
%pythonappend Lattice::setUniverses %{
  for i in range(len(args[0])):
    for j in range(len(args[0][i])):
      args[0][i][j].thisown = 0
%}

/* Python must free memory for each Universe that is not in a Lattice */
%pythonappend Lattice::removeUniverse %{
  args[0].thisown = 1
%}

/* A TrackGenerator owns the memory for the Geometry it contains */
%pythonappend TrackGenerator::TrackGenerator %{
  args[0].thisown = 0
%}

/* A TrackGenerator owns the memory for the Geometry it contains */
%pythonappend TrackGenerator::setGeometry %{
  args[0].thisown = 0
%}

/* A Geometry owns the memory for the root Universe it contains */
%pythonappend Geometry::setRootUniverse %{
  args[0].thisown = 0
%}

/* A Solver owns the memory for the PolarQuadrature it contains */
%pythonappend Solver::setPolarQuadrature %{
  args[0].thisown = 0
%}

/* A Solver owns the memory for the TrackGenerator it contains */
%pythonappend Solver::setTrackGenerator %{
  args[0].thisown = 0
%}

/* A Solver owns the memory for the TrackGenerator it contains */
%pythonappend Solver::Solver %{
  track_generator.thisown = 0
%}

/* A CPUSolver owns the memory for the TrackGenerator it contains */
%pythonappend CPUSolver::CPUSolver %{
  track_generator.thisown = 0
%}

/* A GPUSolver owns the memory for the TrackGenerator it contains */
%pythonappend GPUSolver::GPUSolver %{
  track_generator.thisown = 0
%}

/* A VectorizedSolver owns the memory for the TrackGenerator it contains */
%pythonappend VectorizedSolver::VectorizedSolver %{
  track_generator.thisown = 0
%}
