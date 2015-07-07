/** Rules for Python/C++ transferrable memory ownership */

%module thisown

/* A Cell owns the memory for each Surface it contains */
%pythonappend Cell::addSurface %{
  # SWIG 3 
  if 'surface'in locals():
    surface = locals()['surface']
  # SWIG 2
  else:
    surface = locals()['args'][1]
  surface.thisown = 0
%}

/* Python must free memory for each Surface that is not in a Cell */
%pythonappend Cell::removeSurface %{
  # SWIG 3
  if 'surface' in locals():
    surface = locals()['surface']
  # SWIG 2
  else:
    surface = locals()['args'][0]
  surface.thisown = 1
%}

/* A Cell owns the memory for its Material/Universe fill */
%pythonappend Cell::setFill %{
  # SWIG 3
  if 'fill' in locals():
    fill = locals()['fill']
  # SWIG 2
  else:
    fill = locals()['args'][0]
  fill.thisown = 1
%}

/* A Universe owns the memory for each Cell it contains */
%pythonappend Universe::addCell %{
  # SWIG 3
  if 'cell' in locals():
    cell = locals()['cell']
  # SWIG 2
  else:
    cell = locals()['args'][0]
  cell.thisown = 0
%}

/* Python must free memory for each Cell that is not in a Universe */
%pythonappend Universe::removeCell %{
  # SWIG 3
  if 'cell' in locals():
    cell = locals()['cell']
  # SWIG 2
  else:
    cell = locals()['args'][0]
  # SWIG 2
  cell.thisown = 1
%}

/* A Lattice owns the memory for each Universe it contains */
%pythonappend Lattice::setUniverses %{
  # SWIG 3
  if 'universes' in locals():
    universes = locals()['universes']
  # SWIG 2
  else:
    universes = locals()['args'][0]

  for i in range(len(universes)):
    for j in range(len(universes[i])):
      universes[i][j].thisown = 0
%}

/* Python must free memory for each Universe that is not in a Lattice */
%pythonappend Lattice::removeUniverse %{
  # SWIG 3
  if 'universe' in locals():
    universe = locals()['universe']
  # SWIG 2
  else:
    universe = locals()['args'][0]
  universe.thisown = 1
%}

/* A TrackGenerator owns the memory for the Geometry it contains */
%pythonappend TrackGenerator::TrackGenerator %{
  # SWIG 3
  if 'geometry' in locals():
    geometry = locals()['geometry']
  # SWIG 2
  else:
    geometry = locals()['args'][0]
  geometry.thisown = 0
%}

/* A TrackGenerator owns the memory for the Geometry it contains */
%pythonappend TrackGenerator::setGeometry %{
  # SWIG 3
  if 'geometry' in locals():
    geometry = locals()['geometry']
  # SWIG 2
  else:
    geometry = locals()['args'][0]
  geometry.thisown = 0
%}

/* A Geometry owns the memory for the root Universe it contains */
%pythonappend Geometry::setRootUniverse %{
  # SWIG 3
  if 'root_universe' in locals():
    root_universe = locals()['root_universe']
  # SWIG 2
  else:
    root_universe = locals()['args'][0]
  root_universe.thisown = 0
%}

/* A Solver owns the memory for the PolarQuadrature it contains */
%pythonappend Solver::setPolarQuadrature %{
  # SWIG 3
  if 'polar_quad' in locals():
    polar_quad = locals()['polar_quad']
  # SWIG 2
  else:
    polar_quad = locals()['args'][0]
  polar_quad.thisown = 0
%}

/* A Solver owns the memory for the TrackGenerator it contains */
%pythonappend Solver::setTrackGenerator %{
  # SWIG 3
  if 'track_generator' in locals():
    track_generator = locals()['track_generator']
  # SWIG 2
  else:
    track_generator = locals()['args'][0]
  track_generator.thisown = 0
%}

/* A Solver owns the memory for the TrackGenerator it contains */
%pythonappend Solver::Solver %{
  # SWIG 3
  if 'track_generator' in locals():
    track_generator = locals()['track_generator']
  # SWIG 2
  else:
    track_generator = locals()['args'][0]
  track_generator.thisown = 0
%}

/* A CPUSolver owns the memory for the TrackGenerator it contains */
%pythonappend CPUSolver::CPUSolver %{
  # SWIG 3
  if 'track_generator' in locals():
    track_generator = locals()['track_generator']
  # SWIG 2
  else:
    track_generator = locals()['args'][0]
  track_generator.thisown = 0
%}

/* A GPUSolver owns the memory for the TrackGenerator it contains */
%pythonappend GPUSolver::GPUSolver %{
  # SWIG 3
  if 'track_generator' in locals():
    track_generator = locals()['track_generator']
  # SWIG 2
  else:
    track_generator = locals()['args'][0]
  track_generator.thisown = 0
%}

/* A VectorizedSolver owns the memory for the TrackGenerator it contains */
%pythonappend VectorizedSolver::VectorizedSolver %{
  # SWIG 3
  if 'track_generator' in locals():
    track_generator = locals()['track_generator']
  # SWIG 2
  else:
    track_generator = locals()['args'][0]
  track_generator.thisown = 0
%}
