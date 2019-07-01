/** Rules for Python/C++ transferrable memory ownership */

%module thisown

/* A Cell owns the memory for each Surface it contains */
%pythonappend Cell::addSurface %{
        # SWIG 3
        if 'surface' in locals():
            surface = locals()['surface']
        elif 'args' in locals() and 'surface' in locals()['args']:
            surface = locals()['args']['surface']
        elif 'kwargs' in locals() and 'surface' in locals()['kwargs']:
            surface = locals()['kwargs']['surface']

        # SWIG 2
        else:
            surface = locals()['args'][1]

        surface.thisown = False
%}

/* Python must free memory for each Surface that is not in a Cell */
%pythonappend Cell::removeSurface %{
        # SWIG 3
        if 'surface' in locals():
            surface = locals()['surface']
        elif 'args' in locals() and 'surface' in locals()['args']:
            surface = locals()['args']['surface']
        elif 'kwargs' in locals() and 'surface' in locals()['kwargs']:
            surface = locals()['kwargs']['surface']

        # SWIG 2
        else:
            surface = locals()['args'][0]

        surface.thisown = True
%}

/* A Cell owns the memory for its Material/Universe fill */
%pythonappend Cell::setFill %{
        # SWIG 3
        if 'fill' in locals():
            fill = locals()['fill']
        elif 'args' in locals() and 'fill' in locals()['args']:
            fill = locals()['args']['fill']
        elif 'kwargs' in locals() and 'fill' in locals()['kwargs']:
            fill = locals()['kwargs']['fill']

        # SWIG 2
        else:
            fill = locals()['args'][0]

        fill.thisown = False
%}

/* A Region owns the memory for each node it contains */
%pythonappend Region::addNode %{
        # SWIG 3
        if 'node' in locals():
            node = locals()['node']
        elif 'args' in locals() and 'node' in locals()['args']:
            node = locals()['args']['node']
        elif 'kwargs' in locals() and 'node' in locals()['kwargs']:
            node = locals()['kwargs']['node']

        # SWIG 2
        else:
            node = locals()['args'][0]

        node.thisown = False
%}

/* Python must free memory for each Node that is not in a Region  */
%pythonappend Region::removeNode %{
        # SWIG 3
        if 'node' in locals():
            node = locals()['node']
        elif 'args' in locals() and 'node' in locals()['args']:
            node = locals()['args']['node']
        elif 'kwargs' in locals() and 'node' in locals()['kwargs']:
            node = locals()['kwargs']['node']

        # SWIG 2
        else:
            node = locals()['args'][0]

        node.thisown = True
%}

 /* A Halfspace owns the memory for the Surface it contains */
%pythonappend Halfspace::Halfspace %{
        # SWIG 3
        if 'surface' in locals():
            surface = locals()['surface']
        elif 'args' in locals() and 'surface' in locals()['args']:
            surface = locals()['args']['surface']
        elif 'kwargs' in locals() and 'surface' in locals()['kwargs']:
            surface = locals()['kwargs']['surface']

        # SWIG 2
        else:
            surface = locals()['args'][1]

        surface.thisown = False
%}

/* A Universe owns the memory for each Cell it contains */
%pythonappend Universe::addCell %{
        # SWIG 3
        if 'cell' in locals():
            cell = locals()['cell']
        elif 'args' in locals() and 'cell' in locals()['args']:
            cell = locals()['args']['cell']
        elif 'kwargs' in locals() and 'cell' in locals()['kwargs']:
            cell = locals()['kwargs']['cell']

        # SWIG 2
        else:
            cell = locals()['args'][0]

        cell.thisown = False
%}

/* Python must free memory for each Cell that is not in a Universe */
%pythonappend Universe::removeCell %{
        # SWIG 3
        if 'cell' in locals():
            cell = locals()['cell']
        elif 'args' in locals() and 'cell' in locals()['args']:
            cell = locals()['args']['cell']
        elif 'kwargs' in locals() and 'cell' in locals()['kwargs']:
            cell = locals()['kwargs']['cell']

        # SWIG 2
        else:
            cell = locals()['args'][0]

        cell.thisown = True
%}

/* A Lattice owns the memory for each Universe it contains */
%pythonappend Lattice::setUniverses %{
        # SWIG 3
        if 'num_z' in locals():
            universes = locals()['num_z']
        elif 'args' in locals() and 'num_z' in locals()['args']:
            universes = locals()['args']['num_z']
        elif 'kwargs' in locals() and 'num_z' in locals()['kwargs']:
            universes = locals()['kwargs']['num_z']

        # SWIG 2
        else:
            universes = locals()['args'][0]

        for i in range(len(universes)):
            for j in range(len(universes[i])):
                for k in range(len(universes[i][j])):
                    universes[i][j][k].thisown = False
%}

/* A Lattice owns the memory for each Universe it contains */
%pythonappend Lattice::updateUniverse %{
        # SWIG 3
        if 'universe' in locals():
            universe = locals()['universe']
        elif 'args' in locals() and 'universe' in locals()['args']:
            universe = locals()['args']['universe']
        elif 'kwargs' in locals() and 'universe' in locals()['kwargs']:
            universe = locals()['kwargs']['universe']

        # SWIG 2
        else:
            universe = locals()['args'][3]

        universe.thisown = False
%}

/* Python must free memory for each Universe that is not in a Lattice */
%pythonappend Lattice::removeUniverse %{
        # SWIG 3
        if 'universe' in locals():
            universe = locals()['universe']
        elif 'args' in locals() and 'universe' in locals()['args']:
            universe = locals()['args']['universe']
        elif 'kwargs' in locals() and 'universe' in locals()['kwargs']:
            universe = locals()['kwargs']['universe']

        # SWIG 2
        else:
            universe = locals()['args'][0]

        universe.thisown = True
%}

/* A TrackGenerator owns the memory for the Geometry it contains */
%pythonappend TrackGenerator::TrackGenerator %{
        # SWIG 3
        if 'geometry' in locals():
            geometry = locals()['geometry']
        elif 'args' in locals() and 'geometry' in locals()['args']:
            geometry = locals()['args']['geometry']
        elif 'kwargs' in locals() and 'geometry' in locals()['kwargs']:
            geometry = locals()['kwargs']['geometry']

        # SWIG 2
        else:
            geometry = locals()['args'][0]

        geometry.thisown = False
%}

/* A TrackGenerator owns the memory for the Geometry it contains */
%pythonappend TrackGenerator::setGeometry %{
        # SWIG 3
        if 'geometry' in locals():
            geometry = locals()['geometry']
        elif 'args' in locals() and 'geometry' in locals()['args']:
            geometry = locals()['args']['geometry']
        elif 'kwargs' in locals() and 'geometry' in locals()['kwargs']:
            geometry = locals()['kwargs']['geometry']

        # SWIG 2
        else:
            geometry = locals()['args'][0]

        geometry.thisown = False
%}

/* A Geometry owns the memory for the root Universe it contains */
%pythonappend Geometry::setRootUniverse %{
        # SWIG 3
        if 'root_universe' in locals():
            root_universe = locals()['root_universe']
        elif 'args' in locals() and 'root_universe' in locals()['args']:
            root_universe = locals()['args']['root_universe']
        elif 'kwargs' in locals() and 'root_universe' in locals()['kwargs']:
            root_universe = locals()['kwargs']['root_universe']

        # SWIG 2
        else:
            root_universe = locals()['args'][0]

        root_universe.thisown = False
%}

/* A Solver owns the memory for the Quadrature it contains */
%pythonappend TrackGenerator::setQuadrature %{
        # SWIG 3
        if 'quadrature' in locals():
            quadrature = locals()['quadrature']
        elif 'args' in locals() and 'quadrature' in locals()['args']:
            quadrature = locals()['args']['quadrature']
        elif 'kwargs' in locals() and 'quadrature' in locals()['kwargs']:
            quadrature = locals()['kwargs']['quadrature']

        # SWIG 2
        else:
            quadrature = locals()['args'][0]

        quadrature.thisown = False
%}

/* A Solver owns the memory for the TrackGenerator it contains */
%pythonappend Solver::setTrackGenerator %{
        # SWIG 3
        if 'track_generator' in locals():
            track_generator = locals()['track_generator']
        elif 'args' in locals() and 'track_generator' in locals()['args']:
            track_generator = locals()['args']['track_generator']
        elif 'kwargs' in locals() and 'track_generator' in locals()['kwargs']:
            track_generator = locals()['kwargs']['track_generator']

        # SWIG 2
        else:
            track_generator = locals()['args'][0]

        track_generator.thisown = False
%}

/* A Solver owns the memory for the TrackGenerator it contains */
%pythonappend Solver::Solver %{
        # SWIG 3
        if 'track_generator' in locals():
            track_generator = locals()['track_generator']
        elif 'args' in locals() and 'track_generator' in locals()['args']:
            track_generator = locals()['args']['track_generator']
        elif 'kwargs' in locals() and 'track_generator' in locals()['kwargs']:
            track_generator = locals()['kwargs']['track_generator']

        # SWIG 2
        else:
            track_generator = locals()['args'][0]

        track_generator.thisown = False
%}

/* A CPUSolver owns the memory for the TrackGenerator it contains */
%pythonappend CPUSolver::CPUSolver %{
        # SWIG 3
        if 'track_generator' in locals():
            track_generator = locals()['track_generator']
        elif 'args' in locals() and 'track_generator' in locals()['args']:
            track_generator = locals()['args']['track_generator']
        elif 'kwargs' in locals() and 'track_generator' in locals()['kwargs']:
            track_generator = locals()['kwargs']['track_generator']

        # SWIG 2
        else:
            track_generator = locals()['args'][0]

        # Only disown if user specified optional track_generator parameter
        if track_generator:
            track_generator.thisown = False
%}

/* A GPUSolver owns the memory for the TrackGenerator it contains */
%pythonappend GPUSolver::GPUSolver %{
        # SWIG 3
        if 'track_generator' in locals():
            track_generator = locals()['track_generator']
        elif 'args' in locals() and 'track_generator' in locals()['args']:
            track_generator = locals()['args']['track_generator']
        elif 'kwargs' in locals() and 'track_generator' in locals()['kwargs']:
            track_generator = locals()['kwargs']['track_generator']

        # SWIG 2
        else:
            track_generator = locals()['args'][0]

        track_generator.thisown = False
%}

/* A VectorizedSolver owns the memory for the TrackGenerator it contains */
%pythonappend VectorizedSolver::VectorizedSolver %{
        # SWIG 3
        if 'track_generator' in locals():
            track_generator = locals()['track_generator']
        elif 'args' in locals() and 'track_generator' in locals()['args']:
            track_generator = locals()['args']['track_generator']
        elif 'kwargs' in locals() and 'track_generator' in locals()['kwargs']:
            track_generator = locals()['kwargs']['track_generator']

        # SWIG 2
        else:
            track_generator = locals()['args'][0]

        track_generator.thisown = False
%}


/* A Geometry owns the memory for its Cmfd (if any) */
%pythonappend Geometry::setCmfd %{
        # SWIG 3
        if 'cmfd' in locals():
            cmfd = locals()['cmfd']
        elif 'args' in locals() and 'cmfd' in locals()['args']:
            cmfd = locals()['args']['cmfd']
        elif 'kwargs' in locals() and 'cmfd' in locals()['kwargs']:
            cmfd = locals()['kwargs']['cmfd']

        # SWIG 2
        else:
            cmfd = locals()['args'][0]

        cmfd.thisown = False
%}
