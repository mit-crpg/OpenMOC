#include "TrackGenerator.h"


/**
 * @brief Constructor.
 * @param geometry a pointer to a geometry object
 * @param num_azim number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @param spacing track spacing (cm)
 */
TrackGenerator::TrackGenerator(Geometry* geometry, const int num_azim, 
                               const double spacing) {
    _geometry = geometry;
    setNumAzim(num_azim);
    setTrackSpacing(spacing);
    _tot_num_tracks = 0;
    _tot_num_segments = 0;
    _num_segments = NULL;
    _contains_tracks = false;
    _use_input_file = false;
    _tracks_filename = "";
}


/**
 * @brief Destructor frees memory for all tracks.
 */
TrackGenerator::~TrackGenerator() {

    /* Deletes tracks arrays if tracks have been generated */
    if (_contains_tracks) {
        delete [] _num_tracks;
        delete [] _num_segments;
        delete [] _num_x;
        delete [] _num_y;
        delete [] _azim_weights;

        for (int i = 0; i < _num_azim; i++)
            delete [] _tracks[i];

        delete [] _tracks;
    }
}


/**
 * @brief Return the number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @return the number of azimuthal angles in \f$ 2\pi \f$
 */
int TrackGenerator::getNumAzim() {
    return _num_azim * 2.0;
}


/**
 * @brief Return the user-specified (or default) suggested track spacing (cm).
 * @return the suggested track spacing (cm)
 */
double TrackGenerator::getTrackSpacing() {
    return _spacing;
}



/**
 * @brief Return the geometry for this trackgenerator if one has been set.
 * @return a pointer to the geometry
 */
Geometry* TrackGenerator::getGeometry() {
    if (_geometry == NULL)
        log_printf(ERROR, "Unable to return the TrackGenerator's geometry "
                 "since it has not yet been set");

    return _geometry;
}


/**
 * @brief Return the total number of tracks across the geometry.
 * @return the total number of tracks
 */
int TrackGenerator::getNumTracks() {

    if (!_contains_tracks)
        log_printf(ERROR, "Unable to return the total number of tracks since "
                   "tracks have not yet been generated.");

    return _tot_num_tracks;
}

/**
 * @brief Return an array of the number of tracks for each azimuthal angle.
 * @return array with the number of tracks
 */
int* TrackGenerator::getNumTracksArray() {
    if (!_contains_tracks)
        log_printf(ERROR, "Unable to return the array of the number of "
                   "tracks per azimuthal angle since tracks have not yet "
                   "been generated.");

    return _num_tracks;
}


/**
 * @brief Return the total number of track segments across the geometry.
 * @return the total number of track segments
 */
int TrackGenerator::getNumSegments() {

    if (!_contains_tracks)
        log_printf(ERROR, "Unable to return the total number of segments since "
                   "tracks have not yet been generated.");

    return _tot_num_segments;
}

/**
 * @brief Return an array of the number of segments per track.
 * @return array with the number of segments per track
 */
int* TrackGenerator::getNumSegmentsArray() {
    if (!_contains_tracks)
        log_printf(ERROR, "Unable to return the array of the number of "
                   "segments per track since tracks have not yet "
                   "been generated.");

    return _num_segments;
}


/**
 * @brief Returns a 2D jagged array of the tracks.
 * @details The first index into the array is the azimuthal angle and the
 *          second index is the track number for a given azimuthal angle.
 * @return the 2D jagged array of tracks
 */
Track **TrackGenerator::getTracks() {
    if (!_contains_tracks)
        log_printf(ERROR, "Unable to return the 2D ragged array of the tracks "
                   "since tracks have not yet been generated.");

    return _tracks;
}


/**
 * @brief Return a pointer to the array of azimuthal weights.
 * @return the array of azimuthal weights
 */
double* TrackGenerator::getAzimWeights() {
    if (!_contains_tracks)
        log_printf(ERROR, "Unable to return track azimuthal weights since "
                 "tracks have not yet been generated.");

    return _azim_weights;
}


/**
 * @brief Returns whether or not the trackgenerator contains track that are
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the trackgenerator conatains tracks; false otherwise
 */
bool TrackGenerator::containsTracks() {
    return _contains_tracks;
}


/**
 * @brief Fills an array with the x,y coordinates for each track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires on due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_tracks = track_generator.getNumTracks()
 *          coords = track_generator.retrieveTrackCoords(num_tracks*4)
 * @endcode
 *
 * @param coords an array of coords of length 4 times the number of tracks
 * @param num_tracks the total number of tracks
 */
void TrackGenerator::retrieveTrackCoords(double* coords, int num_tracks) {

    if (num_tracks != 4*getNumTracks())
        log_printf(ERROR, "Unable to retrieve the track coordinates since "
                 "the track generator contains %d tracks with %d coordinates "
                 "but an array of length %d was input", getNumTracks(),
                 4*getNumTracks(), num_tracks);

    int counter = 0;
    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_tracks[i]; j++) {
            coords[counter] = _tracks[i][j].getStart()->getX();
            coords[counter+1] = _tracks[i][j].getStart()->getY();
            coords[counter+2] = _tracks[i][j].getEnd()->getX();
            coords[counter+3] = _tracks[i][j].getEnd()->getY();
            counter += 4;
        }
    }

    return;
}


/**
 * @brief Fills an array with the x,y coordinates for each track segment.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          segments. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_segments = track_generator.getNumSegments()
 *          coords = track_generator.retrieveSegmentCoords(num_segments*4)
 * @endcode
 *
 * @param coords an array of coords of length 5 times the number of segments
 * @param num_segments the total number of track segments
 */
void TrackGenerator::retrieveSegmentCoords(double* coords, int num_segments) {

    if (num_segments != 5*getNumSegments())
        log_printf(ERROR, "Unable to retrieve the track segment coordinates "
                   "since the track generator contains %d segments with %d "
                   "coordinates but an array of length %d was input", 
                   getNumSegments(), 5*getNumSegments(), num_segments);

    segment* curr_segment = NULL;
    double x0, x1, y0, y1;
    double phi;
    segment* segments;    

    int counter = 0;

    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_tracks[i]; j++) {

            x0 = _tracks[i][j].getStart()->getX();
            y0 = _tracks[i][j].getStart()->getY();
            phi = _tracks[i][j].getPhi();

            segments = _tracks[i][j].getSegments();

	    for (int s=0; s < _tracks[i][j].getNumSegments(); s++) {
	        curr_segment = &segments[s];

                coords[counter] = curr_segment->_region_id;

                coords[counter+1] = x0;
                coords[counter+2] = y0;

                x1 = x0 + cos(phi) * curr_segment->_length;
                y1 = y0 + sin(phi) * curr_segment->_length;

                coords[counter+3] = x1;
                coords[counter+4] = y1;

                x0 = x1;
                y0 = y1;

                counter += 5;
            }
        }
    }

    return;
}


/**
 * @brief Set the number of azimuthal angles in \f$ [0, 2\pi] \f$.
 * @param num_azim the number of azimuthal angles in \f$ 2\pi \f$
 */
void TrackGenerator::setNumAzim(int num_azim) {
    if (num_azim < 0)
        log_printf(ERROR, "Unable to set a negative number of azimuthal angles "
                 "%d for the TrackGenerator.", num_azim);
    if (num_azim % 4 != 0)
        log_printf(ERROR, "Unable to set the number of azimuthal angles to %d "
                   "since it is not a multiple of 4", num_azim);

    _num_azim = num_azim / 2.0;    /* Subdivide out angles in [pi,2pi] */
    _contains_tracks = false;
    _use_input_file = false;
    _tracks_filename = "";
}


/**
 * @brief Set the suggested track spacing (cm).
 * @param spacing the suggested track spacing
 */
void TrackGenerator::setTrackSpacing(double spacing) {
    if (spacing < 0)
        log_printf(ERROR, "Unable to set a negative number of track spacing "
                 "%f for the TrackGenerator.", spacing);

    _spacing = spacing;
    _tot_num_tracks = 0;
    _tot_num_segments = 0;
    _contains_tracks = false;
    _use_input_file = false;
    _tracks_filename = "";
}


/**
 * @brief Set a pointer to the geometry to use for track generation.
 * @param geometry a pointer to the geometry
 */
void TrackGenerator::setGeometry(Geometry* geometry) {
    _geometry = geometry;
    _tot_num_tracks = 0;
    _tot_num_segments = 0;
    _contains_tracks = false;
    _use_input_file = false;
    _tracks_filename = "";
}

/**
 * @brief Generates tracks for some number of azimuthal angles and track spacing * @details Computes the effective angles and track spacings. Computes the 
 *          number of tracks for each azimuthal angle, allocates memory for 
 *          all tracks at each angle and sets each track's starting and ending 
 *          points, azimuthal weight, and azimuthal angle.
 */
void TrackGenerator::generateTracks() {

    if (_geometry == NULL)
        log_printf(NORMAL, "Unable to generate tracks since no geometry "
                   "has been set for the TrackGenerator");

    /* Deletes tracks arrays if tracks have been generated */
    if (_contains_tracks) {
        delete [] _num_tracks;
        delete [] _num_segments;
        delete [] _num_x;
        delete [] _num_y;
        delete [] _azim_weights;

        for (int i = 0; i < _num_azim; i++)
            delete [] _tracks[i];

        delete [] _tracks;
    }

    /** Create tracks directory if one does not yet exist */
    std::stringstream directory;
    directory << getOutputDirectory() << "/tracks";
    struct stat st;
    if (!stat(directory.str().c_str(), &st) == 0)
        mkdir(directory.str().c_str(), S_IRWXU);    

    struct stat buffer;
    std::stringstream test_filename;
    test_filename << directory.str() << "/tracks_" 
		  <<  _num_azim*2.0 << "_angles_" 
		  << _spacing << "_cm_spacing.data";

    _tracks_filename = test_filename.str();

    if (!stat(_tracks_filename.c_str(), &buffer)) {
        if (readTracksFromFile()) {
            _use_input_file = true;
            _contains_tracks = true;
	}
    }
    
    /* If not tracks input file exists, generate tracks */
    if (_use_input_file == false) {

        /* Allocate memory for the tracks */
        try {
            _num_tracks = new int[_num_azim];
            _num_x = new int[_num_azim];
            _num_y = new int[_num_azim];
            _azim_weights = new double[_num_azim];
            _tracks = new Track*[_num_azim];
        }
        catch (std::exception &e) {
            log_printf(ERROR, "Unable to allocate memory for TrackGenerator. "
                       "Backtrace:\n%s", e.what());
        }

        /* Check to make sure that height, width of the geometry are nonzero */
        if (_geometry->getHeight() <= 0 || _geometry->getHeight() <= 0)
            log_printf(ERROR, "The total height and width of the geometry "
                       "must be nonzero for track generation. Please specify "
                       "the height and width in the geometry input file.");

        try {
            log_printf(NORMAL, "Computing azimuthal angles and track "
                       "spacings...");

            /* Each element in arrays corresponds to a track angle in phi_eff */
            /* Track spacing along x,y-axes, and perpendicular to each track */ 
            double* dx_eff = new double[_num_azim];
            double* dy_eff = new double[_num_azim];
            double* d_eff = new double[_num_azim];

            /* Effective azimuthal angles with respect to positive x-axis */
            double* phi_eff = new double[_num_azim];

            double x1, x2;
            double iazim = _num_azim*2.0;
            double width = _geometry->getWidth();
            double height = _geometry->getHeight();

            /* Determine azimuthal angles and track spacing */
            for (int i = 0; i < _num_azim; i++) {

                /* desired angle */
                double phi = 2.0 * M_PI / iazim * (0.5 + i);

                /* num intersections with x,y-axes */
                _num_x[i] = (int) (fabs(width / _spacing * sin(phi))) + 1;
                _num_y[i] = (int) (fabs(height / _spacing * cos(phi))) + 1;

                /* total num of tracks */
                _num_tracks[i] = _num_x[i] + _num_y[i];
            
                /* effective/actual angle (not the angle we desire, but close */
                phi_eff[i] = atan((height * _num_x[i]) / (width * _num_y[i]));

                /* fix angles in range(pi/2, pi) */
                if (phi > M_PI / 2)
                    phi_eff[i] = M_PI - phi_eff[i];

                /* Effective track spacing (not spacing we desire, but close */
                dx_eff[i] = (width / _num_x[i]);
                dy_eff[i] = (height / _num_y[i]);
                d_eff[i] = (dx_eff[i] * sin(phi_eff[i]));
            }

            /* Compute azimuthal weights */
            for (int i = 0; i < _num_azim; i++) {
          
                if (i < _num_azim - 1)
                    x1 = 0.5 * (phi_eff[i+1] - phi_eff[i]);
                else
                    x1 = 2 * M_PI / 2.0 - phi_eff[i];

                if (i >= 1)
                    x2 = 0.5 * (phi_eff[i] - phi_eff[i-1]);
                else
                    x2 = phi_eff[i];

                /* Multiply weight by 2 because angles are in [0, Pi] */
                _azim_weights[i] = (x1 + x2) / (2 * M_PI) * d_eff[i] * 2;
            }

            log_printf(NORMAL, "Generating track start and end points...");

            /* Compute track starting and end points */
            for (int i = 0; i < _num_azim; i++) {

                /* Tracks for azimuthal angle i */
                _tracks[i] = new Track[_num_tracks[i]];

                /* Compute start points for tracks starting on x-axis */
                for (int j = 0; j < _num_x[i]; j++)
                    _tracks[i][j].getStart()->setCoords(dx_eff[i] * (0.5+j), 0);

                /* Compute start points for tracks starting on y-axis */
                for (int j = 0; j < _num_y[i]; j++) {

                    /* If track points to the upper right */
                    if (sin(phi_eff[i]) > 0 && cos(phi_eff[i]) > 0)
                        _tracks[i][_num_x[i]+j].getStart()->setCoords(0, 
                                                        dy_eff[i] * (0.5 + j));

                    /* If track points to the upper left */
                    else if (sin(phi_eff[i]) > 0 && cos(phi_eff[i]) < 0)
                        _tracks[i][_num_x[i]+j].getStart()->setCoords(width,
                                                       dy_eff[i] * (0.5 + j));
                }

                /* Compute the exit points for each track */
                for (int j = 0; j < _num_tracks[i]; j++) {

                    /* Set the track's end point */
                    Point* start = _tracks[i][j].getStart();
                    Point* end = _tracks[i][j].getEnd();
                    computeEndPoint(start, end, phi_eff[i], width, height);

                    /* Set the tracks azimuthal angle */
                    _tracks[i][j].setPhi(phi_eff[i]);
                }
            }

            /* Recalibrate track start / end points to geometry global origin */
            int uid = 0;

            for (int i = 0; i < _num_azim; i++) {
	        _tot_num_tracks += _num_tracks[i];

                for (int j = 0; j < _num_tracks[i]; j++) {

		    _tracks[i][j].setUid(uid);
                    uid++;

                    double x0 = _tracks[i][j].getStart()->getX();
                    double y0 = _tracks[i][j].getStart()->getY();
                    double x1 = _tracks[i][j].getEnd()->getX();
                    double y1 = _tracks[i][j].getEnd()->getY();
                    double new_x0 = x0 - _geometry->getWidth()/2.0;
                    double new_y0 = y0 - _geometry->getHeight()/2.0;
                    double new_x1 = x1 - _geometry->getWidth()/2.0;
                    double new_y1 = y1 - _geometry->getHeight()/2.0;
                    double phi = _tracks[i][j].getPhi();

                    _tracks[i][j].setValues(new_x0, new_y0, new_x1,new_y1, phi);
                    _tracks[i][j].setAzimAngleIndex(i);
                
                    if (i != _tracks[i][j].getAzimAngleIndex())
                        log_printf(NORMAL, "i = %d, but azim index = %d", i, 
                                   _tracks[i][j].getAzimAngleIndex());
                }
            }

            delete [] dx_eff;
            delete [] dy_eff;
            delete [] d_eff;
            delete [] phi_eff;

            segmentize();
            _contains_tracks = true;
            dumpTracksToFile();
            _use_input_file = true;
        }

        catch (std::exception &e) {
            log_printf(ERROR, "Unable to allocate memory needed to generate "
                     "tracks. Backtrace:\n%s", e.what());
        }
    }

    initializeBoundaryConditions();
    return;
}


/**
 * @brief This helper method for generateTracks() method finds the end point 
 *        of a track with a defined start point and an angle from x-axis. 
 * @details This function does not return a value but instead saves the x/y 
 *          coordinates of the end point directly within the track's end point.
 * @param start pointer to the track start point
 * @param end pointer to a point to store the end point coordinates
 * @param phi the azimuthal angle
 * @param width the width of the geometry (cm)
 * @param height the height of the geometry (cm)
 */
void TrackGenerator::computeEndPoint(Point* start, Point* end,
                                     const double phi, const double width, 
                                     const double height) {

    double m = sin(phi) / cos(phi);             /* slope */
    double yin = start->getY();                 /* y-coord */
    double xin = start->getX();                 /* x-coord */

    try {
        Point *points = new Point[4];

        /* Determine all possible points */
        points[0].setCoords(0, yin - m * xin);
        points[1].setCoords(width, yin + m * (width - xin));
        points[2].setCoords(xin - yin / m, 0);
        points[3].setCoords(xin - (yin - height) / m, height);

        /* For each of the possible intersection points */
        for (int i = 0; i < 4; i++) {
            /* neglect the trivial point (xin, yin) */
            if (points[i].getX() == xin && points[i].getY() == yin) { }

            /* The point to return will be within the bounds of the cell */
            else if (points[i].getX() >= 0 && points[i].getX() <= width
                     && points[i].getY() >= 0 && points[i].getY() <= height) {
                end->setCoords(points[i].getX(), points[i].getY());
            }
        }

        delete[] points;
        return;
    }

    catch (std::exception &e) {
        log_printf(ERROR, "Unable to allocate memory for intersection points "
                   "in computeEndPoint method of TrackGenerator. "
                   "Backtrace:\n%s", e.what());
    }
}


/**
 * @brief Initializes boundary conditions for each track. 
 * @details Sets boundary conditions by setting the incoming and outgoing tracks *          for each track using a special indexing scheme into the 2d jagged 
 *          array of tracks
 */
void TrackGenerator::initializeBoundaryConditions() {

    log_printf(NORMAL, "Initializing track boundary conditions...");

    int nxi, nyi, nti;   /* nx, ny, nt for a particular angle */
    Track *curr;
    Track *refl;

    /* Loop over only half the angles since we will set the pointers for
     * connecting tracks at the same time
     */
    for (int i = 0; i < floor(_num_azim / 2); i++) {
        nxi = _num_x[i];
        nyi = _num_y[i];
        nti = _num_tracks[i];
        curr = _tracks[i];
        refl = _tracks[_num_azim - i - 1];

        /* Loop over all of the tracks for this angle */
        for (int j = 0; j < nti; j++) {
          
            /* More tracks starting along x-axis than y-axis */
            if (nxi <= nyi) {
                /* Bottom to right hand side */
                if (j < nxi) {
                    curr[j].setTrackIn(&refl[j]);
                    curr[j].setTrackInI(_num_azim - i - 1);
                    curr[j].setTrackInJ(j);
                    
                    refl[j].setTrackIn(&curr[j]);
                    refl[j].setTrackInI(i);
                    refl[j].setTrackInJ(j);

                    curr[j].setReflIn(false);
                    refl[j].setReflIn(false);
                    
		    if (_geometry->getBCBottom() == REFLECTIVE) {
		        curr[j].setBCIn(1);
			refl[j].setBCIn(1);
		    }
		    else {
		        curr[j].setBCIn(0);
			refl[j].setBCIn(0);
		    }

                    curr[j].setTrackOut(&refl[2 * nxi - 1 - j]);
                    curr[j].setTrackOutI(_num_azim - i - 1);
                    curr[j].setTrackOutJ(2 * nxi - 1 - j);

                    refl[2 * nxi - 1 - j].setTrackIn(&curr[j]);
                    refl[2 * nxi - 1 - j].setTrackInI(i);
                    refl[2 * nxi - 1 - j].setTrackInJ(j);

                    curr[j].setReflOut(false);
                    refl[2 * nxi - 1 - j].setReflIn(true);

		    if (_geometry->getBCRight() == REFLECTIVE) {
		        curr[j].setBCOut(1);
			refl[2 * nxi - 1 - j].setBCIn(1);
		    }
		    else {
		        curr[j].setBCOut(0);
			refl[2 * nxi - 1 - j].setBCIn(0);
		    }
                }

                /* Left hand side to right hand side */
                else if (j < nyi) {
                    curr[j].setTrackIn(&refl[j - nxi]);
                    curr[j].setTrackInI(_num_azim - i - 1);
                    curr[j].setTrackInJ(j - nxi);

                    refl[j - nxi].setTrackOut(&curr[j]);
                    refl[j - nxi].setTrackOutI(i);
                    refl[j - nxi].setTrackOutJ(j);

                    curr[j].setReflIn(true);
                    refl[j - nxi].setReflOut(false);

		    if (_geometry->getBCLeft() == REFLECTIVE) {
		        curr[j].setBCIn(1);
			refl[j - nxi].setBCOut(1);
		    }
		    else {
		        curr[j].setBCIn(0);
			refl[j - nxi].setBCOut(0);
		    }
                    
                    curr[j].setTrackOut(&refl[j + nxi]);
                    curr[j].setTrackOutI(_num_azim - i - 1);
                    curr[j].setTrackOutJ(j + nxi);

                    refl[j + nxi].setTrackIn(&curr[j]);
                    refl[j + nxi].setTrackInI(i);
                    refl[j + nxi].setTrackInJ(j);

                    curr[j].setReflOut(false);
                    refl[j + nxi].setReflIn(true);

		    if (_geometry->getBCRight() == REFLECTIVE) {
		        curr[j].setBCOut(1);
		        refl[j + nxi].setBCIn(1);
		    }
		    else {
		        curr[j].setBCOut(0);
		        refl[j + nxi].setBCIn(0);
		    }
                }

                /* Left hand side to top (j > ny) */
                else {
                    curr[j].setTrackIn(&refl[j - nxi]);
                    curr[j].setTrackInI(_num_azim - i - 1);
                    curr[j].setTrackInJ(j - nxi);

                    refl[j - nxi].setTrackOut(&curr[j]);
                    refl[j - nxi].setTrackOutI(i);
                    refl[j - nxi].setTrackOutJ(j);

                    curr[j].setReflIn(true);
                    refl[j - nxi].setReflOut(false);

		    if (_geometry->getBCLeft() == REFLECTIVE) {
		        curr[j].setBCIn(1);
		        refl[j - nxi].setBCOut(1);
		    }
		    else {
		        curr[j].setBCIn(0);
		        refl[j - nxi].setBCOut(0);
		    }

                    curr[j].setTrackOut(&refl[2 * nti - nxi - j - 1]);
                    curr[j].setTrackOutI(_num_azim - i - 1);
                    curr[j].setTrackOutJ(2 * nti - nxi - j - 1);

                    refl[2 * nti - nxi - j - 1].setTrackOut(&curr[j]);
                    refl[2 * nti - nxi - j - 1].setTrackOutI(i);
                    refl[2 * nti - nxi - j - 1].setTrackOutJ(j);

                    curr[j].setReflOut(true);
                    refl[2 * nti - nxi - j - 1].setReflOut(true);

		    if (_geometry->getBCTop() == REFLECTIVE) {
		        curr[j].setBCOut(1);
			refl[2 * nti - nxi - j - 1].setBCOut(1);
		    }
		    else {
		        curr[j].setBCOut(0);
			refl[2 * nti - nxi - j - 1].setBCOut(0);
		    }
                }
            }

            /* More tracks starting on y-axis than on x-axis */
            else {
                /* Bottom to top */
                if (j < nxi - nyi) {
                    curr[j].setTrackIn(&refl[j]);
                    curr[j].setTrackInI(_num_azim - i - 1);
                    curr[j].setTrackInJ(j);

                    refl[j].setTrackIn(&curr[j]);
                    refl[j].setTrackInI(i);
                    refl[j].setTrackInJ(j);
                    
                    curr[j].setReflIn(false);
                    refl[j].setReflIn(false);

		    if (_geometry->getBCBottom() == REFLECTIVE) {
		        curr[j].setBCIn(1);
			refl[j].setBCIn(1);
		    }
		    else {
		        curr[j].setBCIn(0);
			refl[j].setBCIn(0);
		    }

                    curr[j].setTrackOut(&refl[nti - (nxi - nyi) + j]);
                    curr[j].setTrackOutI(_num_azim - i - 1);
                    curr[j].setTrackOutJ(nti - (nxi - nyi) + j);
                    
                    refl[nti - (nxi - nyi) + j].setTrackOut(&curr[j]);
                    refl[nti - (nxi - nyi) + j].setTrackOutI(i);
                    refl[nti - (nxi - nyi) + j].setTrackOutJ(j);

                    curr[j].setReflOut(true);
                    refl[nti - (nxi - nyi) + j].setReflOut(true);

		    if (_geometry->getBCTop() == REFLECTIVE) {
		        curr[j].setBCOut(1);
			refl[nti - (nxi - nyi) + j].setBCOut(1);
		    }
		    else {
		        curr[j].setBCOut(0);
			refl[nti - (nxi - nyi) + j].setBCOut(0);
		    }
                }

                /* Bottom to right hand side */
                else if (j < nxi) {
                    curr[j].setTrackIn(&refl[j]);
                    curr[j].setTrackInI(_num_azim - i - 1);
                    curr[j].setTrackInJ(j);

                    refl[j].setTrackIn(&curr[j]);
                    refl[j].setTrackInI(i);
                    refl[j].setTrackInJ(j);

                    curr[j].setReflIn(false);
                    refl[j].setReflIn(false);

		    if (_geometry->getBCBottom() == REFLECTIVE) {
		        curr[j].setBCIn(1);
		        refl[j].setBCIn(1);
		    }
		    else {
		        curr[j].setBCIn(0);
		        refl[j].setBCIn(0);
		    }

                    curr[j].setTrackOut(&refl[nxi + (nxi - j) - 1]);
                    curr[j].setTrackOutI(_num_azim - i - 1);
                    curr[j].setTrackOutJ(nxi + (nxi - j) - 1);
                                
                    refl[nxi + (nxi - j) - 1].setTrackIn(&curr[j]);
                    refl[nxi + (nxi - j) - 1].setTrackInI(i);
                    refl[nxi + (nxi - j) - 1].setTrackInJ(j);

                    curr[j].setReflOut(false);
                    refl[nxi + (nxi - j) - 1].setReflIn(true);

		    if (_geometry->getBCRight() == REFLECTIVE) {
		        curr[j].setBCOut(1);
			refl[nxi + (nxi - j) - 1].setBCIn(1);
		    }
		    else {
		        curr[j].setBCOut(0);
			refl[nxi + (nxi - j) - 1].setBCIn(0);
		    }
                }

                /* Left-hand side to top (j > nx) */
                else {
                    curr[j].setTrackIn(&refl[j - nxi]);
                    curr[j].setTrackInI(_num_azim - i - 1);
                    curr[j].setTrackInJ(j - nxi);

                    refl[j - nxi].setTrackOut(&curr[j]);
                    refl[j - nxi].setTrackOutI(i);
                    refl[j - nxi].setTrackOutJ(j);

                    curr[j].setReflIn(true);
                    refl[j - nxi].setReflOut(false);

		    if (_geometry->getBCLeft() == REFLECTIVE) {
		        curr[j].setBCIn(1);
			refl[j - nxi].setBCOut(1);
		    }
		    else {
		        curr[j].setBCIn(0);
			refl[j - nxi].setBCOut(0);
		    }

                    curr[j].setTrackOut(&refl[nyi + (nti - j) - 1]);
                    curr[j].setTrackOutI(_num_azim - i - 1);
                    curr[j].setTrackOutJ(nyi + (nti - j) - 1);

                    refl[nyi + (nti - j) - 1].setTrackOut(&curr[j]);
                    refl[nyi + (nti - j) - 1].setTrackOutI(i);
                    refl[nyi + (nti - j) - 1].setTrackOutJ(j);

                    curr[j].setReflOut(true);
                    refl[nyi + (nti - j) - 1].setReflOut(true);

		    if (_geometry->getBCTop() == REFLECTIVE) {
		        curr[j].setBCOut(1);
			refl[nyi + (nti - j) - 1].setBCOut(1);
		    }
		    else {
		        curr[j].setBCOut(0);
			refl[nyi + (nti - j) - 1].setBCOut(0);
		    }
                }
            }
        }
    }

    return;
}


/**
 * @brief Generate segments for each track across the geometry.
 */
void TrackGenerator::segmentize() {

    log_printf(NORMAL, "Segmenting tracks...");

    Track* track;

    if (_num_segments != NULL)
        delete [] _num_segments;

    /* This section loops over all track and segmentizes each one if the
     * tracks were not read in from an input file */
    if (!_use_input_file) {

        /* Loop over all tracks */
        #pragma omp parallel for private(track)
        for (int i=0; i < _num_azim; i++) {
            for (int j=0; j < _num_tracks[i]; j++){
                track = &_tracks[i][j];
 	        log_printf(DEBUG, "Segmenting track %d/%d with i = %d, j = %d",
			   track->getUid(), _tot_num_tracks, i, j);
                _geometry->segmentize(track);
            }
        }

	/* Compute the total number of segments in the simulation */
	_num_segments = new int[_tot_num_tracks];
	_tot_num_segments = 0;
	for (int i=0; i < _num_azim; i++) {
	    for (int j=0; j < _num_tracks[i]; j++) {
	        track = &_tracks[i][j];
		_num_segments[track->getUid()] = track->getNumSegments();
	        _tot_num_segments += _num_segments[track->getUid()];
	    }
	}
    }

    return;
}


/**
 * @brief Writes all track and segment data to a *.tracks file.
 * @details Storing tracks in a file saves time by eliminating segmentation
 *          for commonly simulated geometries.
 */
void TrackGenerator::dumpTracksToFile() {

    if (!_contains_tracks) 
        log_printf(ERROR, "Unable to dump tracks to a file since no tracks "
                 "have been generated for %d azimuthal angles and %f "
                 "track spacing", _num_azim, _spacing);

    FILE* out;
    out = fopen(_tracks_filename.c_str(), "w");

    std::string geometry_to_string = _geometry->toString();
    int string_length = geometry_to_string.length();

    fwrite(&string_length, sizeof(int), 1, out);
    fwrite(geometry_to_string.c_str(), sizeof(char)*string_length, 1, out);

    fwrite(&_num_azim, sizeof(int), 1, out);
    fwrite(&_spacing, sizeof(double), 1, out);
    fwrite(_num_tracks, sizeof(int), _num_azim, out);
    fwrite(_num_x, sizeof(int), _num_azim, out);
    fwrite(_num_y, sizeof(int), _num_azim, out);
    fwrite(_azim_weights, sizeof(double), _num_azim, out);
    
    Track* curr_track;
    double x0, y0, x1, y1;
    double phi;
    int azim_angle_index;
    int num_segments;
    std::vector<segment*> _segments;

    segment* curr_segment;
    double length;
    int material_id;
    int region_id;

    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_tracks[i]; j++) {
            curr_track = &_tracks[i][j];
            x0 = curr_track->getStart()->getX();
            y0 = curr_track->getStart()->getY();
            x1 = curr_track->getEnd()->getX();
            y1 = curr_track->getEnd()->getY();
            phi = curr_track->getPhi();
            azim_angle_index = curr_track->getAzimAngleIndex();
            num_segments = curr_track->getNumSegments();

            fwrite(&x0, sizeof(double), 1, out);
            fwrite(&y0, sizeof(double), 1, out);
            fwrite(&x1, sizeof(double), 1, out);
            fwrite(&y1, sizeof(double), 1, out);
            fwrite(&phi, sizeof(double), 1, out);
            fwrite(&azim_angle_index, sizeof(int), 1, out);
            fwrite(&num_segments, sizeof(int), 1, out);
            
            for (int s=0; s < num_segments; s++) {
                curr_segment = curr_track->getSegment(s);
                length = curr_segment->_length;
                material_id = curr_segment->_material->getId();
                region_id = curr_segment->_region_id;

                fwrite(&length, sizeof(double), 1, out);
                fwrite(&material_id, sizeof(int), 1, out);
                fwrite(&region_id, sizeof(int), 1, out);
            }
        }
    }

    fclose(out);
}


/**
 * @brief Reads tracks in from a file.
 * @details Storing tracks in a file saves time by eliminating segmentation
 *          for commonly simulated geometries.
 */
bool TrackGenerator::readTracksFromFile() {

    /* Deletes tracks arrays if tracks have been generated */
    if (_contains_tracks) {
        delete [] _num_tracks;
        delete [] _num_segments;
        delete [] _num_x;
        delete [] _num_y;
        delete [] _azim_weights;

        for (int i = 0; i < _num_azim; i++)
            delete [] _tracks[i];

        delete [] _tracks;
    }

    int ret;
    FILE* in;
    in = fopen(_tracks_filename.c_str(), "r");

    int string_length;

    ret = fread(&string_length, sizeof(int), 1, in);
    char* geometry_to_string = new char[string_length];
    ret = fread(geometry_to_string, sizeof(char)*string_length, 1, in);

    if (_geometry->toString().compare(geometry_to_string) != 0) {
        log_printf(NORMAL, "Returning from readTracksFromFile");
        return false;
    }

    log_printf(NORMAL, "Reading tracks in from file");

    ret = fread(&_num_azim, sizeof(int), 1, in);
    ret = fread(&_spacing, sizeof(double), 1, in);

    _num_tracks = new int[_num_azim];
    _num_x = new int[_num_azim];
    _num_y = new int[_num_azim];
    _azim_weights = new double[_num_azim];
    _tracks = new Track*[_num_azim];

    ret = fread(_num_tracks, sizeof(int), _num_azim, in);
    ret = fread(_num_x, sizeof(int), _num_azim, in);
    ret = fread(_num_y, sizeof(int), _num_azim, in);
    ret = fread(_azim_weights, sizeof(double), _num_azim, in);

    Track* curr_track;
    double x0, y0, x1, y1;
    double phi;
    int azim_angle_index;
    int num_segments;

    double length;
    int material_id;
    int region_id;

    for (int i=0; i < _num_azim; i++)
        _tot_num_tracks += _num_tracks[i];
    _num_segments = new int[_tot_num_tracks];
    int uid = 0;
    _tot_num_segments = 0;
    
    for (int i=0; i < _num_azim; i++) {

        _tracks[i] = new Track[_num_tracks[i]];

        for (int j=0; j < _num_tracks[i]; j++) {

            ret = fread(&x0, sizeof(double), 1, in);
            ret = fread(&y0, sizeof(double), 1, in);
            ret = fread(&x1, sizeof(double), 1, in);
            ret = fread(&y1, sizeof(double), 1, in);
            ret = fread(&phi, sizeof(double), 1, in);
            ret = fread(&azim_angle_index, sizeof(int), 1, in);
            ret = fread(&num_segments, sizeof(int), 1, in);

            _tot_num_segments += num_segments;
            _num_segments[uid] += num_segments;

            curr_track = &_tracks[i][j];
            curr_track->setValues(x0, y0, x1, y1, phi);
            curr_track->setUid(uid);
            curr_track->setAzimAngleIndex(azim_angle_index);
  
            for (int s=0; s < num_segments; s++) {
                ret = fread(&length, sizeof(double), 1, in);
                ret = fread(&material_id, sizeof(int), 1, in);
                ret = fread(&region_id, sizeof(int), 1, in);

                segment* curr_segment = new segment;
                curr_segment->_length = length;
                curr_segment->_material = _geometry->getMaterial(material_id);
                curr_segment->_region_id = region_id;
                curr_track->addSegment(curr_segment);
            }

            uid++;
        }
    }

    if (ret)
        _contains_tracks = true;
    fclose(in);    

    delete [] geometry_to_string;
    
    return true;
}
