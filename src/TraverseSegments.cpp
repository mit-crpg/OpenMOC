#include "TraverseSegments.h"

/**
 * @brief Constructor for the TraverseSegments class assigns the TrackGenerator
 *        and pulls relevant information from it.
 */
TraverseSegments::TraverseSegments(TrackGenerator* track_generator) {

  /* Save the track generator */
  _track_generator = track_generator;

  /* Determine the type of segment formation used */
  _segment_formation = track_generator->getSegmentFormation();

  /* Determine if a global z-mesh is used for 3D calculations */
  _track_generator_3D = dynamic_cast<TrackGenerator3D*>(track_generator);
  if (_track_generator_3D != NULL) {
    _track_generator_3D->retrieveGlobalZMesh(_global_z_mesh, _mesh_size);
  }
}


/**
 * @brief Destructor for TraverseSegments.
 */
TraverseSegments::~TraverseSegments() {
}


/**
 * @brief Loops over Tracks, applying the provided kernel to all segments and
 *        the functionality described in onTrack(...) to all Tracks.
 * @details The segment formation method imported from the TrackGenerator
 *          during construction is used to redirect to the appropriate looping
 *          scheme. If a kernel is provided (not NULL) then it is deleted at
 *          the end of the looping scheme.
 * @param kernel MOCKernel to apply to all segments
 */
void TraverseSegments::loopOverTracks(MOCKernel* kernel) {

  switch (_segment_formation) {
    case EXPLICIT_2D:
      loopOverTracks2D(kernel);
      break;
    case EXPLICIT_3D:
      loopOverTracksExplicit(kernel);
      break;
    case OTF_TRACKS:
      loopOverTracksByTrackOTF(kernel);
      break;
    case OTF_STACKS:
      loopOverTracksByStackOTF(kernel);
      break;
  }

  if (kernel != NULL)
    delete kernel;
}


/**
 * @brief Loops over all explicit 2D Tracks.
 * @details The onTrack(...) function is applied to all 2D Tracks and the
 *          specified kernel is applied to all segments. If NULL is provided
 *          for the kernel, only the onTrack(...) functionality is applied.
 * @param kernel The MOCKernel dictating the functionality to apply to
 *        segments
 */
void TraverseSegments::loopOverTracks2D(MOCKernel* kernel) {

  /* Loop over all parallel tracks for each azimuthal angle */
  Track** tracks_2D = _track_generator->get2DTracksArray();
  long num_tracks = _track_generator->getNum2DTracks();

#pragma omp for schedule(dynamic)
  for (long t=0; t < num_tracks; t++) {

    Track* track_2D = tracks_2D[t];
    segment* segments = track_2D->getSegments();

    /* Operate on segments if necessary */
    if (kernel != NULL) {
      kernel->newTrack(track_2D);
      traceSegmentsExplicit(track_2D, kernel);
    }

    /* Operate on the Track */
    onTrack(track_2D, segments);
  }
}


/**
 * @brief Loops over all explicit 3D Tracks.
 * @details The onTrack(...) function is applied to all 3D Tracks and the
 *          specified kernel is applied to all segments. If NULL is provided
 *          for the kernel, only the onTrack(...) functionality is applied.
 * @param kernel The MOCKernel dictating the functionality to apply to
 *        segments
 */
void TraverseSegments::loopOverTracksExplicit(MOCKernel* kernel) {

  Track3D**** tracks_3D = _track_generator_3D->get3DTracks();
  int num_azim = _track_generator_3D->getNumAzim();
  int num_polar = _track_generator_3D->getNumPolar();
  int*** tracks_per_stack = _track_generator_3D->getTracksPerStack();

  /* Loop over all tracks, parallelizing over parallel 2D tracks */
  for (int a=0; a < num_azim/2; a++) {
    int num_xy = _track_generator->getNumX(a) + _track_generator->getNumY(a);
#pragma omp for schedule(dynamic) collapse(2)
    for (int i=0; i < num_xy; i++) {

      /* Loop over polar angles */
      for (int p=0; p < num_polar; p++) {

        /* Loop over tracks in the z-stack */
        for (int z=0; z < tracks_per_stack[a][i][p]; z++) {

          /* Extract 3D track */
          Track* track_3D = &tracks_3D[a][i][p][z];

          /* Operate on segments if necessary */
          if (kernel != NULL) {

            /* Reset kernel for a new Track */
            kernel->newTrack(track_3D);

            /* Trace the segments on the track */
            traceSegmentsExplicit(track_3D, kernel);
          }

          /* Operate on the Track */
          segment* segments = track_3D->getSegments();
          onTrack(track_3D, segments);
        }
      }
    }
  }
}


/**
 * @brief Loops over all 3D Tracks using axial on-the-fly ray tracing by Track.
 * @details The onTrack(...) function is applied to all 3D Tracks and the
 *          specified kernel is applied to all segments. If NULL is provided
 *          for the kernel, only the onTrack(...) functionality is applied.
 * @param kernel The MOCKernel dictating the functionality to apply to
 *        segments
 */
void TraverseSegments::loopOverTracksByTrackOTF(MOCKernel* kernel) {

  int num_2D_tracks = _track_generator_3D->getNum2DTracks();
  Track** tracks_2D = _track_generator_3D->get2DTracksArray();
  int*** tracks_per_stack = _track_generator_3D->getTracksPerStack();
  int num_azim = _track_generator->getNumAzim();
  int num_polar = _track_generator_3D->getNumPolar();
  int tid = omp_get_thread_num();

  /* Loop over flattened 2D tracks */
#pragma omp for schedule(dynamic)
  for (int ext_id=0; ext_id < num_2D_tracks; ext_id++) {

    /* Extract indices of 3D tracks associated with the flattened track */
    Track* flattened_track = tracks_2D[ext_id];
    TrackStackIndexes tsi;
    tsi._azim = flattened_track->getAzimIndex();
    tsi._xy = flattened_track->getXYIndex();

    /* Loop over polar angles */
    for (int p=0; p < num_polar; p++) {

      /* Loop over tracks in the z-stack */
      for (int z=0; z < tracks_per_stack[tsi._azim][tsi._xy][p]; z++) {

        /* Extract 3D track and retrieve its information */
        Track3D track_3D;
        tsi._polar = p;
        tsi._z = z;
        _track_generator_3D->getTrackOTF(&track_3D, &tsi);

        /* Operate on segments if necessary */
        if (kernel != NULL) {

          /* Reset kernel for a new Track */
          kernel->newTrack(&track_3D);
          double theta = track_3D.getTheta();
          Point* start = track_3D.getStart();

          /* Trace the segments on the track */
          traceSegmentsOTF(flattened_track, start, theta, kernel);
          track_3D.setNumSegments(kernel->getCount());
        }

        /* Operate on the Track */
        segment* segments = _track_generator_3D->getTemporarySegments(tid);
        onTrack(&track_3D, segments);
      }
    }
  }
}


/**
 * @brief Loops over all 3D Tracks using axial on-the-fly ray tracing by
 *        z-stack.
 * @details The onTrack(...) function is applied to all 3D Tracks and the
 *          specified kernel is applied to all segments. If NULL is provided
 *          for the kernel, only the onTrack(...) functionality is applied.
 * @param kernel The MOCKernel dictating the functionality to apply to
 *        segments
 */
void TraverseSegments::loopOverTracksByStackOTF(MOCKernel* kernel) {

  int num_2D_tracks = _track_generator_3D->getNum2DTracks();
  Track** flattened_tracks = _track_generator_3D->get2DTracksArray();
  int*** tracks_per_stack = _track_generator_3D->getTracksPerStack();
  int num_polar = _track_generator_3D->getNumPolar();
  int tid = omp_get_thread_num();

  /* Allocate array of current Tracks */
  Track3D* current_stack = _track_generator_3D->getTemporary3DTracks(tid);

  /* Loop over flattened 2D tracks */
#pragma omp for schedule(dynamic)
  for (int ext_id=0; ext_id < num_2D_tracks; ext_id++) {

    /* Extract indices of 3D tracks associated with the flattened track */
    TrackStackIndexes tsi;
    Track* flattened_track = flattened_tracks[ext_id];
    tsi._azim = flattened_track->getAzimIndex();
    tsi._xy = flattened_track->getXYIndex();

    /* Loop over polar angles */
    for (int p=0; p < num_polar; p++) {

      /* Retrieve information for the first 3D Track in the z-stack */
      tsi._polar = p;
      int stack_size = tracks_per_stack[tsi._azim][tsi._xy][tsi._polar];
      for (int z=0; z < stack_size; z++) {
        tsi._z = z;
        _track_generator_3D->getTrackOTF(&current_stack[z], &tsi);
      }

      if (kernel != NULL) {

        /* Reset kernel to for the new base Track */
        kernel->newTrack(&current_stack[0]);

        /* Trace all segments in the z-stack */
        traceStackOTF(flattened_track, p, kernel);
        current_stack[0].setNumSegments(kernel->getCount());
      }

      /* Operate on the Track */
      segment* segments = _track_generator_3D->getTemporarySegments(tid);
      onTrack(&current_stack[0], segments);
    }
  }
}


/**
 * @brief Loops over segments in a Track when segments are explicitly generated.
 * @details All segments in the provided Track are looped over and the provided
 *          MOCKernel is applied to them.
 * @param track The Track whose segments will be traversed
 * @param kernel The kernel to apply to all segments
 */
void TraverseSegments::traceSegmentsExplicit(Track* track, MOCKernel* kernel) {

  /* Get direction of the track */
  double phi = track->getPhi();
  double theta = M_PI_2;
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  if (track_3D != NULL)
    theta = track_3D->getTheta();

  for (int s=0; s < track->getNumSegments(); s++) {
    segment* seg = track->getSegment(s);
    kernel->execute(seg->_length, seg->_material, seg->_region_id, 0,
                    seg->_cmfd_surface_fwd, seg->_cmfd_surface_bwd,
                    seg->_starting_position[0], seg->_starting_position[1],
                    seg->_starting_position[2], phi, theta);
  }
}


/**
 * @brief Computes 3D segment lengths on-the-fly for a single 3D track given an
 *        associated 2D Track with a starting point and a polar angle. The
 *        computed segments are passed to the provided kernel.
 * @details Segment lengths are computed on-the-fly using 2D segment lengths
 *          stored in a 2D Track object and 1D meshes from the extruded
 *          FSRs. Note: before calling this function with a SegmentationKernel,
 *          the memory for the segments should be allocated and referenced by
 *          the kernel using the setSegments routine in the kernel.
 * @param flattened_track the 2D track associated with the 3D track for which
 *        3D segments are computed
 * @param start the starting coordinates of the 3D track
 * @param theta the polar angle of the 3D track
 * @param kernel An MOCKernel object to apply to the calculated 3D segments
 */
void TraverseSegments::traceSegmentsOTF(Track* flattened_track, Point* start,
                                        double theta, MOCKernel* kernel) {

  /* Create unit vector */
  double phi = flattened_track->getPhi();
  double cos_phi = cos(phi);
  double sin_phi = sin(phi);
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  int sign = (cos_theta > 0) - (cos_theta < 0);

  /* Extract starting coordinates */
  double x_start_3D = start->getX();
  double x_start_2D = flattened_track->getStart()->getX();
  double x_coord = x_start_3D;
  double y_coord = start->getY();
  double z_coord = start->getZ();

  /* Find 2D distance from 2D edge to start of track */
  double start_dist_2D = (x_start_3D - x_start_2D) / cos_phi;

  /* Find starting 2D segment */
  int seg_start = 0;
  segment* segments_2D = flattened_track->getSegments();
  for (int s=0; s < flattened_track->getNumSegments(); s++) {

    /* Determine if start point of track is beyond current 2D segment */
    double seg_len_2D = segments_2D[s]._length;
    if (start_dist_2D > seg_len_2D) {
      start_dist_2D -= seg_len_2D;
      seg_start++;
    }
    else {
      break;
    }
  }

  Geometry* geometry = _track_generator_3D->getGeometry();
  Cmfd* cmfd = geometry->getCmfd();

  /* For very short tracks, it's possible no significant segments will be
   * traversed */
  if (seg_start == flattened_track->getNumSegments()) {
    log_printf(WARNING, "Track of zero length encountered at starting point "
               "%s traveling on 2D Track: %s at polar angle cos %3.2f "
               "degrees on domain with z-bounds %3.2f and %3.2f",
               start->toString().c_str(), flattened_track->toString().c_str(),
               cos(theta), geometry->getMinZ(), geometry->getMaxZ());
    return;
  }

  /* Extract the appropriate starting mesh */
  int num_fsrs;
  double* axial_mesh;
  bool contains_global_z_mesh;
  if (_global_z_mesh != NULL) {
    contains_global_z_mesh = true;
    num_fsrs = _mesh_size;
    axial_mesh = _global_z_mesh;
  }
  else {
    contains_global_z_mesh = false;
    int extruded_fsr_id = segments_2D[seg_start]._region_id;
    ExtrudedFSR* extruded_FSR = geometry->getExtrudedFSR(extruded_fsr_id);
    num_fsrs = extruded_FSR->_num_fsrs;
    axial_mesh = extruded_FSR->_mesh;
  }

  /* Get the starting z index */
  int z_ind = findMeshIndex(axial_mesh, num_fsrs+1, z_coord, sign);

  /* Loop over 2D segments */
  bool first_segment = true;
  bool segments_complete = false;
  for (int s=seg_start; s < flattened_track->getNumSegments(); s++) {

    /* Extract extruded FSR */
    int extruded_fsr_id = segments_2D[s]._region_id;
    ExtrudedFSR* extruded_FSR = geometry->getExtrudedFSR(extruded_fsr_id);

    /* Determine new mesh and z index */
    if (first_segment || contains_global_z_mesh) {
      first_segment = false;
    }
    else {
      /* Determine the axial region */
      num_fsrs = extruded_FSR->_num_fsrs;
      axial_mesh = extruded_FSR->_mesh;
      z_ind = findMeshIndex(axial_mesh, num_fsrs+1, z_coord, sign);
    }

    /* Extract 2D segment length */
    double remaining_length_2D = segments_2D[s]._length - start_dist_2D;
    start_dist_2D = 0;

    /* Transport along the 2D segment until it is completed */
    while (remaining_length_2D > 0) {

      /* Calculate 3D distance to z intersection */
      double z_dist_3D;
      if (sign > 0)
        z_dist_3D = (axial_mesh[z_ind+1] - z_coord) / cos_theta;
      else
        z_dist_3D = (axial_mesh[z_ind] - z_coord) / cos_theta;

      /* Calculate 3D distance to end of segment */
      double seg_dist_3D = remaining_length_2D / sin_theta;

      /* Calcualte shortest distance to intersection */
      double dist_2D;
      double dist_3D;
      int z_move;
      if (z_dist_3D <= seg_dist_3D) {
        dist_2D = z_dist_3D * sin_theta;
        dist_3D = z_dist_3D;
        z_move = sign;
      }
      else {
        dist_2D = remaining_length_2D;
        dist_3D = seg_dist_3D;
        z_move = 0;
      }

      /* Get the 3D FSR */
      long fsr_id = extruded_FSR->_fsr_ids[z_ind];

      /* Calculate CMFD surface */
      int cmfd_surface_bwd = -1;
      int cmfd_surface_fwd = -1;
      if (cmfd != NULL && dist_3D > TINY_MOVE) {

        /* Determine if this is the first 3D segment handled for the flattened
           2D segment. If so, get the 2D cmfd surface. */
        if (segments_2D[s]._length - remaining_length_2D <= TINY_MOVE)
          cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;

        /* Determine if this is the last 3D segment handled for the flattened
           2D segment. If so, get the 2D cmfd surface. */
        double next_dist_3D = (remaining_length_2D - dist_2D) / sin_theta;
        if (z_move == 0 || next_dist_3D <= TINY_MOVE)
          cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;

        /* Get CMFD cell */
        int cmfd_cell = geometry->getCmfdCell(fsr_id);

        /* Find the backwards surface */
        cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_coord,
                                                    cmfd_surface_bwd);

        /* Find forward surface */
        double z_coord_end = z_coord + dist_3D * cos_theta;
        cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_coord_end,
                                                    cmfd_surface_fwd);
      }

      /* Operate on segment */
      if (dist_3D > TINY_MOVE) {
        double x_centroid = 0;
        double y_centroid = 0;
        double z_centroid = 0;
        if (geometry->containsFSRCentroids()) {
          Point* centroid = geometry->getFSRCentroid(fsr_id);
          x_centroid = centroid->getX();
          y_centroid = centroid->getY();
          z_centroid = centroid->getZ();
        }

        kernel->execute(dist_3D, extruded_FSR->_materials[z_ind], fsr_id, 0,
                        cmfd_surface_fwd, cmfd_surface_bwd,
                        x_coord - x_centroid, y_coord - y_centroid,
                        z_coord - z_centroid, phi, theta);
      }

      /* Move axial height to end of segment */
      x_coord += dist_3D * sin_theta * cos_phi;
      y_coord += dist_3D * sin_theta * sin_phi;
      z_coord += dist_3D * cos_theta;

      /* Shorten remaining 2D segment length and move axial level */
      remaining_length_2D -= dist_2D;
      z_ind += z_move;

      /* Check if the track has crossed a Z boundary */
      if (z_ind < 0 or z_ind >= num_fsrs) {

        /* Reset z index */
        if (z_ind < 0)
          z_ind = 0;
        else
          z_ind = num_fsrs - 1;

        /* Mark the 2D segment as complete */
        segments_complete = true;
        break;
      }
    }

    /* Check if the track is completed due to an axial boundary */
    if (segments_complete)
      break;
  }
}


/**
 * @brief Computes 3D segment lengths on-the-fly for all tracks in a z-stack
 *        for a given associated 2D Track and a polar index on-the-fly and
 *        passes the computed segments to the provided kernel.
 * @details Segment lengths are computed on-the-fly using 2D segment lengths
 *          stored in a 2D Track object and 1D meshes from the extruded
 *          FSRs. Note: before calling this function with SegmentationKernels,
 *          the memory for the segments should be allocated and referenced by
 *          the kernel using the setSegments routine.
 * @param flattened_track the 2D track associated with the z-stack for which
 *        3D segments are computed
 * @param polar_index the index into the polar angles which is associated with
 *        the polar angle of the z-stack
 * @param kernel The MOCKernel to apply to the calculated 3D segments
 */
void TraverseSegments::traceStackOTF(Track* flattened_track, int polar_index,
                                     MOCKernel* kernel) {

  /* Extract information about the z-stack */
  int azim_index = flattened_track->getAzimIndex();
  int track_index = flattened_track->getXYIndex();
  int*** tracks_per_stack = _track_generator_3D->getTracksPerStack();
  int num_z_stack = tracks_per_stack[azim_index][track_index][polar_index];
  double z_spacing = _track_generator_3D->getZSpacing(azim_index, polar_index);

  /* Get information for the first Track in the z-stack */
  TrackStackIndexes tsi;
  Track3D first;
  tsi._azim = azim_index;
  tsi._xy = track_index;
  tsi._polar = polar_index;
  tsi._z = 0;
  _track_generator_3D->getTrackOTF(&first, &tsi);
  double theta = first.getTheta();

  /* Create unit vector */
  double phi = flattened_track->getPhi();
  double cos_phi = cos(phi);
  double sin_phi = sin(phi);
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  double tan_theta = sin_theta / cos_theta;
  int sign = (cos_theta > 0) - (cos_theta < 0);
  double track_spacing_3D = z_spacing / std::abs(cos_theta);

  /* Find 2D distance from 2D edge to start of track */
  double x_start_3D = first.getStart()->getX();
  double x_start_2D = flattened_track->getStart()->getX();
  double y_start_2D = flattened_track->getStart()->getY();
  double start_dist_2D = (x_start_3D - x_start_2D) / cos_phi;

  /* Calculate starting intersection of lowest track with z-axis */
  double z0 = first.getStart()->getZ();
  double start_z = z0 - start_dist_2D / tan_theta;

  /* Adapt for traceStackTwoWay reverse direction */
  //NOTE If more applications for this arise, make 'reverse' an argument
  if (dynamic_cast<TransportKernel*>(kernel) &&
      !dynamic_cast<TransportKernel*>(kernel)->getDirection()) {
    phi += M_PI;
    cos_phi *= -1;
    sin_phi *= -1;

    theta = M_PI - theta;
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    tan_theta = sin_theta / cos_theta;
    sign = (cos_theta > 0) - (cos_theta < 0);

    x_start_3D = first.getEnd()->getX();
    x_start_2D = flattened_track->getEnd()->getX();
    y_start_2D = flattened_track->getEnd()->getY();
    start_dist_2D = (x_start_3D - x_start_2D) / cos_phi;

    z0 = first.getEnd()->getZ();
    start_z = z0 - start_dist_2D / tan_theta;
  }


  /* Get the Geometry and CMFD mesh */
  Geometry* geometry = _track_generator_3D->getGeometry();
  Cmfd* cmfd = geometry->getCmfd();

  /* Extract the appropriate starting mesh */
  int num_fsrs;
  double* axial_mesh;
  if (_global_z_mesh != NULL) {
    num_fsrs = _mesh_size;
    axial_mesh = _global_z_mesh;
  }

  /* Set the current x and y coordinates */
  double x_curr = x_start_2D;
  double y_curr = y_start_2D;

  /* Loop over 2D segments */
  double first_start_z = start_z;
  segment* segments_2D = flattened_track->getSegments();
  for (int s=0; s < flattened_track->getNumSegments(); s++) {

    /* Get segment length and extruded FSR */
    double seg_length_2D = segments_2D[s]._length;
    int extruded_fsr_id = segments_2D[s]._region_id;
    ExtrudedFSR* extruded_FSR = geometry->getExtrudedFSR(extruded_fsr_id);

    /* Determine new mesh and z index */
    if (_global_z_mesh == NULL) {
      num_fsrs = extruded_FSR->_num_fsrs;
      axial_mesh = extruded_FSR->_mesh;
    }

    /* Calculate the end z coordinate of the first track */
    double first_end_z = first_start_z + seg_length_2D / tan_theta;

    /* Find the upper and lower z coordinates of the first track */
    double first_track_lower_z;
    double first_track_upper_z;
    if (sign > 0) {
      first_track_lower_z = first_start_z;
      first_track_upper_z = first_end_z;
    }
    else {
      first_track_lower_z = first_end_z;
      first_track_upper_z = first_start_z;
    }

    /* Loop over all 3D FSRs in the Extruded FSR to find intersections */
    double first_seg_len_3D;
    for (int z_iter = 0; z_iter < num_fsrs; z_iter++) {

      /* If traveling in negative-z direction, loop through FSRs from top */
      int z_ind = z_iter;
      if (sign < 0)
        z_ind = num_fsrs - z_iter - 1;

      /* Extract the FSR ID and Material ID of this 3D FSR */
      long fsr_id = extruded_FSR->_fsr_ids[z_ind];
      Material* material = extruded_FSR->_materials[z_ind];

      /* Find CMFD cell if necessary */
      int cmfd_cell;
      if (cmfd != NULL)
        cmfd_cell = geometry->getCmfdCell(fsr_id);

      /* Get boundaries of the current mesh cell */
      double z_min = axial_mesh[z_ind];
      double z_max = axial_mesh[z_ind+1];

      /* Calculate the local x and y centroid of the Extruded FSR */
      double fsr_x_start = 0;
      double fsr_y_start = 0;
      double z_cent = 0.0;
      if (geometry->containsFSRCentroids()) {
        Point* centroid = geometry->getFSRCentroid(fsr_id);
        fsr_x_start = x_curr - centroid->getX();
        fsr_y_start = y_curr - centroid->getY();
        z_cent = geometry->getFSRCentroid(fsr_id)->getZ();
      }

      /* Calculate z-stack track indexes that cross the 3D FSR */
      int start_track = std::ceil((z_min - first_track_upper_z) / z_spacing);
      int start_full = std::ceil((z_min - first_track_lower_z) / z_spacing);
      int end_full = std::ceil((z_max - first_track_upper_z) / z_spacing);
      int end_track = std::ceil((z_max - first_track_lower_z) / z_spacing);

      /* Check track bounds */
      start_track = std::max(start_track, 0);
      end_track = std::min(end_track, num_z_stack);

      /* Treat lower tracks that do not cross the entire 2D length */
      int min_lower = std::min(start_full, end_full);
      first_seg_len_3D = (first_track_upper_z - z_min) / std::abs(cos_theta);
      for (int i = start_track; i < min_lower; i++) {

        /* Calculate distance traveled in 3D FSR */
        double seg_len_3D = first_seg_len_3D + i * track_spacing_3D;

        /* Determine if segment length is large enough to operate on */
        if (seg_len_3D > TINY_MOVE) {

          /* Initialize CMFD surfaces to none (-1) */
          int cmfd_surface_fwd = -1;
          int cmfd_surface_bwd = -1;

          /* Get CMFD surface if necessary */
          double lower_z = first_track_lower_z + i * z_spacing;
          double upper_z = first_track_upper_z + i * z_spacing;
          double dist_to_corner = std::abs((z_min - lower_z) / cos_theta);
          if (cmfd != NULL) {
            if (sign > 0) {
              cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;
              cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, upper_z,
                                                          cmfd_surface_fwd);
              if (dist_to_corner <= TINY_MOVE)
                cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;
              cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_min,
                                                          cmfd_surface_bwd);
            }
            else {
              if (dist_to_corner <= TINY_MOVE)
                cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;
              cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_min,
                                                          cmfd_surface_fwd);
              cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;
              cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, upper_z,
                                                          cmfd_surface_bwd);
            }
          }

          /* Calculate the entry point of the segment into the FSR */
          double x_entry = fsr_x_start;
          double y_entry = fsr_y_start;
          double z_entry = 0;
          if (sign > 0) {
            double partial_2D = dist_to_corner * sin_theta;
            x_entry += partial_2D * cos_phi;
            y_entry += partial_2D * sin_phi;
            z_entry = z_min - z_cent;
          }
          else {
            z_entry = upper_z - z_cent;
          }

          /* Operate on segment */
          kernel->execute(seg_len_3D, material, fsr_id, i,
                          cmfd_surface_fwd, cmfd_surface_bwd,
                          x_entry, y_entry, z_entry, phi, theta);
        }
      }

      /* Find if there are tracks that traverse the entire 2D length */
      if (end_full > start_full) {

        /* Calculate distance traveled in 3D FSR */
        double seg_len_3D = seg_length_2D / sin_theta;

        /* Determine if segment length is large enough to operate on */
        if (seg_len_3D > TINY_MOVE) {

          /* Treat tracks that do cross the entire 2D length */
          for (int i = start_full; i < end_full; i++) {

            /* Initialize CMFD surfaces to 2D CMFD surfaces */
            int cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;
            int cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;

            /* Get CMFD surfaces if necessary */
            double start_z = first_start_z + i * z_spacing;
            if (cmfd != NULL) {

              /* Calculate start and end z */
              double end_z = first_end_z + i * z_spacing;
              cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, end_z,
                                                          cmfd_surface_fwd);
              cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, start_z,
                                                          cmfd_surface_bwd);
            }

            /* Calculate the entry point of the segment into the FSR */
            double x_entry = fsr_x_start;
            double y_entry = fsr_y_start;
            double z_entry = start_z - z_cent;

            /* Operate on segment */
            kernel->execute(seg_len_3D, material, fsr_id, i,
                            cmfd_surface_fwd, cmfd_surface_bwd,
                            x_entry, y_entry, z_entry, phi, theta);
          }
        }
      }

      /* Find if there are tracks that cross both upper and lower boundaries
         NOTE: this will only be true if there are no tracks that cross the
         entire 2D length in the FSR */
      else if (start_full > end_full) {

        /* Calculate distance traveled in 3D FSR */
        double seg_len_3D = (z_max - z_min) / std::abs(cos_theta);

        /* Determine if segment length is large enough to operate on */
        if (seg_len_3D > TINY_MOVE) {

          /* Treat tracks that cross through both the upper and lower axial
             boundaries */
          for (int i = end_full; i < start_full; i++) {

            /* Initialize CMFD surfaces to none (-1) */
            int cmfd_surface_bwd = -1;
            int cmfd_surface_fwd = -1;

            /* Determine start and end z */
            double enter_z;
            double exit_z;
            if (sign > 0) {
              enter_z = z_min;
              exit_z = z_max;
            }
            else {
              enter_z = z_max;
              exit_z = z_min;
            }

            /* Get CMFD surfaces if necessary */
            double track_start_z = first_start_z + i * z_spacing;
            double dist_to_corner_bwd = (enter_z - track_start_z) / cos_theta;
            if (cmfd != NULL) {

              /* Determine if any corners in the s-z plane are hit */
              if (dist_to_corner_bwd <= TINY_MOVE)
                cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;

              double track_end_z = first_end_z + i * z_spacing;
              double dist_to_corner_fwd = (track_end_z - exit_z) / cos_theta;
              if (dist_to_corner_fwd <= TINY_MOVE)
                cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;

              /* Find CMFD surfaces */
              cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, exit_z,
                                                          cmfd_surface_fwd);
              cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, enter_z,
                                                          cmfd_surface_bwd);
            }

            /* Calculate the entry point of the segment into the FSR */
            double partial_2D = dist_to_corner_bwd * sin_theta;
            double x_entry = fsr_x_start + partial_2D * cos_phi;
            double y_entry = fsr_y_start + partial_2D * sin_phi;
            double z_entry = enter_z - z_cent;

            /* Operate on segment */
            kernel->execute(seg_len_3D, material, fsr_id, i,
                            cmfd_surface_fwd, cmfd_surface_bwd,
                            x_entry, y_entry, z_entry, phi, theta);
          }
        }
      }

      /* Treat upper tracks that do not cross the entire 2D length */
      int min_upper = std::max(start_full, end_full);
      first_seg_len_3D = (z_max - first_track_lower_z) / std::abs(cos_theta);
      for (int i = min_upper; i < end_track; i++) {

        /* Calculate distance traveled in 3D FSR */
        double seg_len_3D = first_seg_len_3D - i * track_spacing_3D;

        /* Determine if segment length is large enough to operate on */
        if (seg_len_3D > TINY_MOVE) {

          /* Initialize CMFD surfaces to none (-1) */
          int cmfd_surface_fwd = -1;
          int cmfd_surface_bwd = -1;

          /* Get CMFD surface if necessary */
          double lower_z = first_track_lower_z + i * z_spacing;
          double upper_z = first_track_upper_z + i * z_spacing;
          double dist_to_corner = (upper_z - z_max) / std::abs(cos_theta);
          if (cmfd != NULL) {
            if (sign > 0) {
              if (dist_to_corner <= TINY_MOVE)
                cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;
              cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_max,
                                                          cmfd_surface_fwd);
              cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;
              cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, lower_z,
                                                          cmfd_surface_bwd);
            }
            else {
              cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;
              cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, lower_z,
                                                          cmfd_surface_fwd);
              if (dist_to_corner <= TINY_MOVE)
                cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;
              cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_max,
                                                          cmfd_surface_bwd);
            }
          }

          /* Calculate the entry point of the segment into the FSR */
          double x_entry = fsr_x_start;
          double y_entry = fsr_y_start;
          double z_entry = 0;
          if (sign < 0) {
            double partial_2D = dist_to_corner * sin_theta;
            x_entry += partial_2D * cos_phi;
            y_entry += partial_2D * sin_phi;
            z_entry = z_max - z_cent;
          }
          else {
            z_entry = lower_z - z_cent;
          }

          /* Operate on segment */
          kernel->execute(seg_len_3D, material, fsr_id, i,
                          cmfd_surface_fwd, cmfd_surface_bwd,
                          x_entry, y_entry, z_entry, phi, theta);
        }
      }
    }
    /* Traverse segment on first track */
    first_start_z = first_end_z;
    x_curr += seg_length_2D * cos_phi;
    y_curr += seg_length_2D * sin_phi;
  }
}


/**
 * @brief A function that searches for the index into a values mesh using a
 *        binary search.
 * @details A binary search is used to calculate the index into a mesh of where
 *          the value val resides. If a mesh boundary is hit, the upper region
 *          is selected for positive-z traversing rays and the lower region is
 *          selected for negative-z traversing rays.
 * @param values an array of monotonically increasing values
 * @param size the size of the values array
 * @param val the level to be searched for in the mesh
 * @param sign the direction of the ray in the z-direction
 */
int TraverseSegments::findMeshIndex(double* values, int size,
                                    double val, int sign) {

  /* Initialize indexes into the values array */
  int imin = 0;
  int imax = size-1;

  /* Check if val is outside the range */
  if (val < values[imin] or val > values[imax]) {
    log_printf(ERROR, "Value out of the mesh range in binary search");
    return -1;
  }

  /* Search for interval containing val */
  while (imax - imin > 1) {

    int imid = (imin + imax) / 2;

    if (val > values[imid])
      imin = imid;
    else if (val < values[imid])
      imax = imid;
    else {
      if (sign > 0)
        return imid;
      else
        return imid-1;
    }
  }
  return imin;
}


/**
 * @brief Loops over all 3D Tracks using axial on-the-fly ray tracking by
 *        z-stack, going forward then backward on each 3D Track.
 * @details The onTrack(...) function is applied to all 3D Tracks and the
 *          specified kernel is applied to all segments. If NULL is provided
 *          for the kernel, only the onTrack(...) functionality is applied.
 * @param kernel The TransportKernel dictating the functionality to apply to
 *        segments
 */
void TraverseSegments::loopOverTracksByStackTwoWay(TransportKernel* kernel) {

  if (_segment_formation != OTF_STACKS)
    log_printf(ERROR, "Two way on-the-fly transport has only been implemented "
                      "for ray tracing by z-stack");

  int num_2D_tracks = _track_generator_3D->getNum2DTracks();
  Track** flattened_tracks = _track_generator_3D->get2DTracksArray();
  int*** tracks_per_stack = _track_generator_3D->getTracksPerStack();
  int num_polar = _track_generator_3D->getNumPolar();
  int tid = omp_get_thread_num();

  /* Loop over flattened 2D tracks */
#pragma omp for schedule(dynamic)
  for (int ext_id=0; ext_id < num_2D_tracks; ext_id++) {

    /* Extract indices of 3D tracks associated with the flattened track */
    TrackStackIndexes tsi;
    Track* flattened_track = flattened_tracks[ext_id];
    tsi._azim = flattened_track->getAzimIndex();
    tsi._xy = flattened_track->getXYIndex();

    /* Loop over polar angles */
    for (int p=0; p < num_polar; p++) {

      /* Retrieve information for the first 3D Track in the z-stack */
      tsi._polar = p;
      tsi._z = 0;
      Track3D track_3D;
      _track_generator_3D->getTrackOTF(&track_3D, &tsi);

      if (kernel != NULL) {

        /* Reset kernel for a new base Track */
        kernel->newTrack(&track_3D);

        /* Trace all segments in the z-stack */
        traceStackTwoWay(flattened_track, p, kernel);
        track_3D.setNumSegments(kernel->getCount());
      }

      /* Operate on the Track */
      segment* segments = _track_generator_3D->getTemporarySegments(tid);
      onTrack(&track_3D, segments);
    }
  }
}


/**
 * @brief Traces the 3D segments of 3D Tracks in a z-stack both forward and
 *        backward across the geometry, applying the kernel provided by the
 *        user when the segment information is calculated.
 * @details This function copies information of the 3D z-stack, ray traces the
 *          z-stack forward using TrackGenerator::traceStackOTF, then reverses
 *          the tracks so that they point backwards, and ray traces in the
 *          reverse direction. This allows segments to be applied to
 *          TransportKernels during the on-the-fly ray tracing process.
 * @param flattened_track the 2D track associated with the z-stack for which
 *        3D segments are computed
 * @param polar_index the polar index of the 3D Track z-stack
 * @param kernel The TransportKernel applied to the calculated 3D segments
 */
void TraverseSegments::traceStackTwoWay(Track* flattened_track, int polar_index,
                                        TransportKernel* kernel) {

  /* Get segments from flattened track */
  segment* segments = flattened_track->getSegments();
  MOCKernel* moc_kernel = dynamic_cast<MOCKernel*>(kernel);

  /* Trace stack forwards */
  kernel->setDirection(true);
  traceStackOTF(flattened_track, polar_index, moc_kernel);
  kernel->post();

  /* Reverse segments in flattened track */
  int num_segments = flattened_track->getNumSegments();
  for (int s = 0; s < num_segments/2; s++) {
    segment tmp_segment = segments[num_segments-s-1];
    segments[num_segments-s-1] = segments[s];
    segments[s] = tmp_segment;
  }

  /* Flip CMFD surfaces on segments in flattened track */
  for (int s = 0; s < num_segments; s++) {
    int tmp_surface = segments[s]._cmfd_surface_fwd;
    segments[s]._cmfd_surface_fwd = segments[s]._cmfd_surface_bwd;
    segments[s]._cmfd_surface_bwd = tmp_surface;
  }

  /* Trace stack backwards */
  kernel->setDirection(false);
  traceStackOTF(flattened_track, polar_index, moc_kernel);
  kernel->post();

  /* Reverse segments in flattened track */
  for (int s = 0; s < num_segments/2; s++) {
    segment tmp_segment = segments[num_segments-s-1];
    segments[num_segments-s-1] = segments[s];
    segments[s] = tmp_segment;
  }

  /* Flip CMFD surfaces on segments in flattened track */
  for (int s = 0; s < num_segments; s++) {
    int tmp_surface = segments[s]._cmfd_surface_fwd;
    segments[s]._cmfd_surface_fwd = segments[s]._cmfd_surface_bwd;
    segments[s]._cmfd_surface_bwd = tmp_surface;
  }
}
