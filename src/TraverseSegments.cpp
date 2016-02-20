#include "TraverseSegments.h"

/**
 * @brief Constructor for the TravseSegments class assigns the TrackGenerator
 */
TraverseSegments::TraverseSegments(TrackGenerator* track_generator) {

  /* Initialize kernel information to NULL */
  _kernels = NULL;

  /* Save the track generator */
  _track_generator = track_generator;

  /* Determine the type of segment formation used */
  _segment_formation = track_generator->getSegmentFormation();
}


// description
TraverseSegments::~TraverseSegments() {};


// TODO: DELETE: all this will be application specific 
/*
void TraverseSegments::execute() {
  pre(); //make kernels, etc
  loopOverTracks();
  post(); //delete kernels, etc
}
*/
void TraverseSegments::execute() {}
void TraverseSegments::onTrack(Track* track, segment* segments) {}


// description
void TraverseSegments::loopOverTracks() {

  switch (_segment_formation) {
    case TWO_DIM:
      loopOverTracks2D();
      break;
    case EXPLICIT:
      loopOverTracksExplicit();
      break;
    case OTF_TRACKS:
      loopOverTracksByTrackOTF();
      break;
    case OTF_STACKS:
      loopOverTracksByStackOTF();
      break;
  }
}


// description
void TraverseSegments::loopOverTracks2D() { 

  //FIXME: loop over all 2D tracks?
  Track2D** tracks_2D = _track_generator->get2DTracks();
  int num_azim = _track_generator->getNumAzim();
  for (int a=0; a < num_azim/2; a++) {
    int num_xy = _track_generator->getNumX(a) + _track_generator->getNumY(a);
#pragma omp for
    for (int i=0; i < num_xy; i++) {
      Track* track_2D = &tracks_2D[a][i];
      _kernels[0]->newTrack(track_2D);
      _track_generator->traceSegmentsExplicit(track_2D, _kernels[0]);
      segment* segments = track_2D->getSegments();
      onTrack(track_2D, segments);
    }
  }
}


// description
void TraverseSegments::loopOverTracksExplicit() {

  int num_2D_tracks = _track_generator->getNum2DTracks();
  Track** flattened_tracks = _track_generator->getFlattenedTracksArray();
  Track3D**** tracks_3D = _track_generator->get3DTracks();
  int*** tracks_per_stack = _track_generator->getTracksPerStack();
  int num_polar = _track_generator->getNumPolar();

#pragma omp for
  /* Loop over flattened 2D tracks */
  for (int ext_id=0; ext_id < num_2D_tracks; ext_id++) {
    
    /* Extract indices of 3D tracks associated with the flattened track */
    Track* flattened_track = flattened_tracks[ext_id];
    int a = flattened_track->getAzimIndex();
    int i = flattened_track->getXYIndex();

    /* Loop over polar angles */
    for (int p=0; p < num_polar; p++) {

      /* Loop over tracks in the z-stack */
      for (int z=0; z < tracks_per_stack[a][i][p]; z++) {

        /* Extract 3D track and initialize segments pointer */
        Track* track_3D = &tracks_3D[a][i][p][z];
        _kernels[0]->newTrack(track_3D);

        /* Trace the segments on the track */
        _track_generator->traceSegmentsExplicit(track_3D, _kernels[0]);
        segment* segments = track_3D->getSegments();

        /* Apply kernel to track */
        onTrack(track_3D, segments);
      }
    }
  }
}


// description
void TraverseSegments::loopOverTracksByTrackOTF() {

  int num_2D_tracks = _track_generator->getNum2DTracks();
  Track** flattened_tracks = _track_generator->getFlattenedTracksArray();
  Track3D**** tracks_3D = _track_generator->get3DTracks();
  int*** tracks_per_stack = _track_generator->getTracksPerStack();
  int num_polar = _track_generator->getNumPolar();
  int tid = omp_get_thread_num();

#pragma omp for
  /* Loop over flattened 2D tracks */
  for (int ext_id=0; ext_id < num_2D_tracks; ext_id++) {
    
    /* Extract indices of 3D tracks associated with the flattened track */
    Track* flattened_track = flattened_tracks[ext_id];
    int a = flattened_track->getAzimIndex();
    int i = flattened_track->getXYIndex();

    /* Loop over polar angles */
    for (int p=0; p < num_polar; p++) {

      /* Loop over tracks in the z-stack */
      for (int z=0; z < tracks_per_stack[a][i][p]; z++) {

        /* Extract 3D track and initialize segments pointer */
        Track* track_3D = &tracks_3D[a][i][p][z];
        _kernels[0]->newTrack(track_3D);
        double theta = tracks_3D[a][i][p][z].getTheta();
        Point* start = track_3D->getStart();

        /* Trace the segments on the track */
        _track_generator->traceSegmentsOTF(flattened_track, start, theta,
                                            _kernels[0]);
        track_3D->setNumSegments(_kernels[0]->getCount());
        segment* segments = _track_generator->getTemporarySegments(tid, 0);

        /* Apply kernel to track */
        onTrack(track_3D, segments);
      }
    }
  }
}


// description
void TraverseSegments::loopOverTracksByStackOTF() {

  int num_2D_tracks = _track_generator->getNum2DTracks();
  Track** flattened_tracks = _track_generator->getFlattenedTracksArray();
  Track3D**** tracks_3D = _track_generator->get3DTracks();
  int*** tracks_per_stack = _track_generator->getTracksPerStack();
  int num_polar = _track_generator->getNumPolar();
  int tid = omp_get_thread_num();

#pragma omp for
  /* Loop over flattened 2D tracks */
  for (int ext_id=0; ext_id < num_2D_tracks; ext_id++) {
    
    /* Extract indices of 3D tracks associated with the flattened track */
    Track* flattened_track = flattened_tracks[ext_id];
    int a = flattened_track->getAzimIndex();
    int i = flattened_track->getXYIndex();

    /* Loop over polar angles */
    for (int p=0; p < num_polar; p++) {
      
      /* Trace all tracks in the z-stack */      
      for (int z = 0; z < tracks_per_stack[a][i][p]; z++)
        _kernels[z]->newTrack(&tracks_3D[a][i][p][z]);
      _track_generator->traceStackOTF(flattened_track, p, _kernels);

      /* Loop over tracks in the z-stack */
      for (int z=0; z < tracks_per_stack[a][i][p]; z++) {

        /* Extract 3D track and initialize segments pointer */
        Track* track_3D = &tracks_3D[a][i][p][z];
        track_3D->setNumSegments(_kernels[z]->getCount());
        segment* segments = _track_generator->getTemporarySegments(tid, z);

        /* Apply kernel to track */
        onTrack(track_3D, segments);
      }
    }
  }
}

