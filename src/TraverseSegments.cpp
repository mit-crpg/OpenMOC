#include "TraverseSegments.h"

/**
 * @brief Constructor for the TravseSegments class assigns the TrackGenerator
 */
TraverseSegments::TraverseSegments(TrackGenerator* track_generator) {

  /* Initialize kernel information to NULL */
  _kernels = NULL;
  _segments = NULL;
  _FSR_volumes = NULL;
  _FSR_locks = NULL;

  /* Save the track generator */
  _track_generator = track_generator;

  /* Extract tracking information from TrackGenerator */
  _tracks_2D = _track_generator->get2DTracks();
  _segment_formation = track_generator->getSegmentFormation();
  _max_optical_length = track_generator->retrieveMaxOpticalLength();
}


// description
TraverseSegments~TraverseSegments() {};


// TODO: DELETE: all this will be application specific 
void TraverseSegments::execute() {
  pre(); //make kernels, etc
  loopOverTracks();
  post(); //delete kernels, etc
}


// description
void TraverseSegments::allocateTemporarySegmentStorage() {

  int max_num_segments = _track_generator->getMaxNumSegments();
  
  /* Determine the number of segment arrays to allocate */
  int num_rows = 1;
  if (_segment_formation == OTF_STACKS)
    num_rows = _track_generator->getMaxNumTracksPerStack();

  /* Allocate memory */
  _segments = new segment*[num_rows];
  for (int z=0; z < num_rows; z++)
    _segments[z] = new segment[max_num_segments];
}



// description
void TraverseSegments::deallocateTemporarySegmentStorage() {

  int max_num_segments = _track_generator->getMaxNumSegments();

  /* Determine the number of segment arrays to deallocate */
  int num_rows = 1;
  if (_segment_formation == OTF_STACKS)
    num_rows = _track_generator->getMaxNumTracksPerStack();

  /* Deallocate memory */
  for (int z=0; z < num_rows; z++)
    delete [] _segments[z];
  delete [] _segments;
}


// description
//TODO: template
void TraverseSegments::allocateKernels<KernelType>() {

  /* Allocate temporary segments for SegmentationKernels */
  if (typeid(KernelType) == typeid(SegmentationKernel))
    allocateTemporarySegmentStorage();

  /* Determine the number of kernels to allocate */
  int num_rows = 1;
  if (_segment_formation == OTF_STACKS)
    num_rows = _track_generator->getMaxNumTracksPerStack();

  /* Allocate memory */
  _kernels = new MOCKernel*[num_rows];
  for (int z=0; z < num_rows; z++)
    kernels[z] = initializeKernel<KernelType>();
}


// description
CounterKernel* TraverseSegments::initializeKernel<CounterKernel>() {
  CounterKernel* kernel = new CounterKernel;
  kernel->setMaxVal(_max_optical_length);
  return kernel;
}


// description
VolumeKernel* TraverseSegments::initializeKernel<VolumeKernel>() {
  VolumeKernel* kernel = new VolumeKernel;
  kernel->setLocks(_FSR_locks);
  kernel->setBuffer(_FSR_volumes);
  return kernel;
}


// description
SegmentationKernel* TraverseSegments::initializeKernel<SegmentationKernel>() {
  SegmentationKernel* kernel = new SegmentationKernel;
  kernel->setMaxVal(_max_optical_length);
  kernel->setSegments(_segments);
  return kernel;
}


// description
//TODO: template
void TraverseSegments::deallocateKernels<KernelType>() {

  /* Deallocate temporary segments for SegmentationKernels */
  if (typeid(KernelType) == typeid(SegmentationKernel))
    deallocateTemporarySegmentStorage();

  /* Determine the number of kernels to deallocate */
  int num_rows = 1;
  if (_segment_formation == OTF_STACKS)
    num_rows = _track_generator->getMaxNumTracksPerStack();

  /* Deallocate memory */
  for (int z=0; z < num_rows; z++)
    delete _kernels[z];
  delete [] _kernels;
}


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
  for (int a=0; a < _num_azim/2; a++) {    
#pragma omp for
    int num_xy = _track_generator->getNumX(a) + _track_generator->getNumY(a);
    for (int i=0; i < num_xy; i++) {
      Track* track_2D = tracks_2D[a][i];
      onTrack(track_2D, 0);
    }
  }
}


// description
void TraverseSegments::loopOverTrackExplicit() { 


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

        /* Trace the segments on the track */ //FIXME
        _segments = track_3D->getSegments();

        /* Apply kernel to track */
        onTrack(track_3D, 0);
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

        /* Trace the segments on the track */
        traceSegmentsOTF(flattened_track, start, theta, _kernels[0]);

        /* Apply kernel to track */
        onTrack(track_3D, 0); // This will be application specific
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
      for (int z = 0; z < max_num_tracks_per_stack; z++)
        kernels[z]->resetCount();
      traceStackOTF(flattened_track, p, kernels);

      /* Loop over tracks in the z-stack */
      for (int z=0; z < tracks_per_stack[a][i][p]; z++) {

        /* Extract 3D track and initialize segments pointer */
        Track* track_3D = &tracks_3D[a][i][p][z];

        /* Apply kernel to track */
        onTrack(track_3D, z);
      }
    }
  }
}

