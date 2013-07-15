/**
 * @file VectorizedSolver.h
 * @brief The VectorizedSolver class.
 * @date May 28, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef VECTORIZEDSOLVER_H_
#define VECTORIZEDSOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <mm_malloc.h>
#include <mkl.h>
#include "CPUSolver.h"
#endif

#define VEC_LENGTH 8


/**
 * @class VectorizedSolver VectorizedSolver.h "openmoc/src/host/VectorizedSolver.h"
 * @brief This is a subclass of the CPUSolver class which uses memory-aligned
 *        data structures and Intel's auto-vectorization.
 */
class VectorizedSolver : public CPUSolver {

private:

    /** The vector length for vector aligned data arrays */
    int _vector_length;

    /** The vector alignment (power of 2) for aligned data arrays */
    int _vector_alignment;

    /** Number of energy groups divided by vector widths (VEC_LENGTH) */
    int _num_vector_lengths;

    void initializeFluxArrays();
    void initializeSourceArrays();

    void normalizeFluxes();
    FP_PRECISION computeFSRSources();
    void scalarFluxTally(segment* curr_segment, 
			 FP_PRECISION* track_flux,
			 FP_PRECISION* fsr_flux);
    void transferBoundaryFlux(int track_id, bool direction,
			      FP_PRECISION* track_flux);
    void addSourceToScalarFlux();
    void computeKeff();
    void transportSweep();

public:
    VectorizedSolver(Geometry* geom=NULL, TrackGenerator* track_generator=NULL);
    virtual ~VectorizedSolver();
 
    int getVectorLength();
    int getVectorAlignment();
    int getNumVectorWidths();

    void setGeometry(Geometry* geometry);
};


#endif /* VECTORIZEDSOLVER_H_ */
