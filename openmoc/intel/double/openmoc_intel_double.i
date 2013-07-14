%module openmoc_intel_double

%{
    #define SWIG_FILE_WITH_INIT
    #include "../../../src/Cell.h"
    #include "../../../src/Geometry.h"
    #include "../../../src/LocalCoords.h"
    #include "../../../src/log.h"
    #include "../../../src/Material.h"
    #include "../../../src/Point.h"
    #include "../../../src/Quadrature.h"
    #include "../../../src/Solver.h"
    #include "../../../src/CPUSolver.h"
    #include "../../../src/ThreadPrivateSolver.h"
    #include "../../../src/VectorizedSolver.h"
    #include "../../../src/Surface.h"
    #include "../../../src/Timer.h"
    #include "../../../src/Track.h" 
    #include "../../../src/TrackGenerator.h"
    #include "../../../src/Universe.h"

    /* Exception helpers */
    static int swig_c_error_num = 0;
    static char swig_c_err_msg[512];

    const char* err_occurred(void) {
        if (swig_c_error_num) {
            swig_c_error_num = 0;
            return (const char*)swig_c_err_msg;
        }
        return NULL;
    }

    void set_err(const char *msg) {
        swig_c_error_num = 1;
        strncpy(swig_c_err_msg, msg, 256);
    }
%}

%warnfilter(506) log_printf(logLevel level, const char *format, ...);

%exception {
    try {
        $function
    } catch (const std::runtime_error &e) {
        SWIG_exception(SWIG_RuntimeError, err_occurred());
        return NULL;
    } catch (const std::exception &e) {
        SWIG_exception(SWIG_RuntimeError, e.what()); 
    }
}


%include "../../numpy.i"


%init %{
     import_array();
%}


%apply (int DIM1, int DIM2, short* IN_ARRAY2) {(int num_x, int num_y, short* universes)}
%apply (double* IN_ARRAY1, int DIM1) {(double* sigma_t, int num_groups)}
%apply (double* IN_ARRAY1, int DIM1) {(double* sigma_a, int num_groups)}
%apply (double* IN_ARRAY1, int DIM1) {(double* sigma_s, int num_groups)}
%apply (double* IN_ARRAY1, int DIM1) {(double* sigma_f, int num_groups)}
%apply (double* IN_ARRAY1, int DIM1) {(double* nu_sigma_f, int num_groups)}
%apply (double* IN_ARRAY1, int DIM1) {(double* chi, int num_groups)}
%apply (float* IN_ARRAY1, int DIM1) {(float* sigma_t, int num_groups)}
%apply (float* IN_ARRAY1, int DIM1) {(float* sigma_a, int num_groups)}
%apply (float* IN_ARRAY1, int DIM1) {(float* sigma_s, int num_groups)}
%apply (float* IN_ARRAY1, int DIM1) {(float* sigma_f, int num_groups)}
%apply (float* IN_ARRAY1, int DIM1) {(float* nu_sigma_f, int num_groups)}
%apply (float* IN_ARRAY1, int DIM1) {(float* chi, int num_groups)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* coords, int num_tracks)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* coords, int num_segments)}

%include <exception.i> 
%include ../../../src/Cell.h
%include ../../../src/Geometry.h
%include ../../../src/LocalCoords.h
%include ../../../src/log.h
%include ../../../src/Material.h
%include ../../../src/Point.h
%include ../../../src/Quadrature.h
%include ../../../src/Solver.h
%include ../../../src/CPUSolver.h
%include ../../../src/ThreadPrivateSolver.h
%include ../../../src/VectorizedSolver.h
%include ../../../src/Surface.h
%include ../../../src/Timer.h
%include ../../../src/Track.h
%include ../../../src/TrackGenerator.h
%include ../../../src/Universe.h

typedef double FP_PRECISION;

