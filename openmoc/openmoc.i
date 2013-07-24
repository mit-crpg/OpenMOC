%module openmoc

%{
    #define SWIG_FILE_WITH_INIT
    #include "../src/Cell.h"
    #include "../src/Geometry.h"
    #include "../src/LocalCoords.h"
    #include "../src/log.h"
    #include "../src/Material.h"
    #include "../src/Point.h"
    #include "../src/Quadrature.h"
    #include "../src/Solver.h"
    #include "../src/CPUSolver.h"
    #include "../src/ThreadPrivateSolver.h"
    #include "../src/Surface.h"
    #include "../src/Timer.h"
    #include "../src/Track.h" 
    #include "../src/TrackGenerator.h"
    #include "../src/Universe.h"

    #define printf PySys_WriteStdout

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



/* If the user uses the --no-numpy flag, then NumPy typemaps will not be used
 * and the NumPy C API will not be embedded in the source code. The NumPy
 * typemaps are used to allow users to pass NumPy arrays to/from the C++ source
 * code. This does not work well on BGQ, and some users may prefer not to embed
 * this into their code, so if --no-numpy is passed in we use SWIG typemaps to
 * allow users to pass in arrays of data as Python lists. */
#ifdef NO_NUMPY

%include typemaps.i

/* Typemap for Lattice::setLatticeCells(int num_x, int num_y, short* universes) 
 * method - allows users to pass in Python lists of Universe IDs for each 
 * lattice cell */
%typemap(in) (int num_x, int num_y, short* universes) {

    if (!PyList_Check($input)) {
        PyErr_SetString(PyExc_ValueError,"Expected a Python list of integers "
			"lattice cells");
	return NULL;
    }

    $1 = PySequence_Length($input);  // num_x
    $2 = PySequence_Length(PyList_GetItem($input,0)); // num_y
    $3 = (short int*) malloc(($1 * $2) * sizeof(short));  // universes

    /* Loop over x */
    for (int i = 0; i < $2; i++) {

      /* Get the inner list in the nested list for the lattice */
        PyObject* outer_list = PySequence_GetItem($input,i);

	/* Check that the length of this list is the same as the length
	 * of the first list */
	if (PySequence_Length(outer_list) != $2) {
	    PyErr_SetString(PyExc_ValueError, "Size mismatch. Expected $1 x $2 "
			    "elements for lattice\n");
	    return NULL;
	}

	/* Loop over y */
        for (int j =0; j < $1; j++) {

	    /* Extract the value from the list at this location */
	    PyObject *o = PySequence_GetItem(outer_list,j);

	    /* If the value is a number, cast it as an int and set the
	     * input array value */
	    if (PyNumber_Check(o)) {
	        $3[i*$1 + j] = (short int) PyInt_AS_LONG(o);
	    } 
	    else {
	        free($3);
	        PyErr_SetString(PyExc_ValueError,"Expected a list of numbers as "
				"universe IDs when constructing lattice cells\n");
		return NULL;
	    }
	}
    }
}


/* If the user did not pass in the --no-numpy flag, then NumPy typemaps will be
 * used and the NumPy C API will be embedded in the source code. This will allow
 * users to pass arrays of data to/from the C++ source code (ie, setting group
 * cross-section values or extracting the scalar flux). */
#else

%include "numpy.i"


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

#endif


%include <exception.i> 
%include ../src/Cell.h
%include ../src/Geometry.h
%include ../src/LocalCoords.h
%include ../src/log.h
%include ../src/Material.h
%include ../src/Point.h
%include ../src/Quadrature.h
%include ../src/Solver.h
%include ../src/CPUSolver.h
%include ../src/ThreadPrivateSolver.h
%include ../src/Surface.h
%include ../src/Timer.h
%include ../src/Track.h
%include ../src/TrackGenerator.h
%include ../src/Universe.h


#ifdef DOUBLE
typedef double FP_PRECISION;
#else
typedef float FP_PRECISION;
#endif

