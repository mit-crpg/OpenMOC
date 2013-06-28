%module openmoc_mic

%{
    #define SWIG_FILE_WITH_INIT
    #include "../../src/dev/mic/MICQuery.h"
    #include "../../src/dev/mic/MICSolver.h"

%}


%include <exception.i> 
%include ../../src/dev/mic/MICSolver.h
%include ../../src/dev/mic/MICQuery.h

#ifdef DOUBLE
typedef double FP_PRECISION
#else
typedef float FP_PRECISION;
#endif
