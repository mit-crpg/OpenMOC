/** Helper typemaps to enable user input to functions that require pointers
    as input */

%module argument_typemaps

/* Use swig typemaps */
%include <typemaps.i>

/* Typemaps for the intrinsic exponentials */
void expF1_poly(FP_PRECISION INPUT, FP_PRECISION* OUTPUT);
void expF2_poly(FP_PRECISION INPUT, FP_PRECISION* OUTPUT);
void expH_poly(FP_PRECISION INPUT, FP_PRECISION* OUTPUT);
void expF1_fractional(FP_PRECISION INPUT, FP_PRECISION* OUTPUT);
void expF2_fractional(FP_PRECISION INPUT, FP_PRECISION* OUTPUT);
void expH_fractional(FP_PRECISION INPUT, FP_PRECISION* OUTPUT);
void expG_fractional(FP_PRECISION INPUT, FP_PRECISION* OUTPUT);
void expG2_fractional(FP_PRECISION INPUT, FP_PRECISION* OUTPUT);
void cram7(FP_PRECISION INPUT, FP_PRECISION* OUTPUT);
