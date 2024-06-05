/*
** This file preloads the entry points of .C routines
** Added as per the request of R-core
** It adds an extra layer of protection by declaring the number of parameters,
** and may bring a slight speed improvement
*/
#include "survS.h" // Include custom header file
#include "survproto.h" // Include custom header file
#include "R_ext/Rdynload.h" // Include R dynamic loading header file
#include "Rversion.h" // Include R version information header file

static const R_CMethodDef Centries[] = {
    // {"Csurvdiff2",  (DL_FUNC) &survdiff2, 14},
    {"Cwcelogrank",  (DL_FUNC) &wcelogrank, 14},
    {"CsurvWCEKM",  (DL_FUNC) &survWCEKM, 10},
    {NULL, NULL, 0}
};

// Static constant array of R_CallMethodDef, defining the .C routines to register
// This is a static constant array of type R_CallMethodDef, used to define the .Call routines to register
static const R_CallMethodDef Callentries[] = {
    {"Csurvfitkm",    (DL_FUNC) &survfitkm,   11},  // Define the "Csurvfitkm" function, pointing to survfitkm, accepting 11 parameters
    {NULL, NULL, 0}  // End marker for the array
};

// Function called when initializing the WCESurvZ package
void R_init_WCESurvZ(DllInfo *dll){
    // Register R routines, the second parameter being NULL indicates no .C routines, the third parameter indicates the .Call routines to register
    R_registerRoutines(dll, Centries,  Callentries, NULL, NULL);

    /* 
    ** The following line makes it so that only the routines defined above can be used by external packages,
    ** meaning internal functions like dmatrix() are now invisible.
    */
    R_useDynamicSymbols(dll, FALSE); 

    /*
    ** This line makes them accessible only through the symbols defined above,
    ** meaning .Call("tmerge", ) will not work, but .Call(Ctmerge, ) will.
    ** This feature was added in version 2.16
    */
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
