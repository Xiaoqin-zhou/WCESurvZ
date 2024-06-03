/* 
** This file was originally intended to support Splus and gradually evolved to support either R or Splus (based on ifdef lines). 
** It now supports only R.
*/
#include "R.h"                 // Include R's header file, which defines the basic R interface
#include "Rinternals.h"        // Include the header file for R's internal structures
#include <R_ext/Utils.h>       // Include the header file for R's extension utilities

/*
** Memory defined by ALLOC is automatically removed by S.
** Memory defined by "Calloc" needs to be removed by myself. The latter is used for objects that need to persist between calls.
*/
#define ALLOC(a,b)  R_alloc(a,b)  // Define the ALLOC macro, mapping it to R's memory allocation function R_alloc
#define CALLOC(a,b) R_Calloc(a,b) // Define the CALLOC macro, mapping it to R's memory allocation function R_Calloc
#define FREE(a)     R_Free(a)     // Define the FREE macro, mapping it to R's memory release function R_Free
