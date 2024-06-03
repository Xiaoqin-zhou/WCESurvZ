#include <math.h>     // Import the math library
#include "survS.h"    // Import the survS header file
#include "survproto.h"// Import the survproto header file

// Define the function wcelogrank
void wcelogrank(int *nn, int *nngroup,    // Total sample size and total number of groups. Pointers are used to import data from R to C due to different data structures.
               double *time, int *status, int *group, double *weight,    // Arrays of length *nn for time, status, and weight. *group is an array representing group labels.
               double *nrisk1, double *nrisk2, double *nriskall,    // Arrays for number at risk in each group and overall.
               double *obs, double *exp, double *tempt, double *var, double *risk   // Arrays of length nngroup for observed values, expected values, variance, and risk.
               )
{
    // Register variables for fast access
    register int i, j, k; // Loop iterator variables. The register keyword suggests storing them in CPU registers for faster access.
    int kk;               // Index variable for accessing the 2D array var
    int n;                // Number of observations in the current stratum
    int ngroup;           // Number of groups
    int ntot;             // Total number of observations
    double km;            // Kaplan-Meier estimate
    double nrisk;         // Number of samples in the current risk set
    double wt;            // Weight value calculated based on the Kaplan-Meier estimate and rho
    double tmp;           // Temporary variable for variance calculation
    double deaths;        // Number of deaths at the current time point

    // Assign variables using pointers to convert R data format to C data format. For array pointers, assignment is not needed.
    n = *nn;  // n = total sample size
    ngroup = *nngroup;  // Number of groups
    var[0] = 0;  // Initialize the variance array

    // Initialize observed values, expected values, and risk arrays to 0
    for (i = 0; i < ngroup; i++) {
        obs[i] = 0;  // Initialize observed values to 0
        exp[i] = 0;  // Initialize expected values to 0
        risk[i] = 0;  // Initialize risk array to 0
    }

    // Perform the actual test
    for (i = n - 1; i >= 0; i--) {  
        
        k = group[i] - 1;  // Group index
        obs[k] += status[i] * weight[i];  // Update observed values

        exp[0] +=  (weight[i] * nrisk1[i]) / (nriskall[i]);  // Update expected values for group 1
        exp[1] +=  (weight[i] * nrisk2[i]) / (nriskall[i]);  // Update expected values for group 2

        tmp =  (weight[i] * nrisk1[i] * nrisk2[i] * (nriskall[i] - 1)) / ((nriskall[i] * nriskall[i]) * (nriskall[i] - weight[i]));
        var[0] += tmp;  // Update variance
    }
}
