#include <math.h>     // Import the math library
#include "survS.h"    // Import the survS header file
#include "survproto.h"// Import the survproto header file

// Define the function survdiff2
void survdiff2(int *nn, int *nngroup, int *nstrat, 
               double *rho, double *time, int *status, 
               int *group, int *strata, double *obs, 
               double *exp, double *tempt, double *var, double *risk, 
               double *kaplan)
{
    // Register variables for fast access
    register int i, j, k; // Variables for loop iteration, register keyword suggests storing them in CPU registers for faster access
    int kk;               // Index variable for accessing the 2D array var
    int n;                // Number of observations in the current stratum
    int ngroup;           // Number of groups
    int ntot;             // Total number of observations
    int istart;           // Starting index of the current stratum
    int koff;             // Offset for the current stratum group, used for accessing obs and exp arrays
    double km;            // Kaplan-Meier estimate
    double nrisk;         // Number of samples in the current risk set
    double wt;            // Weight value, calculated based on Kaplan-Meier estimate and rho
    double tmp;           // Temporary variable for variance calculation
    double deaths;        // Number of deaths at the current time point
    double temp_i;        // Observed - expected value

    // Initialize variables
    ntot = *nn;             // Total number of samples
    ngroup = *nngroup;      // Number of groups
    istart = 0; koff = 0;   // Initial start index and offset

    // Initialize obs, exp, and tempt arrays to 0
    for (i = 0; i < *nstrat * ngroup; i++) {
        obs[i] = 0;
        exp[i] = 0;
        tempt[i] = 0;
    }

    // Loop through each stratum
    for (i = 0; i < ngroup; i++) risk[i] = 0; // Initialize risk array to 0

    // Find the last observation in the current stratum
    n = ntot;

    // Perform the actual test
    for (i = n - 1; i >= istart; i--) {
        // Calculate weight
        wt = 1;

        // Calculate the number of deaths
        deaths = 0;
        for (j = i; j >= istart && time[j] == time[i]; j--) {
            k = group[j] - 1;
            deaths += status[j];
            risk[k] += 1;
            obs[k + koff] += status[j] * wt;
        }

        i = j + 1;

        nrisk = n - i;

        // If there are deaths, calculate expected values and variance
        if (deaths > 0) {
            for (k = 0; k < ngroup; k++) {
                exp[k + koff] += wt * deaths * risk[k] / nrisk;
                temp_i = wt * (obs[k + koff] - exp[k + koff]);
                tempt[k + koff] += temp_i;
            }

            if (nrisk == 1) continue;  // No need to calculate variance if there's only one sample
            kk = 0;
            wt = wt * wt;
            for (j = 0; j < ngroup; j++) {
                tmp = wt * deaths * risk[j] * (nrisk - deaths) / (nrisk * (nrisk - 1));
                var[kk + j] += tmp;
                for (k = 0; k < ngroup; k++) {
                    var[kk] -= tmp * risk[k] / nrisk;
                    kk++;
                }
            }
        }
    }
}
