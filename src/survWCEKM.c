#include <math.h>
#include "survS.h"
#include "survproto.h"
#include <stdio.h>

// Define the function survWCEKM
void survWCEKM(int *nn, int *nngroup,    
	          double *time, int *group,
              double *nevent, double *nrisk, 
              double *estimate, double *stderr2, double *cumhaz2, double *stdchaz2)
{
    int i, k;
    int n = *nn;  // Number of time points
    int ngroup = *nngroup;  // Number of groups
    double km;            // Kaplan-Meier estimate
    double v1;    // Initialize variance for standard error calculation
    double nelson;// Initialize cumulative hazard
    double v2;    // Initialize variance for cumulative hazard standard error
    int istart = 0;

    // Initialize variables
    for (i = 0; i < n; i++) {
        estimate[i] = 1.0;    // Initialize survival probability to 1
        stderr2[i] = 0.0;  // Initialize standard error to 0
        cumhaz2[i] = 0.0;  // Initialize cumulative hazard to 0
        stdchaz2[i] = 0.0; // Initialize standard error of cumulative hazard to 0
    }

    
    // Calculate survival, standard error, cumulative hazard, and its standard error for each group
    for (k = 0; k < ngroup; k++) {
        km = 1.0;    // Initialize survival probability to 1 for each group
        v1 = 0.0;    // Initialize variance for standard error calculation
        nelson = 0.0;// Initialize cumulative hazard
        v2 = 0.0;    // Initialize variance for cumulative hazard standard error

        
        for (i = istart; i < n; i++) {
            
            if (group[i] != k + 1) break;
            if (nrisk[i] > 0) {
                km *= (nrisk[i] - nevent[i]) / nrisk[i];  // Update survival probability
            }
            estimate[i] = km;  // Store survival probability

            v1 += (nrisk[i] > 0) ? nevent[i] / (nrisk[i] * (nrisk[i] - nevent[i])) : 0; // Update variance sum
            stderr2[i] = sqrt(v1); // Store standard error

            nelson += (nrisk[i] > 0) ? nevent[i] / nrisk[i] : 0; // Update cumulative hazard
            cumhaz2[i] = nelson; // Store cumulative hazard

            v2 += (nrisk[i] > 0) ? nevent[i] / (nrisk[i] * nrisk[i]) : 0; // Update variance sum for cumulative hazard
            stdchaz2[i] = sqrt(v2); // Store standard error of cumulative hazard
        }
        istart =i;
    }
}
