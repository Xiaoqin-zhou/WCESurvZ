/*
** Routine to handle raw data for survfit, not applicable for multi-state data or Cox model
** Reference: survival:survfitkm
**
** Variables passed to R are named "zed2", with variable names as "zed" in the code
** y: A matrix with 2 or 3 columns for survival times
** weight: A vector of case weights
** sort1: A vector of sorted indices for start times (if y has 3 columns), otherwise null
** sort2: A vector of sorted indices for end times
** type: Type of survival calculation
** id: Group identifier for clustered variance. 0,1,2,...
** nid: Number of unique groups
** position: Marker for multiple observations, e.g., (1,10), (10,14), (14,20)
**           1*(first observation) + 2*(last observation)
**           Value of 3 for individuals with only one observation.
** influence: Specifies which influence matrix to return, 1*(cumulative hazard) + 2*(survival)
** reverse: 0=ordinary KM, 1=estimate censoring distribution (unfinished)
** entry: 1=retain rows with unique entry times
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

SEXP survfitkm(SEXP y2, SEXP weight2, SEXP sort12, SEXP sort22, 
               SEXP type2, SEXP id2, SEXP nid2, SEXP position2,
               SEXP influence2, SEXP reverse2, SEXP entry2) {
              
    // Declare variables
    int i, i1, i2, j, k, person1, person2; // General index variables
    int nused, nid, type, influence; // Number of observations used, group count, survival analysis type, influence variable
    int ny, ntime; // Number of columns in input data and count of unique time points
    int reverse, entry; // Flags for reversing and recording entry time
    double* time1 = 0, * time2, * status, * wt; // Pointers for time, status, and weights
    double v1, v2, dtemp, haz; // Intermediate and hazard variables
    double temp, dtemp2, dtemp3, frac, btemp; // Temporary and fraction variables
    double d0, d1, nrisk; // Death count, weighted death count, and risk count
    int* sort1 = 0, * sort2, * id = 0; // Pointers for sort indices and group identifiers
    static const char *outnames[]={"time", "n", "estimate", "std.err",
                                     "influence1", "influence2", ""};  // Output variable names
    SEXP rlist; // R list of SEXP type
    double* gwt = 0, * inf1 = 0, * inf2 = 0, * inf3 = 0; // Working vectors and influence matrix pointers
    int* gcount = 0; // Pointer for group counter
    int n1, n2, n3, n4; // Risk, event, censor, and add counts
    int* position = 0; // Pointer for position markers
    double wt1, wt2, wt3, wt4; // Weighted versions of risk, event, censor, and add counts
                      
    /* Output variables */
    double  *n[8],  *dtime,  // n array and unique time pointer
            *kvec, *nvec, *std[3], *imat1=0, *imat2=0; // Current estimates, standard errors, and influence matrix pointers /* =0 to avoid -Wall warning */
    double km, nelson; // KM estimate and Nelson-Aalen estimate /* current estimates */

    /* Map input data */
    ny = ncols(y2);     /* Get the number of columns in y2, 2 for regular survival data, 3 for start-stop data */
    nused = nrows(y2);  /* Get the number of rows in y2, which is the amount of data used */
    
    time2 = REAL(y2);        /* If y2 has two columns, map y2 to a C double pointer pointing to the time column */
        
    status = time2 + nused;      /* status points to the status column, located nused elements after time2 */
    wt = REAL(weight2);          /* Map weight2 to a C double pointer pointing to the weight column */
    sort2 = INTEGER(sort22);     /* Map sort22 to a C integer pointer */
    nused = LENGTH(sort22);      /* Get the length of sort22, which is the actual data used */

    type = asInteger(type2);     /* Convert type2 to an integer, representing the type of survival analysis */
    nid = asInteger(nid2);       /* Convert nid2 to an integer, representing the number of unique groups */
    position = INTEGER(position2); /* Map position2 to a C integer pointer pointing to position markers */
    influence = asInteger(influence2); /* Convert influence2 to an integer, representing the type of influence matrix */
    reverse = asInteger(reverse2);     /* Convert reverse2 to an integer, indicating whether to reverse */
    entry = asInteger(entry2);         /* Convert entry2 to an integer, indicating whether to record entry times */

    /* nused serves two purposes. The first is the length of the input data y, used only to set time1, time2, and status.
       The second is the number of observations we actually use, which is the length of sort2.
       This routine can be called multiple times, with sort1/sort2 pointing to different subsets of data, while y, wt, id, and position remain constant.
    */

    /* First step, get the number of unique times, for memory allocation.
       Provided with the number of groups (unique id values).
       I count both unique start and end times.
    */
    // Count the unique time2 values (i.e., end times), ignoring the entry flag
    ntime = 1; // Initialize the count of unique time points to 1
    temp = time2[sort2[0]];  /* Set temp to the first sorted end time, the smallest end time */

    // Traverse all used data points starting from the second point, to get unique time points for constructing intervals or calculating survival probabilities
    for (i = 1; i < nused; i++) {
        i2 = sort2[i]; // Get the index of the ith sorted end time

        // If the current end time is different from temp and the position marker indicates a valid end time or the status is non-zero (indicating an event)
        if ((position[i2] > 1 || status[i2] > 0) && time2[i2] != temp) {
            ntime++; // Increase the count of unique time points
            temp = time2[i2]; // Update temp to the current end time
        }
    }

 
    /* Allocate memory for output
        n has 6 columns representing risk count, event count, censored count,
        followed by 3 weighted versions, and optionally two columns for risk set additions (when entry=1)
    */
    // PROTECT to safeguard the R object rlist from garbage collection
    PROTECT(rlist = mkNamed(VECSXP, outnames));

    // Allocate memory for dtime and add it to rlist
    dtime = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, ntime)));
    j = 6;

    // Allocate memory for n[0] as a ntime * j matrix and add it to rlist
    n[0] = REAL(SET_VECTOR_ELT(rlist, 1, allocMatrix(REALSXP, ntime, j)));

    // Set pointers for each element in the n array
    for (i = 1; i < j; i++) {
        n[i] = n[0] + i * ntime;
    }

    // Allocate memory for kvec as a ntime * 2 matrix and add it to rlist
    kvec = REAL(SET_VECTOR_ELT(rlist, 2, allocMatrix(REALSXP, ntime, 2)));

    // nvec points to the second column of the kvec matrix, representing the Nelson-Aalen estimate
    nvec = kvec + ntime;

    // Allocate memory for std as a ntime * ny matrix and add it to rlist
    // std[0] represents the standard error of survival rate, std[1] represents the standard error of cumulative hazard, std[2] represents the standard error of AUC
    std[0] = REAL(SET_VECTOR_ELT(rlist, 3, allocMatrix(REALSXP, ntime, ny)));
    std[1] = std[0] + ntime;

    // If entry flag is not 1, only calculate end times
    temp = time2[sort2[0]];  /* Set temp to the first sorted end time, the smallest end time */
    dtime[0] = temp; // Set the first element of dtime array to temp
    k = 1; // Initialize the index k for dtime to 1

    // Traverse all used data points starting from the second point, this step mainly gets dtime (unique time points of both groups)
    for (i = 1; i < nused; i++) {
        i2 = sort2[i]; // Get the index of the ith sorted end time

        // If the current end time is different from temp and the position marker indicates a valid end time or the status is non-zero (indicating an event)
        if ((position[i2] > 1 || status[i2] > 0) && time2[i2] != temp) {
            temp = time2[i2]; // Update temp to the current end time
            dtime[k++] = temp; // Store the current end time in the dtime array and increment k
        }
    }

 
    /*
    ** Next, compute all the counts
    ** Temporary variables n1 = risk count, n2 = event count, n3 = censor count,
    **   n4 = addition count. wt1, wt2, wt3, wt4 = weighted versions of n1 to n4.
    **   Store them in n[][0-5] to return to the R routine.
    ** person1, person2 track sort1 and sort2 respectively,
    **   i1 and i2 likewise.
    */
    // Initialize person1 and person2 to nused-1, the index of the last data point
    person1 = nused - 1;
    person2 = nused - 1;

    // Initialize n1 (risk count) and wt1 (weighted risk count) to 0
    n1 = 0;
    wt1 = 0;

    // printf("ntime:%d",ntime);   ntime is the number of unique time points
    // printf("nused:%d",nused);   nused is the number of time points within the group

    // Traverse time points in reverse, starting from the last time point
    for (k = ntime - 1; k >= 0; k--) {
        // Initialize n2 (event count), n3 (censor count), wt2 (weighted event count), and wt3 (weighted censor count) to 0
        n2 = 0;
        n3 = 0;
        wt2 = 0;
        wt3 = 0;

        // Traverse data points in reverse, sorted by end time
        for (; person2 >= 0; person2--) {
            i2 = sort2[person2]; // Get the index of the person2-th sorted end time

            // If the current end time is less than the current time point, break the loop
            if (time2[i2] < dtime[k]) break;

            // Increase the risk count
            n1++;
            wt1 += wt[i2];    // Increase the weighted risk count

            // If the current status indicates an event
            if (status[i2] == 1) {
                n2++;            // Increase the event count
                wt2 += wt[i2];   // Increase the weighted event count
            }
            else if (position[i2] & 2) {
                /*
                ** If this is the last of the interval (a,b](b,c](c,d].. for
                ** the individual (position[i2]=2 or 3), it is a "real" censor
                */
                n3++;           // Increase the censor count
                wt3 += wt[i2];  // Increase the weighted censor count
            }
        }

        // Store the risk count, event count, censor count, and corresponding weighted values for the current time point in the n array
        n[0][k] = n1;
        n[1][k] = n2;
        n[2][k] = n3;
        n[3][k] = wt1;
        n[4][k] = wt2;
        n[5][k] = wt3;
    }

    // R_CheckUserInterrupt();  /* Check for Ctrl-C interrupt */


    /*
    ** Calculate survival rates and cumulative hazards, simple variance
    ** Reference: survfitkm: estimates
    */
    // Initialize Nelson-Aalen hazard and KM survival rate
    nelson = 0.0;
    km = 1.0;

    // Initialize variances v1 and v2
    v1 = 0;
    v2 = 0;

    printf("reverse:%f",reverse);
    // Traverse all time points
    for (i = 0; i < ntime; i++) {
        d0 = n[1][i];    /* Unweighted event count */
        d1 = n[4][i];    /* Weighted event count */
        nrisk = n[3][i]; /* Risk count */

        // Handle Nelson-Aalen hazard (type is 1 or 3)
        if (type == 1 || type == 3) {
            if (d0 > 0 && d1 > 0) {  /* At least one event, weight > 0 */
                nelson += d1 / nrisk; // Update Nelson-Aalen hazard
                v2 += d1 / (nrisk * nrisk); // Update variance v2
            }
            nvec[i] = nelson; // Store Nelson-Aalen estimate
            std[1][i] = sqrt(v2); // Store standard error
        }

        // Handle KM survival rate (type less than 3)
        if (type < 3) {
            if (d0 > 0 && d1 > 0) {  /* At least one event */
                km *= (nrisk - d1) / nrisk; // Update KM survival rate
                v1 += d1 / (nrisk * (nrisk - d1)); /* Greenwood's formula for variance */
            }
            kvec[i] = km; // Store KM survival rate
            std[0][i] = sqrt(v1); // Store variance
        }
    }

    UNPROTECT(1); // Release protection for rlist
    return (rlist); // Return the result list
}

