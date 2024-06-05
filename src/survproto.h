/*
** Prototypes for all survival functions
** Including them in each routine helps prevent parameter errors
* These function declarations help prevent errors caused by parameter mismatches between function definitions and calls
*/

// Define a function to convert a one-dimensional double array to a two-dimensional array
double **dmatrix(double *array, int nrow, int ncol);

// Define a function to convert a one-dimensional int array to a two-dimensional array
int    **imatrix(int *array, int nrow, int ncol);

// Function to initialize a loop, receiving min and max values as parameters
void init_doloop(int min, int max);

// Function to execute a loop, receiving the number of loops and an index array as parameters
int doloop(int nloops, int *index);

// Function to calculate the number of individuals at risk
// Receives parameters: number of individuals, start time, end time, status, sort index 1, sort index 2, stratification information
int *norisk(int n, double *time1, double *time2, double *status, 
            int *sort1, int *sort2, int *strata);

// Function to calculate the P value step
// Receives parameters: number of columns, index array, index array 2, weight array, data array, factor array, dimension array, cut array, step size, edge
double pystep(int nc, int *index, int *index2, double *wt, 
              double *data, int *fac, int *dims, double **cuts, 
              double step, int edge);

// Define a structure snode
typedef struct snode {
    int value;              // Node value
    int depth;              // Number of forward links
    struct snode *forward[1];  // Set of forward links
} snode;

// Define the survfitkm function for calculating survival curves
// This function calculates survival curves and uses R's SEXP type for parameters to interact with R
// Receives multiple SEXP type parameters representing data passed from R
// This function receives multiple parameters representing survival time, weight, sort index, type, ID, number of groups, position, influence, reverse flag, and entry flag
SEXP survfitkm(SEXP y2, SEXP weight2, SEXP sort12, SEXP sort22, 
               SEXP type2, SEXP id2, SEXP nid2, SEXP position2, 
               SEXP influence2, SEXP reverse2, SEXP entry2);

// void survdiff2(int   *nn,     int   *nngroup,    int   *nstrat, 
// 	       double *rho,    double *time,       int   *status, 
// 	       int   *group,  int   *strata,	   double *obs, 
// 	       double *exp, double *tempt,    double *var,        double *risk, 
// 	       double *kaplan);

void wcelogrank(int   *nn,     int   *nngroup,    
	          double *time,       int   *status, int   *group, double *weight, 
              double *nrisk1, double *nrisk2,  double *nriskall,
              double *obs, double *exp, double *tempt,   double *var,   double *risk);

// void survWCEKM2(int *nn, int *nngroup,    
// 	          double *time, int *group,
//               double *nevent, double *nrisk, 
//               double *estimate, double *stderr, double *cumhaz, double *stdchaz);

void survWCEKM(int *nn, int *nngroup,    
	          double *time, int *group,
              double *nevent, double *nrisk, 
              double *estimate, double *stderr2, double *cumhaz2, double *stdchaz2);


