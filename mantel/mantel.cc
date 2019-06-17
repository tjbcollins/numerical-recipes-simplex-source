/*
  This is an interactive (manual) version of telios (hence mantel), my implementation 
  of Numerical Recipes' downhill simplex method for function minimization.

  Usage: mantel <list of initial parameters>

  So e.g. if you are minimizing y(x1, x2, x3) then you will invoke mantel with:

    mantel x10 x20 x30

  Mantel then gives you a lists of parameters and prompts you for the corresponding
  value of y for this list.

  tjbc 10/8/2014
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "/home/tcol/utils/mylib.h"
#include "/theory/tcol/NR3/code/nr3.h"
#include "/theory/tcol/NR3/code/amoeba.h"

static const int MAXSTR=10000;

Int iteration=0;
Int Ndim;
Doub latest;
FILE *flog;
int UsingScript = 0; 
char script[MAXSTR];

int main(int argc, char **argv) {
    Doub tol = 1.0e-4;
    Doub y(VecDoub v);
    Int i;
    int I;
    char logFile[] = "mantel.log";

    flog = safeopen(logFile, "w");

    if (argc > 1 && !strcmp(argv[1], "-S")) {
	UsingScript = 1;
	strcpy(script, argv[2]);
	fprintf(stderr, "Using script '%s'\n", script);
    }

    if ((!UsingScript && argc < 2) || (UsingScript && argc < 4)) {
	fprintf(stderr, "Error: Incorrect number of arguments.\n");
	fprintf(stderr, "Usage: %s [-S script] x1 x2 ... xn\n", argv[0]);
	fprintf(stderr, "\nYou specify the initial set of parameters {xn}, and %s\n", argv[0]);
	fprintf(stderr, "gives you successive parameter vectors and prompts you for\n");
	fprintf(stderr, "the corresponding metric values y(x1, ..., xn), minimizing y in the process.\n"
		"-S: read in the metric values from mantel.in, write out the new vector to mantel.out,\n"
		"    and run script 'script' after each step\n");
	exit(1);
    }

    if (UsingScript)
	Ndim = argc-3;
    else
	Ndim = argc-1; // Number of dimensions (remember a simplex has N+1 vertices in an N dimensional space)

    printf("Using a %d-point simplex in an %d-dimensional space.\n", Ndim+1, Ndim);
    printf("Steps are recorded in the file '%s'\n", logFile);
    if (!UsingScript) printf("Enter values for y(x1, ..., xn), or 'q' to quit.\n");

    VecDoub v(Ndim);      // This is the 1-D array containing the points, NR3-style
    VecDoub dv(Ndim);     // This is the initial perturbation to these points
    VecDoub vFinal(Ndim); // This is the converged solution

    // Read the commandline parameters
    for (i=0; i < Ndim; i++) {
	I = UsingScript ? i+2 : i;
	if (sscanf(argv[I+1], "%le", &(v[i])) != 1) {
	    fprintf(stderr, "Error interpreting %s as a number.\n", argv[i+1]);
	    exit(1);
	}
	dv[i] = 0.05 * v[i]; // By default, step 5% in each direction
    }

    Amoeba am(tol); // This contructs the instance of the Amoeba class

    vFinal = am.minimize(v, dv, y);
    
    printf("\nSuccess!  Mantel converged to a tolerance of %g in %d steps to a value "
	   "y = %16.8le with this vector:\n\n", 
	   tol, iteration, am.fmin);
    printf("[");
    for (i=0; i < Ndim; i++) 
	printf("%16.8le%s", vFinal[i], i != (Ndim-1) ? ", " : "]\n\n");

    fclose(flog);
    return(0);
}

//-----------------------------------------------------------------------------------

// The function which evaluates y(x1, ..., xn)
Doub y(VecDoub v) {
    Int i, tries=0, maxTries=10, varRead=0;
    Doub metric;
    char metricString[MAXSTR];
    FILE *fp;
    
    if (UsingScript) {
	system("rm mantel.out");
	fp = safeopen("mantel.out", "w");
	for (i=0; i < Ndim; i++) fprintf(fp, "%11.4le\n", v[i]);
	for (i=0; i < Ndim; i++) fprintf(stdout, "v[%d] = %11.4le\n", i, v[i]);
	fprintf(stderr, "Executing script '%s'\n", script);
	fclose(fp);
	system(script);
    }

    fprintf(flog, "%03d ", iteration);
    printf("[%03d] y(", iteration++);
    fflush(flog);
    for (i=0; i < Ndim; i++) {
	// printf("%11.4le%s", v[i], i == Ndim ? ") = ": ", ");
	printf("%11.4le%s", v[i], i == (Ndim-1) ? ") = ": " ");
	fprintf(flog, "%16.8le", v[i]);
    }
    
    if (UsingScript) {
	fp = safeopen("mantel.in", "r");
	if (fscanf(fp, "%le", &metric) != 1) {
	    printf("Error: contents of 'mantel.in' not a number.\n");
	    fprintf(stderr, "metric: %g\n", metric);
	    exit(1);
	}
	fclose(fp);
	printf("metric: %g\n", metric);
    } else {
	scanf("%s", metricString);
	if (metricString[0] == 'q') {
	    printf("\nExiting at user's request.\n");
	    fprintf(flog, "\nExiting at user's request.\n");
	    fclose(flog);
	    exit(0);
	}
	
	if (sscanf(metricString, "%le", &metric) != 1) 
	    do {
		printf("Error: I couldn't understand that as a number. Try again: ");
		scanf("%s", metricString);
		varRead = sscanf(metricString, "%le", &metric);
	    } while (tries++ < maxTries && varRead != 1);
	
	if (tries == maxTries) {
	    fprintf(stderr, "That's %d tries, for goodnes' sake.  I give up.\n", tries);
	    exit(1);
	}
    }

    fprintf(flog, "%16.8le\n", metric);

    return(latest = metric);
}
