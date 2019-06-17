/* ------------------------------------------------------------------ 
   This file contains a set of general library routines. 
   
   Tim Collins; collected first on 4/27/05. 
   ------------------------------------------------------------------
*/

#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

// For fileExists:
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h> // On msdos systems include <io.h> instead

#include "/home/tcol/utils/myConstants.h"
#include "/home/tcol/utils/mylib.h"

// Return pointer to a string containing the username, or NULL if not in the environment list
char *getUserName(char **env) {
    int i;

    for (i=0; env[i] != NULL; i++) 
	if (!strncmp(env[i], "USER", 4))
	    return(&(env[i][5]));
    return(NULL);
}

// Copy the first word of string s into w and return a pointer to
// the next word in s (skipping intermediate whitespace)
char *getword(char *w, char *s) {
    int i;

    for (i=0; !iswhitespace(w[i]=s[i]) && w[i] != '\0'; i++);
    w[i] = '\0';
    while (iswhitespace(s[i])) i++;
    return(&(s[i]));
}

// Reverse a string in place:
char *reverse(char *s) {
    char *temp, c;
    int i, I;

    I = strlen(s);
    temp = (char *)malloc(sizeof(char) * I);
    strcpy(temp, s);
    for (i=0; i < I; i++)
	s[i] = temp[I-i];
    free(temp);

    return(s);
}

// Get a line from fp; assumes the memory is already allocated:
char *getline(FILE *fp, char *l) {
    char *p, c;

    p=l;
    for (; (c=fgetc(fp)) != '\n' && c != '\0' && c != EOF; *(l++)=c);
    *l='\0';
    return(p);
}

// Get a line from fp; assumes the memory is already allocated; sets EOF flag if EOF is encountered
char *getlineEOF(FILE *fp, char *l, bool &eof) {
    char *p, c;

    p=l;
    for (; (c=fgetc(fp)) != '\n' && c != '\0' && c != EOF; *(l++)=c) 
	if (c == EOF) eof = true;
    *l='\0';
    return(p);
}

// Look through file "filename" for the first number following string "string"; return the final instance.
// the variable "found" is set to 1 if the variable is found, and 0 otherwise.
double fetch(char *filename, char *string, bool &found) {
    double x, val=0.0;
    FILE *fp;

    found = false;
    
    fp = fopen(filename, "r");
    if (fp == NULL) {
	fprintf(stderr, "Error opening file \"%s\" while fetching \"%s\": file does not (currently) exist\n", filename, string);
	return(0);
    } else {
	do {
	    queueup(fp, NULL, string);
	    if (!feof(fp))
		if (fscanf(fp, "%le", &x)==1) {
		    val = x;
		    found = true;
		}
	} while (!feof(fp));
	
	fclose(fp);
	
	return(val);
    }
}

// Reads an ascii tecplot data file and puts the data in a structure:
void readTecplot(char *s, struct tecplotData &d) {
    FILE *fp;

    fprintf(stderr, "Reading tecplot file \"%s\"\n", s);
    fp = safeopen(s, "r");
    readTecplotTitle(fp, d);
    readTecplotDatapackingFormat(fp, d);
    readTecplotVarNames(fp, d);
    readTecplotDimensions(fp, d);
    countTecplotZones(fp, d);
    readTecplotData(fp, d);
    fclose(fp);
}

void freeTecplot(struct tecplotData &d) {
    int i;
    
    switch (d.ndim) {
    case 1: d3free(d.d1, d.nzones, d.nvars); break;
    case 2: d4free(d.d2, d.nzones, d.nvars, d.I); break;
    case 3: d5free(d.d3, d.nzones, d.nvars, d.I, d.J); break;
    }

    free(d.title);
    for (i=0; i < d.nvars; i++) free(d.varNames[i]);
    free(d.varNames);
    free(d.filename);

    d.allocated = false;
}

// Allocate t to the size of s and copy s into t
char *alloacAndCp(char *t, char *s) {
    t = (char *)malloc(strlen(s) * sizeof(char));
    strcpy(t, s);
    return(t);
}

void readTecplotTitle(FILE *fp, struct tecplotData &d) {
    int i;
    char c;

    // First, find the title, if any
    rewind(fp);
    if (queueuplc(fp, NULL, "title")) {
	// fprintf(stderr, "Found the title string!\n");
	queueup(fp, NULL, "\"");
	// Now see how long the title is:
	for (i=0; fgetc(fp) != '\"'; i++);
	d.title = (char *)malloc(i * sizeof(char));
	rewind(fp);
	queueuplc(fp, NULL, "title");
	queueup(fp, NULL, "\"");
	for (i=0; (c=fgetc(fp)) != '\"'; i++) d.title[i] = c;
    } else {
	d.title = (char *)malloc(sizeof(char) * strlen("No title"));
	strcpy(d.title, "  No title");
    }
    fprintf(stderr, "  Tecplot title = \"%s\"\n", d.title);
}

void readTecplotDatapackingFormat(FILE *fp, struct tecplotData &d) {
    int i;
    char c;

    rewind(fp);
    if (queueuplc(fp, NULL, "datapacking")) {
	queueup(fp, NULL, "=");
	while (iswhitespace(c=fgetc(fp)));
	if (c=='P' || c=='p')
	    d.point = true;
	else
	    d.point = false;
    } else {
	d.point = true;
    }
    fprintf(stderr, "  Tecplot datapacking = \"%s\"\n", d.point ? "POINT" : "BLOCK");
}

// Read the variable names from an ascii Tecplot file. Assumes double quotes instead of single.
// Also assumes all variables are listed on one line. Which isn't what tecplot always does. You are warned.
void readTecplotVarNames(FILE *fp, struct tecplotData &d) {
    int nquotes, i, j, linelen;
    char c;
    char *tempstr;

    rewind(fp);
    if (queueuplc(fp, NULL, "variables")) {
	for (nquotes=0, linelen=0; (c=fgetc(fp)) != '\n'; linelen++) if (c=='\"') nquotes++;
	d.nvars = nquotes/2;
    } else d.nvars = 0;
    fprintf(stderr, "  %d variables:\n", d.nvars);
    rewind(fp);
    queueuplc(fp, NULL, "variables");
    d.varNames = (char **)malloc(sizeof(char *) * d.nvars);
    for (i=0; i < d.nvars; i++) {
	d.varNames[i] = (char *)malloc(sizeof(char) * (linelen+1));
	queueup(fp, NULL, "\"");
	for (j=0; (c=fgetc(fp)) != '\"'; j++) d.varNames[i][j] = c;
	d.varNames[i][j] = '\0';
    }
    for (i=0; i < d.nvars; i++)
	fprintf(stderr, "  \"%s\"\n", d.varNames[i]);
}

// Read the variable names from an ascii Tecplot file:
void readTecplotDimensions(FILE *fp, struct tecplotData &d) {
    char c;
    int i;
    
    rewind(fp);
    fprintf(stderr, "  Looking for dimenionality: \n"); 
    queueuplc(fp, NULL, "zone");
    // for (; (c=fgetc(fp))!='\n'; ) printf("[%c]\n",c);
    queueuplc(fp, NULL, "i");
    queueup(fp, NULL, "=");
    fscanf(fp, "%d", &(d.I));

    // fprintf(stderr, ">> I = %d\n", d.I);

    for (; c != '\n' && c != 'j' && c != 'J';) /* printf("[%c]\n", */ c=fgetc(fp) /*) */;
    // printf("==> [%c]\n", c);
    if (c == '\n') d.ndim = 1;
    else {
	d.ndim = 2;
	queueup(fp, NULL, "=");
	fscanf(fp, "%d", &(d.J));
	// printf("And j = %d\n", d.J);
	for (; c != '\n' && c != 'k' && c != 'K';) c=fgetc(fp);
	if (c != '\n') {
	    d.ndim = 3;
	    queueup(fp, NULL, "=");
	    fscanf(fp, "%d", &(d.K));
	    if (d.K == 1) d.ndim = 2;
	}
    }
    fprintf(stderr, "  %d dimensional data set of size %d", d.ndim, d.I);
    if (d.ndim>1) fprintf(stderr, "x%d", d.J);
    if (d.ndim>2) fprintf(stderr, "x%d", d.K);
    fprintf(stderr, "\n");
}

// Read the variable names from an ascii Tecplot file:
void countTecplotZones(FILE *fp, struct tecplotData &d) {
    rewind(fp);
    for (d.nzones=0; queueuplc(fp, NULL, "zone"); d.nzones++);
    fprintf(stderr, "  %d zone%c\n\n", d.nzones, d.nzones > 1 ? 's' : ' ');
}

// Read the variable names from an ascii Tecplot file:
// Currently this routine assumes, for BLOCK data packing,  that the first two 
// variables are vertex locations and the rest are cell-centered. If that's not 
// the case, the code will need to be modified.
void readTecplotData(FILE *fp, struct tecplotData &d) {
    int i, j, k, nz, nv, nprinttimes=12, imax, jmax, kmax, n;
    char c;

    rewind(fp);
    queueuplc(fp, NULL, "zone");
    if (queueuplc(fp, NULL, "solutiontime")) {
	d.timeData = true;
	d.times = dalloc(d.nzones);
	fprintf(stderr, "  File includes time-record data at times:\n");
    } else  d.timeData = false;
    rewind(fp);

    d.allocated = true;

    switch (d.ndim) {
    case 1:
	d.d1 = d3alloc(d.nzones, d.nvars, d.I);
	break;
    case 2:
	d.d2 = d4alloc(d.nzones, d.nvars, d.I, d.J);
	break;
    case 3:
	d.d3 = d5alloc(d.nzones, d.nvars, d.I, d.J, d.K);
	break;
    }
    
    if (d.point) {
	for (nz=0; nz < d.nzones; nz++) {
	    queueuplc(fp, NULL, "zone");
	    if (d.timeData) {
		queueuplc(fp, NULL, "solutiontime");
		queueup(fp, NULL, "=");
		fscanf(fp, "%le", &(d.times[nz]));
		fprintf(stderr, "  [%12.4le]%c", d.times[nz], (nz+1) % nprinttimes ? ' ' : '\n');
	    }
	    queueup(fp, NULL, "\n");
	    
	    switch (d.ndim) {
	    case 1:
		for (i=0; i < d.I; i++)
		    for (nv=0; nv < d.nvars; nv++)
			fscanf(fp, "%le", &(d.d1[nz][nv][i]));
		break;
	    case 2:
		for (j=0; j < d.J; j++)
		    for (i=0; i < d.I; i++)
			for (nv=0; nv < d.nvars; nv++)
			    fscanf(fp, "%le", &(d.d2[nz][nv][i][j]));
		break;
	    case 3:
		for (k=0; k < d.K; k++)
		    for (j=0; j < d.J; j++)
			for (i=0; i < d.I; i++)
			    for (nv=0; nv < d.nvars; nv++)
				fscanf(fp, "%le", &(d.d3[nz][nv][i][j][k]));
		break;
	    }
	}
    } else {
	for (nz=0; nz < d.nzones; nz++) {
	    queueuplc(fp, NULL, "zone");
	    queueuplc(fp, NULL, "datapacking");
	    while ((c=fgetc(fp))!='\n');
	    for (nv=0; nv < d.nvars; nv++) {
		if (nv < 2) { imax = d.I; jmax = d.J; kmax = d.K; } 
		else { imax = d.I-1; jmax = d.J-1; kmax = d.K-1; }
		switch (d.ndim) {
		case 1:
		    for (i=0; i < imax; i++)
			fscanf(fp, "%le", &(d.d1[nz][nv][i]));
		    break;
		case 2:
		    for (n=0, j=0; j < jmax; j++)
			for (i=0; i < imax; i++) {
			    fscanf(fp, "%le", &(d.d2[nz][nv][i][j]));
			    // fprintf(stderr, "[%10.2le]%c", d.d2[nz][nv][i][j], 
			    //		((n++)+1)%8==0 ? '\n' : ' ');
			}
		    break;
		case 3:
		    // fprintf(stderr, "(%d, %d, %d)\n", imax, jmax, kmax);
		    for (n=0, k=0; k < kmax; k++)
			for (j=0; j < jmax; j++)
			    for (i=0; i < imax; i++)
				fscanf(fp, "%le", &(d.d3[nz][nv][i][j][k]));
		    break;
		}
	    }
	}
    }
    if (d.timeData) fprintf(stderr, "\n");
}

// This reflects a 2D tecplot datastructure across the equator. Only works for 2D plot. Ignores time info,
// filename. Requires that the original datastructure persist. Doesn't copy polar coordinate info.
void tecplotHalf2FullSphere(struct tecplotData &d, struct tecplotData &e) {
    int nz, nv, i, j, j0;

    e.ndim   = d.ndim;
    e.nvars  = d.nvars;
    e.I      = d.I;
    e.J      = 1 + 2*(d.J-1);
    e.title  = d.title;
    e.nzones = d.nzones;
    
    e.varNames = (char **)malloc(sizeof(char *) * e.nvars);
    for (nv=0; nv < e.nvars; nv++)
	e.varNames[nv] = d.varNames[nv];

    e.d2 = d4alloc(e.nzones, e.nvars, e.I, e.J);

    for (nz=0; nz < e.nzones; nz++)
	for (nv=0; nv < e.nvars; nv++) 
	    for (i=0; i < e.I; i++)
		for (j=0; j < e.J; j++) {
		    if (j < d.J) {
			e.d2[nz][nv][i][j] = d.d2[nz][nv][i][j];
		    } else {
			j0 = e.J - (j+1);
			e.d2[nz][nv][i][j] = d.d2[nz][nv][i][j0];
			if (nv == 0) e.d2[nz][nv][i][j] = -e.d2[nz][nv][i][j]; // reflect the y coordinate
		    }
		}
}

// This expands a 2D tecplot data structure to 3D (symmetric in the 3rd dimension).
// kmax is the number of vertices the 3rd dimension.
// Assumes all variables after the first two are cell-centered.
void expandTo3D(tecplotData &d, int K) {
    int i, j, k, nv, nz, imax, jmax, kmax;
    char **oldVarNames;
    char *newVar = "zl [um]";
    double x, y, z, r, theta, phi, x0, y0;

    if (d.ndim != 2) {
	fprintf(stderr, "Error in expandTo3D: Expecting a 2D tecplot datastructure\n");
	exit(1);
    }

    d.ndim      = d.ndim+1;
    d.nvars     = d.nvars+1;
    d.K         = K;
    oldVarNames = d.varNames;
    d.varNames  = (char **)malloc(sizeof(char *) * d.nvars);
    for (nv=0; nv < 2; nv++) 
	d.varNames[nv] = oldVarNames[nv];
    d.varNames[2] = (char *)malloc(sizeof(char) * (1 + sizeof(newVar)));
    strcpy(d.varNames[2], newVar);
    for (nv=2; nv < d.nvars-1; nv++)
	d.varNames[nv+1] = oldVarNames[nv];

    fprintf(stderr, "  %d variables now\n", d.nvars);
    for (nv=0; nv < d.nvars; nv++)
	fprintf(stderr, "  variable %d: %s\n", nv, d.varNames[nv]);

    fprintf(stderr, "  3rd dimension of size: %d\n", d.K);

    d.d3 = d5alloc(d.nzones, d.nvars, d.I, d.J, d.K);

    for (nz=0; nz < d.nzones; nz++) {
	for (nv=3; nv < d.nvars; nv++) {
	    imax = d.I-1; jmax = d.J-1; kmax = d.K-1; 
	    
	    for (i=0; i < imax; i++)
		for (j=0; j < jmax; j++)
		    for (k=0; k < kmax; k++)
			d.d3[nz][nv][i][j][k] = d.d2[nz][nv-1][i][j];
	}
	imax = d.I; jmax = d.J; kmax = d.K; 
	for (i=0; i < imax; i++)
	    for (j=0; j < jmax; j++)
		for (k=0; k < kmax; k++) {
		    x0    = d.d2[nz][0][i][j];
		    y0    = d.d2[nz][1][i][j];
		    r     = sqrt(x0*x0 + y0*y0);
		    theta = atan2(y0, x0);
		    phi   = ((double)k) * 2.0 * M_PI / ((double)(kmax-1));
		    y     = r * sin(theta) * cos(phi); // Draco uses and unusual definition of x, y, z; x is the axis of rotation
		    z     = r * sin(theta) * sin(phi);
		    x     = r * cos(theta);
		    d.d3[nz][0][i][j][k] = x;
		    d.d3[nz][1][i][j][k] = y;
		    d.d3[nz][2][i][j][k] = z;
		}
    }
}

void interpTecplotToRegularMesh(struct tecplotData &d, struct tecplotData &e, 
					      double Rmin, double Rmax, int I, int J) {
    int i, j;
    static const int X=0, Y=1, RHO=2;
    double ***cc;
    double ***F; 
    double **G;
    double r, theta;
    int i0, j0, i1, j1, i2, i0last;
    int nz = 0; // zone number
    bool found;
    double x1, y1, x2, y2, x3, y3, x4, y4, f1, f2, f3, f4;
    char *s = "Regular Mesh";
    bool verbose=true;

    cc = d3alloc(2, d.I-1, d.J-1);
    for (i=0; i < d.I-1; i++)
	for (j=0; j < d.J-1; j++) {
	    cc[X][i][j] = 0.25 * (d.d2[nz][X][i][j] + d.d2[nz][X][i+1][j] + d.d2[nz][X][i][j+1] + d.d2[nz][X][i+1][j+1]);
	    cc[Y][i][j] = 0.25 * (d.d2[nz][Y][i][j] + d.d2[nz][Y][i+1][j] + d.d2[nz][Y][i][j+1] + d.d2[nz][Y][i+1][j+1]);
	}

    // Now find the values of rho at the vertices. This is assuming that Rmin and Rmax are large/small enough to 
    // restrict the range in radius so that you aren't interpolating below the cell center for i=0 (first zone)
    // or i=d.I-1 (last zone). I.e., I am not handling the lower and upper i boundaries, just the j boundaries. 
    G = d2alloc(d.I, d.J);

    for (i=1; i < d.I-1; i++)
	for (j=0; j < d.J; j++) {
	    if (j == 0) {
		x1 =  cc[X][i-1][j];   y1 = -cc[Y][i-1][j];   f1 = d.d2[nz][RHO][i-1][j];   // modified for bc
		x2 =  cc[X][i][j];     y2 = -cc[Y][i][j];     f2 = d.d2[nz][RHO][i][j];     // modified for bc
		x3 =  cc[X][i][j];     y3 =  cc[Y][i][j];     f3 = d.d2[nz][RHO][i][j]; 
		x4 =  cc[X][i-1][j];   y4 =  cc[Y][i-1][j];   f4 = d.d2[nz][RHO][i-1][j]; 
	    } else if (j == d.J-1) {
		x1 =  cc[X][i-1][j-1]; y1 =  cc[Y][i-1][j-1]; f1 = d.d2[nz][RHO][i-1][j-1]; 
		x2 =  cc[X][i][j-1];   y2 =  cc[Y][i][j-1];   f2 = d.d2[nz][RHO][i][j-1]; 
		x3 = -cc[X][i][j-1];   y3 =  cc[Y][i][j-1];   f3 = d.d2[nz][RHO][i][j-1];   // modified for bc
		x4 = -cc[X][i-1][j-1]; y4 =  cc[Y][i-1][j-1]; f4 = d.d2[nz][RHO][i-1][j-1]; // modified for bc 
	    } else {
		x1 =  cc[X][i-1][j-1]; y1 =  cc[Y][i-1][j-1]; f1 = d.d2[nz][RHO][i-1][j-1]; 
		x2 =  cc[X][i][j-1];   y2 =  cc[Y][i][j-1];   f2 = d.d2[nz][RHO][i][j-1]; 
		x3 =  cc[X][i][j];     y3 =  cc[Y][i][j];     f3 = d.d2[nz][RHO][i][j]; 
		x4 =  cc[X][i-1][j];   y4 =  cc[Y][i-1][j];   f4 = d.d2[nz][RHO][i-1][j]; 
	    }
	    
	    G[i][j] = fourPtInterp(x1, y1, f1, x2, y2, f2, x3, y3, f3, x4, y4, f4, 
				   d.d2[nz][X][i][j], d.d2[nz][Y][i][j]);
	}

    // Allocate and populate the fine, regular mesh:
    e.d2 = d4alloc(1, 3, I, J);
    for (i0last=0, i=0; i < I; i++) {
	fprintf(stderr, "*");
	if ((i+1)%80==0) fprintf(stderr, "\n");
	for (j=0; j < J; j++) {
	    r     = Rmin + (((double)i)/((double)(I-1))) * (Rmax-Rmin);
	    theta = (((double)j)/((double)(J-1))) * M_PI / 2.0; // assuming a half sphere
	    e.d2[nz][X][i][j] = r * cos(theta);
	    e.d2[nz][Y][i][j] = r * sin(theta);
	    // Next determine the zone of the original mesh which contains this point, via exhaustive search (can be improved)
	    for (i0=1, found=false; i0 < d.I-1 && !found; i0++) {
		i2 = 1 + ((i0-1 + i0last))%(d.I-2);
		for (j0=0; j0 < d.J-1 && !found; j0++)
		    // Search each zone in d2; d.d2 contains the vertex locations
		    // to see which contains e.d2, the fine-mesh vertices
		    if (inQuadrilateral(d.d2[nz][X][i2][j0],     d.d2[nz][Y][i2][j0], 
					d.d2[nz][X][i2+1][j0],   d.d2[nz][Y][i2+1][j0], 
					d.d2[nz][X][i2+1][j0+1], d.d2[nz][Y][i2+1][j0+1], 
					d.d2[nz][X][i2][j0+1],   d.d2[nz][Y][i2][j0+1], 
					e.d2[nz][X][i][j],       e.d2[nz][Y][i][j], 
					verbose && false)) {
			// G[i][j] are the values of rho at the low-res mesh vertices. Use these to interpolate to the
			// fine-mesh vertices.
			e.d2[nz][RHO][i][j] = fourPtInterp(d.d2[nz][X][i2][j0],     d.d2[nz][Y][i2][j0],     G[i2][j0], 
							   d.d2[nz][X][i2+1][j0],   d.d2[nz][Y][i2+1][j0],   G[i2+1][j0], 
							   d.d2[nz][X][i2+1][j0+1], d.d2[nz][Y][i2+1][j0+1], G[i2+1][j0+1], 
							   d.d2[nz][X][i2][j0+1],   d.d2[nz][Y][i2][j0+1],   G[i2][j0+1], 
							   e.d2[nz][X][i][j],       e.d2[nz][Y][i][j]); 
			// e.d2[nz][RHO][i][j] = G[i2][j0];
			found = true;
		    }
	    }
	    i0last = i2-1;
	    if (verbose && false) fprintf(stderr, "## <%d, %d> <%d, %d> (%12.4le, %12.4le) (%12.4le, %12.4le) { %12.4le } %s\n",
				 i, j, i0, j0, r, theta * 180.0/M_PI, e.d2[nz][X][i][j], e.d2[nz][Y][i][j], e.d2[nz][RHO][i][j], 
				 found ? "found" : "lost"); 
	}
    }

    fprintf(stderr, "\n");
    e.ndim     = 2;
    e.title    = (char *)malloc(sizeof(char) * (1 + strlen(s)));
    strcpy(e.title, s);
    e.I        = I;
    e.J        = J;
    e.nvars    = 3; // x, y, rho
    e.nzones   = 1;
    e.varNames = (char **)malloc(sizeof(char *) * 3);
    for (i=0; i < 3; i++)
	e.varNames[i] = d.varNames[i];
    e.point    = d.point;
    e.timeData = false;

    // Release memory for the temporary arrays:
    d3free(cc, 2, d.I-1);
    d2free(G, d.I);
}

// Compute the r, theta values for a 2D tecplot data structure
void findTecplotPolarCoords(struct tecplotData &d) {
    int i, j, k, nz;

    d.r0     = dalloc(d.I); // r for j==0
    d.theta0 = dalloc(d.J); // theta for i==0
    if (d.ndim > 2)
	d.phi0   = dalloc(d.K); // phi for j != 0 (say for d.J/2)
    else
	d.K  = 1;

    for (nz = 0; nz < d.nzones; nz++)
	for (k=0; k < d.K; k++)
	    for (j=0; j < d.J; j++)
		for (i=0; i < d.I; i++) {
		    if (d.ndim == 2) {
			if (j == 0) d.r0[i] = sqrt(d.d2[nz][0][i][j] * d.d2[nz][0][i][j] + d.d2[nz][1][i][j] * d.d2[nz][1][i][j]);
			if (i == 0) d.theta0[j] = atan2(d.d2[nz][1][i][j], d.d2[nz][0][i][j]);
			if (d.theta0[j] < 0.0) d.theta0[j] = d.theta0[j] + 2.0 * M_PI;
		    } else if (d.ndim > 2) {
			if (j == 0) d.r0[i] = sqrt(d.d3[nz][0][i][j][k] * d.d3[nz][0][i][j][k] + d.d3[nz][1][i][j][k] * d.d3[nz][1][i][j][k]);
			if (i == 0) d.theta0[j] = atan2(d.d3[nz][1][i][j][k], d.d3[nz][0][i][j][k]);
			if (d.theta0[j] < 0.0) d.theta0[j] = d.theta0[j] + 2.0 * M_PI;
			if (i == 0 && j == d.J/2) d.phi0[k] = atan2(d.d3[nz][2][i][j][k], d.d3[nz][1][i][j][k]);
			if (d.phi0[k] < 0.0) d.phi0[k] = d.phi0[k] + 2.0 * M_PI;
		    }
		}

    if (d.ndim > 2 && false)
	for (nz=0, j=d.J/2, i=0, k=0; k < d.K; k++)
	    fprintf(stderr, "&& k=%d; phi0 = %12.4le (z=%12.4le y = %12.4le) \n", 
		    k, d.phi0[k] * 180.0/M_PI, d.d3[nz][2][i][j][k], d.d3[nz][1][i][j][k]);
}

// Write a tecplot dataset to file:
void writeTecplot(char *s, struct tecplotData &d) {
    int i, j, k, nz, nv, n, imax, jmax, kmax, nperline=5, m;
    FILE *fp;
    bool point = false;
    double tiny = 1.0e-32;
    
    fp = safeopen(s, "w");
    fprintf(fp, "title = \"%s\"\n", d.title);
    fprintf(fp, "variables = ");
    for (nv=0; nv < d.nvars; nv++)
	fprintf(fp, "\"%s\" ", d.varNames[nv]);
    fprintf(fp, "\n");
    for (nz=0; nz < d.nzones; nz++) {
	fprintf(fp, "zone i=%d ", d.I);
	if (d.ndim > 1) fprintf(fp, "j=%d ", d.J);
	if (d.ndim > 2) fprintf(fp, "k=%d ", d.K);
	fprintf(fp, "\n");
	if (point) {
	    fprintf(fp, "datapacking=point\n"); 
	    switch (d.ndim) {
	    case 1:
		for (i=0; i < d.I; i++)
		    for (nv=0; nv < d.nvars; nv++) 
			fprintf(fp, "%12.4le%c", fabs(d.d1[nz][nv][i]) < tiny ? 0.0 : d.d1[nz][nv][i], 
				nv==d.nvars-1 ? '\n' : ' ');
		break;
	    case 2:
		for (j=0; j < d.J; j++)
		    for (i=0; i < d.I; i++)
			for (nv=0; nv < d.nvars; nv++)
			    fprintf(fp, "%12.4le%c", d.d2[nz][nv][i][j], 
				    nv==d.nvars-1 ? '\n' : ' ');
		break;
	    case 3:
 		for (k=0; k < d.K; k++)
		    for (j=0; j < d.J; j++)
			for (i=0; i < d.I; i++)
			    for (nv=0; nv < d.nvars; nv++)
				fprintf(fp, "%12.4le%c", fabs(d.d3[nz][nv][i][j][k]) < tiny ? 0.0 : d.d3[nz][nv][i][j][k], 
					nv==d.nvars-1 ? '\n' : ' ');
		break;
	    }
	} else {
	    m=0;
	    fprintf(fp, "datapacking=block\n"); 
	    if (d.nvars > 3)
		fprintf(fp, "varlocation=([3-%d]=cellcentered)\n", d.nvars); 
	    else
		fprintf(fp, "varlocation=([3]=cellcentered)\n"); 
	    for (nv=0; nv < d.nvars; nv++) {
		// for (m=0, nv=0; nv < 1; nv++) {
		if (nv < d.ndim) { imax = d.I; jmax = d.J; kmax = d.K; } 
		else { imax = d.I-1; jmax = d.J-1; kmax = d.K-1; }
		switch (d.ndim) {
		case 1:
		    for (i=0, n=0; i < imax; i++)
			fprintf(fp, "%12.4le%c", d.d1[nz][nv][i], ((n++)+1)%nperline==0 ? '\n' : ' ');
		    break;
		case 2:
		    for (j=0, n=0; j < jmax; j++)
			for (i=0; i < imax; i++) {
			    m++;
			    fprintf(fp, "%12.4le%c", d.d2[nz][nv][i][j], ((n++)+1)%nperline==0 ? '\n' : ' ');
			}
		    if (((n++))%nperline!=0) fprintf(fp, "\n");
		    break;
		case 3:
		    for (k=0, n=0; k < kmax; k++)
			for (j=0; j < jmax; j++)
			    for (i=0; i < imax; i++) {
				m++;
				fprintf(fp, "%12.4le%c", d.d3[nz][nv][i][j][k], ((n++)+1)%nperline==0 ? '\n' : ' ');
			    }
		    if (((n++))%nperline!=0) fprintf(fp, "\n");
		    break;
		}
	    }
	    fprintf(fp, "\n");
	    fprintf(stderr, "  %d numbers written to Tecplot file \"%s\"\n", (unsigned int)m, s);
	}
    }
    fclose(fp);
}

// Write a tecplot dataset to file:
void writeTecplot3DIslice(char *s, struct tecplotData &d, int i0) {
    int i, j, k, nz, nv, n, imax, jmax, kmax, nperline=5, m;
    FILE *fp;
    bool point = false;
    double tiny = 1.0e-32;
    
    fp = safeopen(s, "w");

    fprintf(fp, "title = \"%s\"\n", d.title);
    fprintf(fp, "variables = ");
    for (nv=0; nv < d.nvars; nv++)
	fprintf(fp, "\"%s\" ", d.varNames[nv]);
    fprintf(fp, "\n");
    for (nz=0; nz < d.nzones; nz++) {
	fprintf(fp, "zone i=%d j=%d k=%d\n", 1, d.J, d.K);
	if (point) {
	    fprintf(fp, "datapacking=point\n"); 
	    for (k=0; k < d.K; k++)
		for (j=0; j < d.J; j++) {
		    i = i0;
		    // for (i=0; i < d.I; i++)
		    for (nv=0; nv < d.nvars; nv++)
			fprintf(fp, "%12.4le%c", fabs(d.d3[nz][nv][i][j][k]) < tiny ? 0.0 : d.d3[nz][nv][i][j][k], 
				nv==d.nvars-1 ? '\n' : ' ');
		}
	} else {
	    m=0;
	    fprintf(fp, "datapacking=block\n"); 
	    if (d.nvars > 4)
		fprintf(fp, "varlocation=([4-%d]=cellcentered)\n", d.nvars); 
	    else
		fprintf(fp, "varlocation=([4]=cellcentered)\n"); 
	    for (nv=0; nv < d.nvars; nv++) {
		// for (m=0, nv=0; nv < 1; nv++) {
		if (nv < d.ndim) { imax = d.I; jmax = d.J; kmax = d.K; } 
		else { imax = d.I-1; jmax = d.J-1; kmax = d.K-1; }
		for (k=0, n=0; k < kmax; k++)
		    for (j=0; j < jmax; j++) {
			// for (i=0; i < imax; i++) {
			i = i0;
			m++;
			fprintf(fp, "%12.4le%c", d.d3[nz][nv][i][j][k], ((n++)+1)%nperline==0 ? '\n' : ' ');
		    }
		if (((n++))%nperline!=0) fprintf(fp, "\n");
	    }
	    fprintf(fp, "\n");
	    fprintf(stderr, "  %d numbers written\n", (unsigned int)m);
	}
    }
    fclose(fp);
}

// Given a string, returns the index of that variable in a tecplot data structure
int findTecplotVarIndex(char *s, struct tecplotData &d) {
    int i, var=-1;

    for (i=0; i < d.nvars; i++)
	if (!strcmp(s, d.varNames[i])) { var = i; i = d.nvars; }
    if (var == -1) fprintf(stderr, "  Variable \"%s\" not present\n", s);
    else fprintf(stderr, "  Variable \"%s\" has index %d\n", s, var);
    return(var);
}

// Finds the max abs value of a list of numbers
double maxAbsList(double *v, int N) {
    double max;
    int i;

    max = fabs(v[0]);
    for (i=0; i < N; i++)
	if (fabs(v[i]) > max) max = fabs(v[i]);
    return(max);
}

// Uses tanh functions to window so it's 1 within [-xmax, xmax], 0 outside,
// with a smooth taper.
double windowFnSymm(double x, double xmax) {
    return(windowFn(x, -xmax, xmax));
}

// Uses tanh functions to window so it's 1 within [-xmax, xmax], 0 outside,
// with a smooth taper.
double windowFn(double x, double xmin, double xmax) {
    return(0.5*(1.0+tanh(xmax - x)) * 0.5*(1.0+tanh(x - xmin)));
}

// Copies a single line from one stream to another. 
// Returns false if end of file has been reached.
bool passLine(FILE *fin, FILE *fout) {
    char c;

    // Interestingly, if you combine these loops and try to direct putc to
    // send the character to NULL, it seg faults. So I special case
    // fout == NULL.
    if (fout==NULL) 
	do c = fgetc(fin), fout; 
	while (c != EOF && c != '\n');
    else
	do putc(c = fgetc(fin), fout); 
	while (c != EOF && c != '\n');
    return(!(c == EOF));
}

/* Given the port angles thetap and phip, return the polar angle of 
   point theta, phi with respect to the port angles. See euler.nb
   for the derivation. Used by detune.cc. */
double beam_theta(double theta, double phi, double thetap, double phip) {
    double x, y, z, r;
    
    double thetanew;

    x = cos(phi) * cos(phip) * sin(theta) -
	sin(phi) * sin(phip) * sin(theta);
    y = cos(thetap) * (cos(phip) * sin(phi) * sin(theta) + cos(phi) * sin(phip) * sin(theta)) -
	cos(theta) * sin(thetap);
    z = cos(theta) * cos(thetap) + 
	sin(thetap) * (cos(phip) * sin(phi) * sin(theta) + cos(phi) * sin(phip) * sin(theta));

    r = sqrt(x*x + y*y);
    thetanew = atan2(r, z);
    return(thetanew);
}



// Write a tecplot file of a spherical projection 2-D array
// If thetaOnly, plot illumination as a function of theta (averaged over phi).
// If append, append to the plot file.

void plotSphericalProjection(struct sphericalProjection &S, char *s, bool thetaOnly, bool append) {
    FILE *fp;
    int i, j, k, ii, jj;
    int nphi, ntheta, nn, mm, nskip, mskip;
    double x, y, z;
    double phiSum, phiAve;

    if (append)
	fp = safeopen(s, "a");
    else
	fp = safeopen(s, "w");

    if (thetaOnly) {
	// Header
	if (!append) 
	    fprintf(fp, "variables=\"theta (deg)\",\"I\"\n");

	fprintf(fp, "zone i=%d f=point t=\"Phi-averaged Spherical Projection\"\n", S.N);
	
	// Find the average over phi:
	for (i=0; i < S.N; i++) { // for each theta...
	    for (phiSum = 0.0, j=0; j < S.M; j++) phiSum += S.f[i][j];
	    phiAve = phiSum / ((double)S.M);
	    // fprintf(stderr, "%d %g %g %g\n", i, S.theta[i][0] * rad, phiAve, phiSum);
	    fprintf(fp, "%8.6le %8.6le\n", S.theta[i][0] * rad, phiAve);
	}
    } else {
	// Header
	if (!append)
	    fprintf(fp, "variables=\"x\",\"y\",\"z\",\"I\"\n");
	
	if (false) {
	    // Aitoff projection
	    fprintf(fp, "zone i=%d j=%d f=point t=\"Aitoff Projection\"\n", S.N, S.M);
	    for (j=0; j < S.M; j++) {
		for (i=0; i < S.N; i++) {
		    aitoff(S.theta[i][j], S.phi[i][j], &x, &y);
		    fprintf(fp, "%8.6e %8.6e %8.6e %8.6e\n", x, y, 0.0, S.f[i][j]);
		}
	    }
	    
	    nskip = 10; 
	    mskip = 10;
	    nn = S.N / nskip;
	    mm = S.M / mskip;
	    
	    // Aitoff grid
	    fprintf(fp, "zone i=%d j=%d f=point t=\"Aitoff grid\"\n", nn+1, mm+1);
	    for (jj=0; jj <= mm; jj++) {
		for (ii=0; ii <= nn; ii++) {
		    i = ii * nskip;
		    j = jj * mskip;
		    if (i >= S.N) i = S.N-1;
		    if (j >= S.M) j = S.M-1;
		    aitoff(S.theta[i][j], S.phi[i][j], &x, &y);
		    fprintf(fp, "%8.6e %8.6e %8.6e %8.6e\n", x, y, 0.0, S.f[i][j]);
		}
	    }
	}
	
	// 3-D Spherical projection
	// fprintf(fp, " VARLOCATION=([6,7,8]=nodal)");
	fprintf(fp, "\nzone i=%d j=%d f=point t=\"Spherical Projection\"\n", S.N, S.M);
	for (j=0; j < S.M; j++) {
	    for (i=0; i < S.N; i++) {
		x = sin(S.theta[i][j]) * cos(S.phi[i][j]);
		y = sin(S.theta[i][j]) * sin(S.phi[i][j]);
		z = cos(S.theta[i][j]);
		fprintf(fp, "%8.6e %8.6e %8.6e %8.6e\n", x, y, z, S.f[i][j]);
	    }
	}
    }
    fclose(fp);
}

// Calculate the spherical projection generated by a port configuration. 
// If spotMask is true in pcf the sim-input file is used to determine spot-masking properties.
// If repointedAngles use the repointed angles rather than the unrepointed.
// If pointingPlot show the port locations rather than the true projection.
// If kruer, multiply not just by the cosine factor but by cosine^3, which lore says 
// Kruer's book advises as a way to account for the effects of the atmosphere.
void calcSphericalProjection(struct pointingFile &pf, struct sphericalProjection &S, bool kruer, bool useEmods) {
    int i, j, l;
    int ring, port;
    int N, M;
    double stdev;
    double thetap, phip;
    double hemisphereSign;
    double dp;
    double xp, yp, zp;
    double xi = 1.0;
    double dS;

    N = S.N; M = S.M;

    for (i=0; i < N; i++) /* loop over theta */
	for (j=0; j < M; j++) { /* loop over phi */
	    S.theta[i][j] = M_PI * ((double)(i))/((double)(N-1));
	    S.phi[i][j]   = M_PI * ((double)(j))/((double)(M-1)) * 2.0;
	    S.f[i][j]     = 0.0;
	    
	    for (ring=0; ring < pf.nrings; ring++) { /* loop over rings */
		for (l=0; l < pf.nports[ring]; l++) {
		    port = pf.rn[ring][l] - 1; // 0-referenced

		    // Note: We are using the unrepointed port angles because the repointing angles
		    // are captured by the offsets represented by pf.xcent and pf.ycent, which must
		    // be pre-calculated by calcRepointingOffsets(pf). That the caller is reponsible
		    // for making sure that routine is called first.
		    thetap = pf.theta0[port] * deg;
		    phip   = pf.phi0[port]   * deg;

		    if (thetap < M_PI / 2.0) hemisphereSign = -1.0; // the secondary ellipse is always toward the equator
		    else hemisphereSign = 1.0;
		    // This is the dot product between the port unit vector and the point unit vector:
		    dp = cos(S.theta[i][j]) * cos(thetap) + 
		 	   cos(S.phi[i][j]) * cos(phip) * sin(S.theta[i][j]) * sin(thetap) + 
	  		   sin(S.phi[i][j]) * sin(phip) * sin(S.theta[i][j]) * sin(thetap);
		    // if (dp >= 0.0 && port == plotPort) { // If the point isn't over the horizon
		    if (dp >= 0.0) { // If the point isn't over the horizon
			unitVector(xp, yp, zp, S.theta[i][j], S.phi[i][j]); // <xp, yp, zp> points toward the point of interest
			// Now rotate so the FF is in the x-y plane:
			// Rotate by phi to bring the plane down to the proper port angle:
			rotPhi  (xp, yp, zp, -phip);
			// Rotate around y by the port angle thetap:
			rotyaxis(xp, yp, zp, -thetap);
			// Now the far-field coordinates are x, y
			dS = generalSpotShape(xp * hemisphereSign, yp, 
						      pf.p[ring][PCF_SGexp], pf.p[ring][PCF_SG2nd], pf.SGSM, // SGexp, 
						                                                             // SGexp2nd, SGSM
						      pf.p[ring][PCF_eta],   pf.p[ring][PCF_eta2nd],         // eta, eta2nd
						      pf.p[ring][PCF_R0]*xi, pf.RSM,                         // R0, RSM
						      pf.xcent[port],        pf.ycent[port],                 // (xcent, ycent)
						   // 0.0,                   0.0,                            // xcent, ycent
						      pf.p[ring][PCF_xoff],                                  // xoff
						      pf.xoffSM,             // xoffSM (a neg offset is equator-ward)
						      pf.p[ring][PCF_f2nd])  // f2nd
			    * dp // Projection effect: the cosine of the angle between the unit vectors
			    * (kruer ? dp * dp : 1);  // According to Kruer the cosine factor should be cosine^3 to 
                                                      // account for the atmosphere
			if (useEmods) dS *= pf.emods[ring];
			S.f[i][j] += dS / pf.spotIntegral[port]; // Normalize so each spot, no matter how distorted, has the same integral
		    }
		}
	    }
	}
}

// Find the illumination uniformity for a given port configuration with file name s. If spotMask is true, the 
// sim-input file is used to determine spot-masking properties.
double illuminUnif(struct pointingFile &pf, struct sphericalProjection &S, bool thetaOnly) {
    double stdev;

    normalizeSphericalProjection(S);                   // normalize the results 
    stdev = sigmarmsSphericalProjection(S, thetaOnly); // find the rms

    return(stdev);
}

// This is John Marozas' general spot-shape formula for a combinary primary and secondary elliptical super-Gaussian
// spot, with circular super-Gaussian spot-masking
double generalSpotShape(double x, double y, 
			double SGexp, double SGexp2nd, double SGSM, 
			double eta, double eta2nd, 
			double R0, double RSM, 
			double xcent, double ycent, 
			double xoff, double xoffSM, 
			double f2nd) {
    // x and y are the far-field coordinates (um)
    // SGexp and SGexp2nd are the super-Gaussian exponents for the primary spot and secondary ellipse
    // SGSM is the super-Gaussian exponent of the spot-masking function
    // eta and eta2nd are the ellipticities for the primary spot and secondary ellipse
    // R0 is the primary spot radius
    // RSM is the spot-mask aperturing radius (um)
    // (xcent, ycent) is the center of the primary spot 
    // xoff is the secondary-ellipse x offset (x being the "vertical" direction)
    // xoffSM is the offset of the spot-masking ellipse (in the vertical direction)
    // f2nd is the secondary-ellipse weight (the total intensity being propto Iprimary + f2 * I2ndary)
    //
    // The offsets are in units of the target radius

    double spotShapeR(double x, double y, double eta, double xoff, double xcent, double ycent);
    double IspotMask, IprimaryEllipse, IsecondaryEllipse;
    static const double log5pct = log(0.05);

    IprimaryEllipse   = exp(log5pct * pow(spotShapeR(x, y, eta,    0.0,    xcent, ycent) / R0,  SGexp));
    IsecondaryEllipse = exp(log5pct * pow(spotShapeR(x, y, eta2nd, xoff,   xcent, ycent) / R0,  SGexp2nd));
    IspotMask         = exp(log5pct * pow(spotShapeR(x, y, 1.0,    xoffSM, 0.0,   0.0  ) / RSM, SGSM));

    return(IspotMask * (IprimaryEllipse + f2nd * IsecondaryEllipse));
}

double spotShapeR(double x, double y, double eta, double xoff, double xcent, double ycent) {
    return(sqrt(pow(eta * (x - xoff - xcent), 2.0) + pow(y - ycent, 2.0)));
}

// Find the unit vector with spherical-polar angles theta and phi (rad):
void unitVector(double &x, double &y, double &z, double theta, double phi) {
    x = cos(phi) * sin(theta);
    y = sin(phi) * sin(theta);
    z = cos(theta);
}

// Rotate vector <x, y, z> about the z axis by angle phi (rad) in a ccw direction i.e. toward the y axis
void rotPhi(double &x, double &y, double &z, double phi) {
    double x1, y1; 
    x1 = x * cos(phi) - y * sin(phi);
    y1 = y * cos(phi) + x * sin(phi);
    x = x1; y = y1;
}

// Rotate vector <x, y, z> about the x axis by angle theta (rad):
void rotTheta(double &x, double &y, double &z, double theta) {
    double y1, z1; 
    y1 = y * cos(theta) + z * sin(theta);
    z1 = z * cos(theta) - y * sin(theta);
    y = y1; z = z1;
}

// Rotate vector <x, y, z> about the y axis by angle theta (rad):
void rotyaxis(double &x, double &y, double &z, double theta) {
    double x1, z1; 
    x1 = x * cos(theta) + z * sin(theta);
    z1 = z * cos(theta) - x * sin(theta);
    x = x1; z = z1;
}

// Related to 2-D spherical projection arrays and 3-D plan projection arrays:

// Allocate a plane projection object and and initialize its contents to zero:
void allocPlaneProjection(int N, int M, struct planeProjection &P) {
    int i, j;

    P.N = N;
    P.M = M;
    P.x = d2alloc(N, M);
    P.y = d2alloc(N, M);
    P.z = d2alloc(N, M);
    P.f = d2alloc(N, M);
    for (i=0; i < N; i++)
	for (j=0; j < M; j++)
	    P.x[i][j] = P.y[i][j] = P.z[i][j] = P.f[i][j] = 0.0;
}

// Allocate and initialize to zero a spherical projection structure:
void allocSphericalProjection(int N, int M, struct sphericalProjection &S) {
    int i, j;

    S.N = N;
    S.M = M;
    S.theta = d2alloc(N, M);
    S.phi   = d2alloc(N, M);
    S.f     = d2alloc(N, M);
    for (i=0; i < N; i++)
	for (j=0; j < M; j++)
	    S.theta[i][j] = S.phi[i][j] = S.f[i][j] = 0.0;
}

// deallocate a sphericalProjection structure
void sphericalProjectionFree(struct sphericalProjection &S) {
    int i, j;

    d2free(S.theta, S.N);
    d2free(S.phi,   S.N);
    d2free(S.f,     S.N);
    // free(&S);
}

double averageSphericalProjection(struct sphericalProjection &S) {
    double Omega;  /* total solid angle */
    double sum;    /* used to normalize the output */
    double dOmega; /* differential solid angle */
    int i, j;

    for (sum=0.0, Omega=0.0, i=1; i < S.N; i++) /* loop over theta */
	for (j=1; j < S.M; j++) { /* loop over phi */
	    dOmega = sin(S.theta[i][j]) 
		* (S.theta[i][j] - S.theta[i-1][j])
		* (S.phi[i][j] - S.phi[i][j-1]);
	    sum += 0.25 * (S.f[i][j] + S.f[i-1][j] + S.f[i][j-1] + S.f[i-1][j-1])
		* dOmega;
	    Omega += dOmega;
	}
    
    sum /= Omega;
    
    return(sum);
}

// Finds a normalized (i.e. fractional) sigma RMS of spherical projection S
double sigmarmsSphericalProjection(struct sphericalProjection &S, bool thetaOnly) {
    double Omega;   // total solid angle 
    double sum;     // used to normalize the output 
    double dOmega;  // differential solid angle 
    double thisp; 
    double ave; 
    double sumdPhi;
    double dphi;
    double dtheta;
    double sumdTheta;
    int i, j, ip, im, jp, jm;
    static int  arraySize = 0;
    static bool allocated = false;
    static double *phiAve; // S.f averaged over phi

    if (thetaOnly) {
	// A 1-D temporary array is used and kept around for future calls.
	// If the current call needs a larger array, it is reallocated.
	if (S.N > arraySize) {
	    arraySize = S.N;
	    if (allocated) free(phiAve);
	    phiAve = dalloc(arraySize);
	    allocated = true;
	}
	
	ave       = 0.0;
	sumdTheta = 0.0;
	// Find the average over phi:
	for (i=0; i < S.N; i++) { // for each theta...
	    phiAve[i] = 0.0;
	    sumdPhi   = 0.0;
	    for (j=0; j < S.M; j++) { // loop over phi 
		if (j == 0) 
		    dphi       = S.phi[i][j+1] - S.phi[i][j];
		else 
		    dphi       = S.phi[i][j] - S.phi[i][j-1];
		phiAve[i] += S.f[i][j] * dphi;
		sumdPhi   += dphi;
	    }
	    phiAve[i] /= sumdPhi;

	    j = 0;
	    if (i == 0)
		dtheta     = S.theta[i+1][j] - S.theta[i][j];
	    else
		dtheta     = S.theta[i][j] - S.theta[i-1][j];
	    sumdTheta += dtheta * sin(S.theta[i][j]);
	    ave       += phiAve[i] * dtheta * sin(S.theta[i][j]);
	}
	ave /= sumdTheta; // This is the ave over theta of phiAve

	// Now find the rms of avePhi wrt theta:
	sumdTheta = 0.0;
	sum       = 0.0;
	j         = 0;
	for (i=0; i < S.N; i++) { // for each theta...
	    if (i == 0)
		dtheta     = S.theta[i+1][j] - S.theta[i][j];
	    else
		dtheta     = S.theta[i][j] - S.theta[i-1][j];
	    sumdTheta += dtheta * sin(S.theta[i][j]);
	    sum       += dtheta * sin(S.theta[i][j]) * (phiAve[i] / ave - 1.0) * (phiAve[i] / ave - 1.0);
	}
	sum /= sumdTheta;

	// free(phiAve);
    } else {
	ave = averageSphericalProjection(S);

	for (sum=0.0, Omega=0.0, i=1; i < S.N; i++) /* loop over theta */
	    for (j=1; j < S.M; j++) { /* loop over phi */
		dOmega = sin(S.theta[i][j]) 
		    * (S.theta[i][j] - S.theta[i-1][j])
		    * (S.phi[i][j] - S.phi[i][j-1]);
		thisp = 0.25 * (S.f[i][j] + S.f[i-1][j] + S.f[i][j-1] + S.f[i-1][j-1]);
		sum = sum + (thisp/ave - 1.0) * (thisp/ave - 1.0) * dOmega;
		Omega += dOmega;
	    }
	sum /= Omega;
    }
    
    return(sqrt(sum));
}

// Normalize the array S.f[i][j], defined over the surface of a unit sphere, and return the stdev. 
void normalizeSphericalProjection(struct sphericalProjection &S) {
    double ave; 
    int i, j;

    ave = averageSphericalProjection(S);

    for (i=0; i < S.N; i++) // loop over theta 
	for (j=0; j < S.M; j++)  // loop over phi 
	    S.f[i][j] /= ave;
}

/*
// Normalize the array S.f[i][j], defined over the surface of a unit sphere, and return the stdev. 
double normalizeSphericalProjection(struct sphericalProjection &S, bool thetaOnly) {
    double Omega;  // total solid angle 
    double sum;    // used to normalize the output 
    double dOmega; // differential solid angle 
    double stdevsq;// square of the standard deviation of the perturbation 
    double phiSum; // sum along a latitude line (for thetaOnly == true)
    double phiOmegaSum; // solid angle of a phi row of cells
    double phiAve; // altitudinal average
    int i, j;

    // First, calculate the surface integral:
    for (sum=0.0, Omega=0.0, i=1; i < S.N; i++) // loop over theta 
	for (j=1; j < S.M; j++) {               // loop over phi   
	    dOmega = sin(S.theta[i][j]) 
		* (S.theta[i][j] - S.theta[i-1][j])
		* (S.phi[i][j] - S.phi[i][j-1]);
	    sum += 0.25 * (S.f[i][j] + S.f[i-1][j] + S.f[i][j-1] + S.f[i-1][j-1])
		* dOmega;
	    Omega += dOmega;
	}
    
    sum /= Omega;
    
    // Normalize the function defined on the surface of the sphere, S.f[i][j]:
    if (fabs(sum) > 1.0e-12)
	for (i=0; i < S.N; i++)     // loop over theta 
	    for (j=0; j < S.M; j++) // loop over phi   
		S.f[i][j] = (S.f[i][j] / sum);

    // Now what's the standard deviation?
    if (!thetaOnly) {
	for (stdevsq=0.0, i=1; i < S.N; i++) // loop over theta 
	    for (j=1; j < S.M; j++) {      // loop over phi   
		dOmega = sin(S.theta[i][j]) 
		    * (S.theta[i][j] - S.theta[i-1][j])
		    * (S.phi[i][j]   - S.phi[i][j-1]);
		stdevsq += pow(S.f[i][j] - 1.0, 2.0) * dOmega; 
	    }
    } else {
	for (stdevsq=0.0, i=1; i < S.N; i++) { // loop over theta 
	    for (j=1; j < S.M; j++) {      // loop over phi   
		phiSum = 0.0;
		phiOmegaSum = 0.0;
		dOmega = sin(S.theta[i][j]) 
		    * (S.theta[i][j] - S.theta[i-1][j])
		    * (S.phi[i][j]   - S.phi[i][j-1]);
		phiSum += S.f[i][j];
		phiOmegaSum + dOmega;
	    }
	    phiAve = phiSum / ((double)S.M);
	    stdevsq += pow(phiSum - 1.0, 2.0) * phiOmegaSum; 
	}
    }
    
    stdevsq /= Omega;

    return(sqrt(stdevsq));
}
*/

// Here we calculate a matric which maximizes the spread in overlapped wavelengths (given in an
// accompnaying pointing file structure), with a windowing width dLambdaMaxIR (A).
// S is only used for its theta and phi grid.
// The metric function "Lambda" is put in S.f[i][j].
// Note, this metric is non-negative. 
double detuningMetric(pointingFile &pf, struct sphericalProjection &S, double dLambdaMaxIR, bool thetaOnly) {
    int i, j, k, l;
    int ring, port, portIndex;
    int N, M;
    double Omega;  // total solid angle 
    double sum;    // used to normalize the output 
    double dOmega; // differential solid angle 
    double thetap, phip;
    double dtheta, dphi;
    double hemisphereSign;
    double dp;
    double xp, yp, zp;
    double metric;
    struct sphericalProjection *Q;
    bool verbose = true;

    Q = (struct sphericalProjection *)malloc(pf.totalNumPorts * sizeof(struct sphericalProjection));
    for (i=0; i < pf.totalNumPorts; i++) 
	allocSphericalProjection(S.N, S.M, Q[i]);
    
    N = S.N; M = S.M;

    // fprintf(stderr, "Got here!\n"); exit(0);

    if (verbose) fprintf(stderr, "Calculating %d spherical projections, one for each beam...\n", pf.totalNumPorts);
    // Create a spherical projection for each beam (port)
    for (i=0; i < N; i++) /* loop over theta */
	for (j=0; j < M; j++) { /* loop over phi */
	    S.theta[i][j] = M_PI * ((double)(i))/((double)(N-1));
	    S.phi[i][j]   = M_PI * ((double)(j))/((double)(M-1)) * 2.0;
	    
	    for (portIndex=0, ring=0; ring < pf.nrings; ring++) { /* loop over rings */
		for (l=0; l < pf.nports[ring]; l++, portIndex++) {
		    port = pf.rn[ring][l] - 1; // 0-referenced

		    // Note: We are using the unrepointed port angles because the repointing angles
		    // are captured by the offsets represented by pf.xcent and pf.ycent, which must
		    // be pre-calculated by calcRepointingOffsets(pf). That the caller is reponsible
		    // for making sure that routine is called first.
		    thetap = pf.theta0[port] * deg;
		    phip   = pf.phi0[port]   * deg;

		    if (thetap < M_PI / 2.0) hemisphereSign = -1.0; // the secondary ellipse is always toward the equator
		    else hemisphereSign = 1.0;
		    // This is the dot product between the port unit vector and the point unit vector:
		    dp = cos(S.theta[i][j]) * cos(thetap) + 
		 	   cos(S.phi[i][j]) * cos(phip) * sin(S.theta[i][j]) * sin(thetap) + 
	  		   sin(S.phi[i][j]) * sin(phip) * sin(S.theta[i][j]) * sin(thetap);
		    // if (dp >= 0.0 && port == plotPort) { // If the point isn't over the horizon
		    if (dp >= 0.0) { // If the point isn't over the horizon
			unitVector(xp, yp, zp, S.theta[i][j], S.phi[i][j]); // <xp, yp, zp> points toward the point of interest
			// Now rotate so the FF is in the x-y plane:
			// Rotate by phi to bring the plane down to the proper port angle:
			rotPhi  (xp, yp, zp, -phip);
			// Rotate around y by the port angle thetap:
			rotyaxis(xp, yp, zp, -thetap);
			// Now the far-field coordinates are x, y
			Q[portIndex].f[i][j] = generalSpotShape(xp * hemisphereSign, yp, 
						      pf.p[ring][PCF_SGexp], pf.p[ring][PCF_SG2nd], pf.SGSM, // SGexp, SGexp2nd, SGSM
						      pf.p[ring][PCF_eta],   pf.p[ring][PCF_eta2nd],         // eta, eta2nd
						      pf.p[ring][PCF_R0],    pf.RSM,                         // R0, RSM
						      pf.xcent[port],        pf.ycent[port],                 // (xcent, ycent)
						   // 0.0,                   0.0,                            // xcent, ycent
						      pf.p[ring][PCF_xoff],                                  // xoff
						      pf.xoffSM,             // xoffSM (a neg offset is equator-ward)
						      pf.p[ring][PCF_f2nd])  // f2nd
			    * dp // Projection effect: the cosine of the angle between the unit vectors
			    * dp * dp // According to Kruer the cosine factor should be cosine^3 to account for the atmosphere
			    / pf.spotIntegral[port]; // Normalize so each spot, no matter how distorted, has the same integral
		    }
		}
	    }
	}

    // Normalize the projections for each beam (the phase plates guarantee this):
    for (i=0; i < pf.totalNumPorts; i++) 
	normalizeSphericalProjection(Q[i]); 

    if (verbose) fprintf(stderr, "Calculating detuned overlap function Lambda...\n");
    for (sum=0.0, Omega=0.0, i=0; i < N; i++) // loop over theta 
	for (j=0; j < M; j++) {               // loop over phi   
	    S.f[i][j] = 0.0;
	    for (k=0; k < pf.totalNumPorts; k++) 
		for (l=0; l < k; l++) 
		    S.f[i][j] += 
			Q[k].f[i][j] * // intensity for port k
			Q[l].f[i][j] * // intensity for port l
			windowFnSymm(pf.dLambdaIR[k], dLambdaMaxIR) * // window fn to prevent shifts outside allowed limit
			windowFnSymm(pf.dLambdaIR[l], dLambdaMaxIR) * // window fn to prevent shifts outside allowed limit
			fabs(pf.dLambdaIR[k] - pf.dLambdaIR[l]); // ...and weighting by the difference in wavelengths!
	}

    if (verbose) fprintf(stderr, "Calculating metric...\n");
    // Now sum Lambda over the sphere:
    for (metric=0.0, Omega=0.0, i=0; i < N; i++) 
	for (j=0; j < M; j++) {               
	    if (i > 0) 
		dtheta = S.theta[i][j] - S.theta[i-1][j];
	    else
		dtheta = S.theta[i+1][j] - S.theta[i][j];
	    if (j > 0)
		S.phi[i][j] - S.phi[i][j-1];
	    else
		S.phi[i][j+1] - S.phi[i][j];
	    dOmega = sin(S.theta[i][j]) * dtheta * dphi;
	    metric += S.f[i][j] * dOmega;
	    Omega += dOmega;
	}
    metric /= Omega;

    // Free up memory!
    for (i=0; i < pf.totalNumPorts; i++) sphericalProjectionFree(Q[i]);
    free(Q);

    return(metric);
}

/* This calculates the coordinates using an Aitoff projection. 
   Assumes theta in [0, pi], phi in [0, 2pi].  The resulting 
   (x, y) coordinates have x in [-1, 1], y in [-0.5, 0.5]. */
void aitoff(double theta, double phi, double *x, double *y) {
    double a;

    theta -= M_PI_2;
    phi   -= M_PI;

    if (fabs(theta) < 1.0e-12) {
	*x = 0.0;
	*y = 1.0;
    } else if (fabs(theta - M_PI) < 1.0e-12) {
	*x = 0.0;
	*y = -1.0;
    } else {
	a = acos(cos(theta) * cos(phi / 2.0));
	*x = (2.0 * a * cos(theta) * sin(phi / 2.0) / sin(a)) / M_PI;
	*y = (a * sin(theta) / sin(a)) / M_PI;
    }
}

// Write a draco PD pointing configuration file
void write_pcf(char *fname, struct pointingFile &p) {
    int i, j;
    FILE *fp;
    time_t rawtime;
    struct tm *timeinfo;
    
    fp = safeopen(fname, "w");

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    fprintf(fp, "! Pointing configuration file generated on %s", 
	    asctime(timeinfo));
    
    fprintf(fp, "%s%d\n", p.header, 11 + p.nrings);
    for (i=0; i < 11 + p.nrings; i++) 
	fprintf(fp, "%d ", p.block_sizes[i]);
    fprintf(fp, "\n\n");

    // The port angles are in deg.

    for (i=0; i < p.block_sizes[0]; i++)
	fprintf(fp, "%9.3lf%c", windowShift(p.theta0[i], 0.0, 360.0), (((i+1)%4)==0) ? '\n' : ' ');
    fprintf(fp, "\n\n");

    for (i=0; i < p.block_sizes[1]; i++)
	fprintf(fp, "%9.3lf%c", windowShift(p.phi0[i], 0.0, 360.0), (((i+1)%4)==0) ? '\n' : ' ');
    fprintf(fp, "\n\n");

    for (i=0; i < p.block_sizes[2]; i++)
	fprintf(fp, "%9.3lf%c", windowShift(p.theta1[i], 0.0, 360.0), (((i+1)%4)==0) ? '\n' : ' ');
    fprintf(fp, "\n\n");

    for (i=0; i < p.block_sizes[3]; i++)
	fprintf(fp, "%9.3lf%c", windowShift(p.phi1[i], 0.0, 360.0), (((i+1)%4)==0) ? '\n' : ' ');
    fprintf(fp, "\n\n");


    for (j=0; j < p.nrings; j++) {
	for (i=0; i < p.nports[j]; i++)
	    fprintf(fp, "%3d%c", p.rn[j][i], (((i+1)%8)==0) ? '\n' : ' ');
	fprintf(fp, "\n\n");
    }

    fprintf(fp, "\n\n");
    
    for (i=0; i < 7; i++) {
	for (j=0; j < p.nrings; j++) 
	    fprintf(fp, "%7.2lf ", p.p[j][i]);
	fprintf(fp, "\n");
    }

    fclose(fp);
}

// shifts x to within [x0, x1) through integratl multiple shifts of x1-x0.
double windowShift(double x, double x0, double x1) {
    double y, z;

    if (fabs(x1-x0) < 1.0e-32) return(x);
    if (x1 < x0) {
	z = x0; x0 = x1; x1 = z;
    }
    
    if (x >= 0.0) 
	return(modf((x-x0)/(x1-x0), &z) * (x1-x0));
    else {
	z = x/(x1-x0);
	y = z-floor(z);
	return(y*(x1-x0));
    }
}    

// Calculate the offsets in x and y corresponding to the repointings in this pointing configuration.
// Used in spherical projections of the pointing configuration.
void calcRepointingOffsets(struct pointingFile &p) {
    int port;
    double x1, y1, z1; // unit vectors for the old and new pointings
    double hemisphereSign;

    for (port=0; port < p.totalNumPorts; port++) {
	if (p.theta0[port] < 90.0) hemisphereSign = -1.0; // In the N hemisphere positive x is north; in the S it's south
	else hemisphereSign = 1.0;
	// To calculate the offsets, rotate the FF plane until it is coincident with x-y, then
	// the offsets are x1-x0 and y1-y0. 
	// unitVector(x0, y0, z0, p.theta0[i], p.phi0[i]); // <x0, y0, z0> points toward the port, fyi
	unitVector(x1, y1, z1, p.theta1[port]*deg, p.phi1[port]*deg); // <x1, y1, z1> points toward the repointed location on target
	// Now rotate so the FF is in the x-y plane--i.e., the rotations which would take <x0, y0, z0> so x0 == y0 == 0.
	// Rotate by phi0 to bring the plane down to the proper port angle:
	rotPhi  (x1, y1, z1, -p.phi0[port]*deg);
	// Rotate around y by the port angle theta0:
	rotyaxis(x1, y1, z1, -p.theta0[port]*deg);
	p.xcent[port] = x1 * hemisphereSign; 
	p.ycent[port] = y1;
    }
}

// This reads in the pointing file. It also writes out the header.
// Note: This code assumes the standard five-ring PD logical structure.
void read_pcf(char *fname, struct pointingFile &p, bool SpotMask) {
    FILE *fp;
    int i, j, port, ring, l, I;
    char c;
    bool present;
    bool verbose = false;

    fp = safeopen(fname, "r");

    // Read the header
    i=0;
    while ((c = fgetc(fp))=='!') {
	p.header[i++] = c;
	while ((c = fgetc(fp))!='\n') p.header[i++] = c;
	p.header[i++] = c;
    }
    p.header[i++] = '\0';
    
    if (c != '\n') ungetc(c, fp);
    
    fscanf(nc(fp), "%d", &I);
    p.nrings = I - 11; // The spot-shape block is 7 lines; there are four blocks of theta, phi
    // pointing angles; so the remaining numbers in this list are the numbers of ports in each
    // ring, and the number of rings is 17 - this number.

    // Actually, the number of rings is same as the last seven entries in the block sizes array,
    // whereas p.nrings-11 may be off by one, depending on whether there's an entry before the 
    // spot properties block indicating how the spot radius is being specified.

    for (i=0; i < I; i++) 
	fscanf(nc(fp), "%d", &(p.block_sizes[i]));

    p.nrings = p.block_sizes[I-1];

    if (verbose) fprintf(stderr, "%d blocks and %d rings\n", i, p.nrings);

    if (verbose) fprintf(stderr, "Block sizes: %d %d %d %d\n", 
	    p.block_sizes[0], p.block_sizes[1], p.block_sizes[2], p.block_sizes[3]);
    
    for (i=0; i < p.block_sizes[0]; i++)
	fscanf(nc(fp), "%lf", &(p.theta0[i]));
    for (i=0; i < p.block_sizes[1]; i++)
	fscanf(nc(fp), "%lf", &(p.phi0[i]));
    for (i=0; i < p.block_sizes[2]; i++) 
	fscanf(nc(fp), "%lf", &(p.theta1[i]));
    for (i=0; i < p.block_sizes[3]; i++) 
	fscanf(nc(fp), "%lf", &(p.phi1[i]));

    p.totalNumPorts = p.block_sizes[0];

    for (j=0; j < p.nrings; j++) p.nports[j] = p.block_sizes[4+j];

    // This reads the ports in each ring--those are rn[ring][n]
    for (j=0; j < p.nrings; j++)
	for (i=0; i < p.nports[j]; i++) {
	    fscanf(nc(fp), "%d", &(p.rn[j][i]));
	    p.ringsByPort[p.rn[j][i]-1] = j; // record for ea. port the ring it's in (0-referenced)
	}

    if (verbose) 
	for (i=0; i < I; i++)
	    fprintf(stderr, "block_sizes[%d] = %d\n", i, p.block_sizes[i]);

    // True if there is a block prior to the spot spec saying how the radii are being specified:
    if (p.block_sizes[I-8] == 1) p.radiusInterpPresent = true; else p.radiusInterpPresent = false;
    if (verbose) fprintf(stderr, "block_sizes[I-8] = %d\n", p.block_sizes[I-8]);
    // Note this can break if the number of ports in the last ring is one!
    // Hypothetically one could add the total number of numbers in all the blocks and see if it's
    // present that way.

    if (p.radiusInterpPresent) {
	fscanf(nc(fp), "%d", &(p.radiusInterp));
	if (verbose) fprintf(stderr, "Radius interpretation switch present in PCF with a value of %d\n",
	    p.radiusInterp);
    }

    // Now read the spot properties.  So e.g. the SG exponents are p.p[ring][0].
    for (i=0; i < 7; i++)
	for (j=0; j < p.nrings; j++)
	    fscanf(nc(fp), "%lf", &(p.p[j][i]));

    if (verbose) {
	fprintf(stderr, "First ring: SGexp = %g, radius = %g, ellipticity = %g\n", 
			 p.p[0][0], p.p[0][1], p.p[0][2]);
	fprintf(stderr, "Last ring: SGexp = %g, radius = %g, ellipticity = %g\n", 
			 p.p[p.nrings-1][0], p.p[p.nrings-1][1], p.p[p.nrings-1][2]);
    }

    fclose(fp);

    // Spot Masking parameters:
    p.RSM      = 1000.0;
    p.xoffSM   = 0.0;
    p.SGSM     = 20.0;
    p.spotMask = SpotMask;
    if (SpotMask && (fp = fopen("the_simulation_input.txt", "r"))!=NULL) {
	// fprintf(stderr, "Searching sim-input file for SM parms\n");
	// open the sim-input file and find any spot-masking parameters
	present = queueupFileStatus(fp, NULL, "FFspotmask_R95"); 
	if (present) {
	    queueuplc(fp, NULL, "=");
	    fscanf(fp, "%lf", &(p.RSM));
	}
	rewind(fp);
	present = queueupFileStatus(fp, NULL, "FFspotmask_SG"); 
	if (present) {
	    queueuplc(fp, NULL, "=");
	    fscanf(fp, "%lf", &(p.SGSM));
	}
	fclose(fp);
    } 

    // NOTE! If the repointed angles get changed, call calcRepointingOffsets and 

    // Calculate the centroids of the primary and secondary ellipses
    // due to the repointing angles:
    calcRepointingOffsets(p);

    // Calculate the integral of each spot in the far-field for use in normalization:
    calcPCFSpotIntegrals(p);

    // Initialize the detuning values:
    for (i=0; i < p.totalNumPorts; i++)
	p.dLambdaIR[i] = 0.0; 
}

// Read a ring emod file and find the last emod record and put it in the pointing file
// structure.
void getEmods(char *filename, pointingFile &pf) {
    FILE *fp;
    int i, ntimes, nrings, time, ring;
    double d;

    fp = safeopen(filename, "r");

    fscanf(fp, "%d %d", &ntimes, &nrings);
    for (i=0; i < ntimes; i++) fscanf(fp, "%lf", &d); // read and discard times
    for (time=0; time < ntimes; time++)
	for (ring=0; ring < nrings; ring++)
	    fscanf(fp, "%lf", &(pf.emods[ring]));

    fclose(fp);
}

void calcPCFSpotIntegrals(struct pointingFile &p) {
    int port, ring, l;

    for (ring=0; ring < p.nrings; ring++) /* loop over rings */
	for (l=0; l < p.nports[ring]; l++) { // loop over ports in each ring
	    port = p.rn[ring][l] - 1; // 0-referenced
	    // fprintf(stderr, "port %02d\n", port);
	    p.spotIntegral[port] = spotIntegral(p.p[ring][PCF_SGexp], p.p[ring][PCF_SG2nd], p.SGSM, 
						p.p[ring][PCF_eta],   p.p[ring][PCF_eta2nd],
						p.p[ring][PCF_R0],    p.RSM,                
						p.xcent[port],        p.ycent[port],        
						p.p[ring][PCF_xoff],  
						p.xoffSM,             
						p.p[ring][PCF_f2nd]);
	}
}

// This performs an integral of the general spot function. 
double spotIntegral(double SGexp, double SGexp2nd, double SGSM, 
		    double eta, double eta2nd, 
		    double R0, double RSM, 
		    double xcent, double ycent, 
		    double xoff, double xoffSM, 
		    double f2nd) {
    double tolerance = 1.0e-7;
    double I=0.0, Ilast=0.0;
    double x, y, dx, dy;
    double limit = 1.5;
    int i, j;
    int N=10;

    do {
	N *= 2;
	dx = 2.0 * limit / ((double) N); // area element size
	dy = dx;
	Ilast = I;
	I = 0.0;
	for (i=0; i < N; i++)
	    for (j=0; j < N; j++) {
		x = (((double)i) + 0.5) * dx - limit;
		y = (((double)j) + 0.5) * dy - limit;
		I += dx * dy * generalSpotShape(x, y, 
						SGexp, SGexp2nd, SGSM, 
						eta, eta2nd, 
						R0, RSM, 
						xcent, ycent, 
						xoff, xoffSM, 
						f2nd);
	    }
	// fprintf(stderr, "[%05d] %16.8le %16.8le %16.8le\n", N, I, Ilast, fabs(I-Ilast)/(I+Ilast));
    } while (fabs(I-Ilast)/(I+Ilast) > tolerance);
    return(I);
}

double distsq(double x, double y, double z) {
    return(x*x+y*y+z*z);
}

double dist(double x, double y, double z) {
    double distsq(double x, double y, double z);

    return(sqrt(distsq(x, y, z)));
}

// Streams and files:

// This function takes a stream, skips any white space or comments
// (signified with a '!') and returns the stream.  "nc" for
// "no comments".
FILE *nc(FILE *f) {
    char c;
    int not_done = 1;
    int iswhitespace(char c);
    
    do {
	do c = fgetc(f);
	while (c != EOF && iswhitespace(c));
	// Comment encountered; read to the end of the line
	if (c=='!') 
	    do c = fgetc(f);
	    while (c != EOF && c != '\n');
	if (c == EOF || !iswhitespace(c)) not_done = 0;
    } while (not_done);
    if (!iswhitespace(c) && c != EOF) ungetc(c, f);
    return(f);
}

/* open a file and quit out on an error */

FILE *safeopen(char *s, char *f)
{
  FILE *temp;

  /* fprintf(stderr, "Safeopening file [%s] with mode [%s].\n", s, f); */
  
  if ((temp = fopen(s, f))==NULL) {
    fprintf(stderr, "Error opening file \"%s\"\n", s);
    exit(0);
  }
  return(temp);
}

/* open a file and quit out on an error; include traceback info */

FILE *safeopentb(char *s, char *f, char *tb)
{
  FILE *temp;

  if ((temp = fopen(s, f))==NULL) {
    fprintf(stderr, "Error opening file \"%s\" in routine {%s}\n", s, tb);
    exit(0);
  }
  return(temp);
}

// Skip white space and return the new char *
char *sskipws(char *s) {
    int iswhitespace(char c);
    while(iswhitespace(*s)) s++;
    return(s);
}

// This function does 2-D linear interpolation on a rectangular grid.
double linint2D(int nx, int ny, double *x, double *y, double **f, 
		double x0, double y0) {

    double linint(double y1, double y2, double x1, double x2, double x);

    int i, j;
    int i0, j0;

    double f0, f1;
    
    // Note that f[][] is taken as the value of the function at the
    // grid points, not at the cell centers.

    // First, find the indices of the grid cell containing the point.

    if (x0 <= x[0]) i0 = 0;
    else if (x0 >= x[nx-1]) i0 = nx-2;
    else for (i=0; i < nx-1; i++)
	if (x0 >= x[i] && x0 <= x[i+1]) i0 = i;

    if (y0 <= y[0]) j0 = 0;
    else if (y0 >= y[ny-1]) j0 = ny-2;
    else for (j=0; j < ny-1; j++)
	if (y0 >= y[j] && y0 <= y[j+1]) j0 = j;

    i0 = (i0 < 0 ? 0 : (i0 >= nx-1 ? nx-1 : i0));
    j0 = (j0 < 0 ? 0 : (j0 >= ny-1 ? ny-1 : j0));

//    fprintf(stderr, "%12.4le %12.4le %12.4le %12.4le\n", 
//	    f[i0][j0], f[i0+1][j0], f[i0][j0+1], f[i0+1][j0+1]);

    f0 = linint(f[i0][j0],   f[i0+1][j0],   x[i0], x[i0+1], x0);
    f1 = linint(f[i0][j0+1], f[i0+1][j0+1], x[i0], x[i0+1], x0);
    /*
    fprintf(stderr, "[%3d %3d|%12.4le < %12.4le < %12.4le|%12.4le < "
	    "%12.4le < %12.4le|%12.4le, %12.4le, %12.4le, %12.4le: "
	    "%12.4le]\n", 
	    i0, j0, 
	    x[i0], x0, x[i0+1], 
	    y[j0], y0, y[j0+1], 
	    f[i0][j0], f[i0+1][j0], f[i0][j0+1], f[i0+1][j0+1], 
	    linint(f0, f1, y[j0], y[j0+1], y0));
    */
    return(linint(f0, f1, y[j0], y[j0+1], y0));
}

// This function does 2-D linear interpolation on a rectangular grid.
double logint2D(int nx, int ny, double *x, double *y, double **f, 
		double x0, double y0) {

    double logint(double y1, double y2, double x1, double x2, double x);

    int i, j;
    int i0, j0;

    double f0, f1;
    
    // Note that f[][] is taken as the value of the function at the
    // grid points, not at the cell centers.

    // First, find the indices of the grid cell containing the point.

    if (x0 <= x[0]) i0 = 0;
    else if (x0 >= x[nx-1]) i0 = nx-2;
    else for (i=0; i < nx-1; i++)
	if (x0 >= x[i] && x0 <= x[i+1]) i0 = i;

    if (y0 <= y[0]) j0 = 0;
    else if (y0 >= y[ny-1]) j0 = ny-2;
    else for (j=0; j < ny-1; j++)
	if (y0 >= y[j] && y0 <= y[j+1]) j0 = j;

    f0 = logint(f[i0][j0], f[i0+1][j0], x[i0], x[i0+1], x0);
    f1 = logint(f[i0][j0+1], f[i0+1][j0+1], x[i0], x[i0+1], x0);

    /*
    fprintf(stderr, "[%3d %3d|%12.4le < %12.4le < %12.4le|%12.4le < "
	    "%12.4le < %12.4le|%12.4le, %12.4le, %12.4le, %12.4le: "
	    "%12.4le]\n", 
	    i0, j0, 
	    x[i0], x0, x[i0+1], 
	    y[j0], y0, y[j0+1], 
	    f[i0][j0], f[i0+1][j0], f[i0][j0+1], f[i0+1][j0+1], 
	    logint(f0, f1, y[j0], y[j0+1], y0));
    */
    return(logint(f0, f1, y[j0], y[j0+1], y0));
}

// Wait for a file to be created
void waitForFileCreation(char *s) {
    time_t rawtime;
    struct tm *timeinfo;

    // Wait for the file to be created
    for (; !fileExists(s); ) {
	fprintf(stderr, "Waiting for %s to be created...", s);
	system("sleep 60");

	// Get start-time info:
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	
	fprintf(stderr, "%s", asctime(timeinfo));
    }
}

// This function tests for file existence.
int fileExists(char *f) {
    int file_descriptor;
    int exists;

    file_descriptor = open(f, O_EXCL);
    
    if (file_descriptor < 0) exists=0; else exists=1;

    close(file_descriptor);

    return(exists);
}

// This function replaces the extension of the string f with
// string n, unless the file name doesn't have an extension, in which
// case it just appends the new extension. In the latter case it won't
// append the new extension if it means the string length will be
// greater than N.
char *newExt(char *f, char *n, int N) {
  char c;
  int i, L, hasExt=0, extIndex;

  L = strlen(f);

  for (i=0; i < L; i++)
    if (f[i]=='.') { hasExt=1; extIndex=i; }

  if (hasExt) {
    // Is there enough room to add the new extension?
    if (extIndex+strlen(n)>N) return(f);
    f[extIndex]='\0';
  } else {
    // Again, is there enough room?
    if (L+strlen(n)>N) return(f);
  }

  sprintf(f, "%s.%s", f, n);
  return(f);
}

// Plots an array in ascii tecplot format
void tecplot3d(char *filename, char *title, 
	       char *xname, char *yname, char *zname, char *fname, 
	       double *x, double *y, double *z, double ***q, 
	       int xN, int yN, int zN) {
    int i, j, k;
    FILE *fp;
    FILE *safeopen(char *s, char *f);
    char s[1000];
    
    fp = safeopen(filename, "w");

    fprintf(fp, "title=\"%s\"\nvariables=\"%s\",\"%s\",\"%s\"\n",
	    title, xname, yname, fname);
    
    for (k=0; k < zN; k++) {
	sprintf(s, "%s = %13.5le", zname, z[k]);

	fprintf(fp, "zone i=%d j=%d f=point T=\"%s\"\n", yN, xN, s);
	
	for (i=0; i < xN; i++)
	    for (j=0; j < yN; j++)
		fprintf(fp, "%16.8le %16.8le %16.8le\n",
			x[i], y[j], q[i][j][k]);
    }

    fclose(fp);
}

// Plots an array in ascii tecplot format
void tecplot2d(char *filename, char *title, 
	       char *xname, char *yname, char *fname, 
	       double *x, double *y, double **q, 
	       int xN, int yN) {
    int i, j;
    FILE *fp;
    FILE *safeopen(char *s, char *f);
    
    fp = safeopen(filename, "w");

    fprintf(fp, "title=\"%s\"\nvariables=\"%s\",\"%s\",\"%s\"\n",
	    title, xname, yname, fname);
    fprintf(fp, "zone i=%d j=%d f=point\n", yN, xN);

    for (i=0; i < xN; i++)
	for (j=0; j < yN; j++)
	    fprintf(fp, "%16.8le %16.8le %16.8le\n",
		    x[i], y[j], q[i][j]);

    fclose(fp);
}

// Plots an array in ascii tecplot format
void tecplot1d(char *filename, char *title, 
	       char *xname, char *fname, 
	       double *x, double *q, 
	       int xN) {
    int i;
    FILE *fp;
    FILE *safeopen(char *s, char *f);
    
    fp = safeopen(filename, "w");

    fprintf(fp, "title=\"%s\"\nvariables=\"%s\",\"%s\"\n",
	    title, xname, fname);
    fprintf(fp, "zone i=%d f=point\n", xN);

    for (i=0; i < xN; i++)
	fprintf(fp, "%16.8le %16.8le\n",
		x[i], q[i]);

    fclose(fp);
}

/* returns the sign of a number: -1 if x < 0 or 1 if x >= 0 */
double sign(double x) {
    if (fabs(x) < 1.0e-20)
	return(0.0);
    else
	return(x / fabs(x));
}

/* Free an array of doubles. */
void d2free(double **p, int n) {
    int i;

    for (i=0; i < n; i++)
	free(p[i]);

    free(p);
}

void d3free(double ***p, int n, int m) {
    int i;

    for (i=0; i < n; i++)
	d2free(p[i], m);

    free(p);
}

void d4free(double ****p, int n, int m, int o) {
    int i;

    for (i=0; i < n; i++)
	d3free(p[i], m, o);

    free(p);
}

void d5free(double *****p, int n1, int n2, int n3, int n4) {
    int i;

    for (i=0; i < n1; i++)
	d4free(p[i], n2, n3, n4);

    free(p);
}

double ***d3alloc(int n, int m, int o) {
    void testptr(char *ptr);
    double **d2alloc(int n, int m);
    double ***p;
    int i;

    p = (double ***)malloc(sizeof(double **) * n);

    testptr((char *)p);
    
    for (i=0; i < n; i++)
	p[i] = d2alloc(m, o);

    return(p);
}

double ****d4alloc(int n, int m, int o, int l) {
    void testptr(char *ptr);
    double ***d3alloc(int n, int m, int o);
    double ****p;
    int i;

    p = (double ****)malloc(sizeof(double ***) * n);

    testptr((char *)p);
    
    for (i=0; i < n; i++)
	p[i] = d3alloc(m, o, l);

    return(p);
}

double *****d5alloc(int n1, int n2, int n3, int n4, int n5) {
    void testptr(char *ptr);
    double ****d4alloc(int n, int m, int o, int q);
    double *****p;
    int i;

    p = (double *****)malloc(sizeof(double ****) * n1);

    testptr((char *)p);
    
    for (i=0; i < n1; i++)
	p[i] = d4alloc(n2, n3, n4, n5);

    return(p);
}

int ***i3alloc(int n, int m, int o) {
    void testptr(char *ptr);
    int **i2alloc(int n, int m);
    int ***p;
    int i;

    p = (int ***)malloc(sizeof(int **) * n);

    testptr((char *)p);
    
    for (i=0; i < n; i++)
	p[i] = i2alloc(m, o);

    return(p);
}

int **i2alloc(int n, int m) {
    void testptr(char *ptr);
    int *ialloc(int n);
    int **p;
    int i;

    p = (int **)malloc(sizeof(int *) * n);

    // p = new double[n][m];

    testptr((char *)p);
    
    for (i=0; i < n; i++)
	p[i] = ialloc(m);

    return(p);
}

double *dalloc(int n) {
    double *p;
    int i;

    p = (double *)malloc(sizeof(double) * n);
    for (i=0; i < n; i++) p[i] = 0.0;

    if (p == NULL) {
	fprintf(stderr, "Error allocating memory.\n");
	exit(1);
    }
    else return(p);
}

double **d2alloc(int n, int m) {
    void testptr(char *ptr);
    double *dalloc(int n);
    double **p;
    int i;

    p = (double **)malloc(sizeof(double *) * n);

    // p = new double[n][m];

    testptr((char *)p);
    
    for (i=0; i < n; i++)
        p[i] = dalloc(m);

    // printf("[[%d]]\n", (int)(sizeof(*p)/sizeof(double)));

    return(p);
}

char *stralloc(int n) {
    char *p;
    int i;

    p = (char *)malloc(sizeof(char) * n);
    for (i=0; i < n; i++) p[i] = ' ';

    if (p == NULL) {
	fprintf(stderr, "Error allocating memory.\n");
	exit(1);
    }
    else return(p);
}

int *ialloc(int n) {
    int *p, i;

    p = (int *)malloc(sizeof(int) * n);

    for (i=0; i < n; i++) p[i] = 0;

    if (p == NULL) {
	fprintf(stderr, "Error allocating memory.\n");
	exit(1);
    }
    else return(p);
}

/* Given an array of x values and another of y values, this finds
   the value of y at a given x. Assumes x is monotonically increasing.  */

double findy(double *x, double *y, int max, double x0) {
    int i, i0 = -1;
    double linint(double y1, double y2, double x1, double x2, double x);
    double ynew=0.0;
/*
    if (x0 < x[0]) return(y[0]);
    if (x0 >= x[max-1]) return(y[max-1]);
*/
    for (i=0; i < max-1; i++) 
	if (x0 >= x[i] && x0 < x[i+1]) {
	    ynew = linint(y[i], y[i+1], x[i], x[i+1], x0);
	    i0 = i;
	}
    
    if (i0 == -1) {
	i--;
	ynew = linint(y[i], y[i+1], x[i], x[i+1], x0);
    }

    return(ynew);
}

/* Given an array of x values and another of y values, this finds
   the value of y at a given x. Assumes x is monotonically increasing.
   With this version of the routine, if it gets a negative answer, 
   and the values are positive, it tries using logarithmic interpolation. */

double findypos(double *x, double *y, int max, double x0) {
    int i, i0 = -1;
    double linint(double y1, double y2, double x1, double x2, double x);
    double logint(double y1, double y2, double x1, double x2, double x);
    double ynew=0.0;
/*
    if (x0 < x[0]) return(y[0]);
    if (x0 >= x[max-1]) return(y[max-1]);
*/
    for (i=0; i < max-1; i++) 
	if (x0 >= x[i] && x0 < x[i+1]) {
	    ynew = linint(y[i], y[i+1], x[i], x[i+1], x0);
	    i0 = i;
	}

    if (i0 == -1) {
	i--;
	ynew = linint(y[i], y[i+1], x[i], x[i+1], x0);
    }
/*
    printf("%12.4le %12.4le %12.4le \n%12.4le %12.4le %12.4le \n",
	   x[i], x0, x[i+1], y[i], ynew, y[i+1]);
*/
    if (ynew < 0.0 && y[i] > 0.0 && y[i+1] > 0.0)
	ynew = logint(y[i], y[i+1], x[i], x[i+1], x0);

    return(ynew);
}

/* Finds the max of a specified number of doubles. */

double dmax(int n, double a, ...) {
    double max = a, b;
    int i;
    va_list ap;

    va_start(ap, a);
    
    for (i = 0; i < n-1; i++) {
	b = va_arg(ap, double);
	if (b > max) max = b;
    }

    return(max);
}

/* Finds the min of a specified number of doubles. */

double dmin(int n, double a, ...) {
    double min = a, b;
    int i;
    va_list ap;

    va_start(ap, a);
    
    for (i = 0; i < n-1; i++) {
	b = va_arg(ap, double);
	if (b < min) min = b;
    }

    return(min);
}

/* prints a Usage error message. */

void usage(int argc, char **argv, int nargs, char *usage_string) {
    if (argc != nargs) {
	fprintf(stderr, "Error: Incorrect number of arguments.\n");
	fprintf(stderr, "Usage: %s %s\n", argv[0], usage_string);
	exit(1);
    }
}

/* Counts the number of records in f and closes the file after. */

int countlines(char *name) {
  char c; 
  int l = 0;
  FILE *f;
  FILE   *safeopen(char *s, char *f);

  f = safeopen(name, "r");
  while ((c = fgetc(f))!=EOF) if (c == '\n') l++;
  fclose(f);
  return(l);
}

/* check to make sure an array was dynamically allocated correctly */

void testptr(char *ptr)
{
  if (ptr == NULL) {
    fprintf(stderr, "Memory allocation error.\n"); exit(0);
  }
}

/* return the max of a and b */

double max(double a, double b) {
  if (a > b) return(a); else return(b);
}

/* return the min of |a| and |b| */

double min(double a, double b) {
  /* if (fabs(a) < fabs(b)) return(a); else return(b); */
  if (fabs(a) < fabs(b)) return(a); else return(b);
}

/* This fn does linear interpolation */

double linint(double y1, double y2, double x1, double x2, double x) 
{
  double m, b, y;

  if (fabs(x2-x1) < MYTINY) m = MYHUGE; 
  else m = (y2-y1)/(x2-x1); 
  b = y1-m*x1;
  y = m*x + b;
  return(y);
}

/* Logarithmic interpolation */

double logint(double y1, double y2, double x1, double x2, double x) 
     /* This fn does logarithmic interpolation
	to find y(x) from x1,x2,y1,y2 */
{
  double A, tau, y;
  /* done by fitting y(x) to y(x) = A*exp(tau*x) */

  /* note that if y1/y2 <= 0 you're headin' for a mess o' trouble */

  if (fabs(x2-x1) < MYTINY) tau = MYHUGE; 
  else tau = log(y1/y2) / (x1-x2);
  A = y1*exp(-tau*x1);
  y = A*exp(tau*x);
  return(y);
}

/* print an error message and quit */

void printerr(char *s) {
  fprintf(stderr, "%s\n", s);
  exit(0);
}

/* read from file f until you see the string s, or until the end of file */

void queueup(FILE *f, FILE *g, char *s) {
  char addchar(char *s, int l, char c);
  char c;
  int i, L;
  char *S;
  unsigned char eof;

  eof=(unsigned char)EOF;

  L = strlen(s);

  /* allocate memory: */
  if ((S=(char *)malloc((unsigned)(L+1)))==NULL) {
    fprintf(stderr, "Error in queueup() allocating memory\n");
    exit(1);
  }

  S[L]='\0';
  /* read in the first L characters: */
  for (i=0; i < L; i++) {
    S[i] = fgetc(f);
    if (g != NULL) fputc(S[i], g);
  }

  if (strcmp(s,S)) /* if the first L chars match S then quit */

/* WARNING: On some systems, EOF should be used; on others, eof, which
   is the unsigned equivalent of EOF, defined above. */

      /* while (strcmp(s,S) && (c=fgetc(f))!=eof) { */
      while (strcmp(s,S) && (c=fgetc(f))!=EOF) { 
	  addchar(S, L, c);
	  if (g != NULL) fputc(c, g);
      }
}

// This version is not sensitive to case:
// Returns true if not at the end of the file; false otherwise.

int queueuplc(FILE *f, FILE *g, char *s) {
  char addchar(char *s, int l, char c);
  char c;
  int i, L;
  char *S;
  unsigned char eof;
  char tolc(char c);

  eof=(unsigned char)EOF;

  L = strlen(s);

  /* allocate memory: */
  if ((S=(char *)malloc((unsigned)(L+1)))==NULL) {
    fprintf(stderr, "Error in queueup() allocating memory\n");
    exit(1);
  }

  S[L]='\0';
  /* read in the first L characters: */
  for (i=0; i < L; i++) {
      S[i] = tolc(fgetc(f));
      if (g != NULL) fputc(S[i], g);
  }
  
  if (strcmp(s,S)) /* if the first L chars match S then quit */

/* WARNING: On some systems, EOF should be used; on others, eof, which
   is the unsigned equivalent of EOF, defined above. */

      /* while (strcmp(s,S) && (c=fgetc(f))!=eof) { */
      while (strcmp(s,S) && (c=tolc(fgetc(f)))!=EOF) { 
	  addchar(S, L, c);
	  if (g != NULL) fputc(c, g);
      }

  c = fgetc(f);
  if (c != EOF) ungetc(c, f);

  return(!(c==EOF));
}

// Returns true if not at the end of the file; false otherwise.

int queueupFileStatus(FILE *f, FILE *g, char *s) {
  char addchar(char *s, int l, char c);
  char c;
  int i, L;
  char *S;
  unsigned char eof;

  eof=(unsigned char)EOF;

  L = strlen(s);

  /* allocate memory: */
  if ((S=(char *)malloc((unsigned)(L+1)))==NULL) {
    fprintf(stderr, "Error in queueup() allocating memory\n");
    exit(1);
  }

  S[L]='\0';
  /* read in the first L characters: */
  for (i=0; i < L; i++) {
      S[i] = fgetc(f);
      if (g != NULL) fputc(S[i], g);
  }
  
  if (strcmp(s,S)) /* if the first L chars match S then quit */

/* WARNING: On some systems, EOF should be used; on others, eof, which
   is the unsigned equivalent of EOF, defined above. */

      /* while (strcmp(s,S) && (c=fgetc(f))!=eof) { */
      while (strcmp(s,S) && (c=fgetc(f))!=EOF) { 
	  addchar(S, L, c);
	  if (g != NULL) fputc(c, g);
      }

  c = fgetc(f);
  if (c != EOF) ungetc(c, f);

  return(!(c==EOF));
}

// Convert a character to lower case.

char tolc(char c) {
    if (c >= 'A' && c <= 'Z') 
	return(c - ('A'-'a'));
    else
	return(c);
}

/* This version uses an unsigned char for EOF; necessary on some
   systems. */

void queueup_unsgn(FILE *f, FILE *g, char *s) {
  char addchar(char *s, int l, char c);
  char c;
  int i, L;
  char *S;
  unsigned char eof;

  eof=(unsigned char)EOF;

  L = strlen(s);

  /* allocate memory: */
  if ((S=(char *)malloc((unsigned)(L+1)))==NULL) {
    fprintf(stderr, "Error in queueup() allocating memory\n");
    exit(1);
  }

  S[L]='\0';
  /* read in the first L characters: */
  for (i=0; i < L; i++) {
    S[i] = fgetc(f);
    if (g != NULL) fputc(S[i], g);
  }

  if (strcmp(s,S)) /* if the first L chars match S then quit */

/* WARNING: On some systems, EOF should be used; on others, eof, which
   is the unsigned equivalent of EOF, defined above. */

      while (strcmp(s,S) && (c=fgetc(f))!=eof) { 
	  addchar(S, L, c);
	  if (g != NULL) fputc(c, g);
      }
}

/* shift the string so that you remove the first character and add c
   on to the end */

char addchar(char *s, int l, char c) {
  int i;

  for (i=0; i < l-1; i++) s[i]=s[i+1];
  s[l-1]=c;
  return(c);
}

/* This uses the system command rand() to generate a (repeatably) random
   double between 0 and 1. */

double myrand(void) {
    static int firsttime=1;

    if (firsttime) {
	srand(987654321);
	firsttime=0;
    }

    return((double)rand()/((double)RAND_MAX));
}
/* This uses the system command rand() to generate a (repeatably) random
   double between 0 and 1. */

double myRandSeed(int seed) {
    static int firsttime=1;

    if (firsttime) {
	srand(seed);
	firsttime=0;
    }

    return((double)rand()/((double)RAND_MAX));
}

/* Read from file f until you see the string s, or until the end of file. */
/* If g != NULL, the characters are echoed to string g. */
void queueupsN(FILE *f, char *g, char *s, int M) {
  char addchar(char *s, int l, char c);
  char c;
  int i, L, gi=0, len=0;
  char *S;
  unsigned char eof;

  eof=(unsigned char)EOF;

  L = strlen(s);

  /* allocate memory: */
  if ((S=(char *)malloc((unsigned)(L+1)))==NULL) {
    fprintf(stderr, "Error in queueupsN() allocating memory\n");
    exit(1);
  }

  S[L]='\0';
  /* read in the first L characters: */
  for (i=0; i < L; i++) {
    S[i] = fgetc(f);
    if (g != NULL) 
      g[gi++] = S[i];
    /* fputc(S[i], g); */
  }

  len = L;

  if (strcmp(s,S)) /* if the first L chars match S then quit */
    while (strcmp(s,S) && (c=fgetc(f))!=eof && len < M) { 
      addchar(S, L, c);
      if (g != NULL) { 
	g[gi++] = c;
	len++;
      }
      /* fputc(c, g); */
    }
  if (g != NULL) g[gi++] = '\0';
  free(S);
}

/* read from file f until you see the string s, or until the end of file */
/* if g != NULL, the characters are echoed to g. */
void queueupN(FILE *f, FILE *g, char *s, int M) {
  char addchar(char *s, int l, char c);
  char c;
  int i, L, len=0;
  char *S;
  unsigned char eof; 

  eof=(unsigned char)EOF;

  L = strlen(s);

  /* allocate memory: */
  if ((S=(char *)malloc((unsigned)(L+1)))==NULL) {
    fprintf(stderr, "Error in queueup() allocating memory\n");
    exit(1);
  }

  S[L]='\0';
  /* read in the first L characters: */
  for (i=0; i < L; i++) {
    S[i] = fgetc(f);
    if (g != NULL) fputc(S[i], g);
  }

  len=L;

  if (strcmp(s,S)) /* if the first L chars match S then quit */
    while (strcmp(s,S) && (c=fgetc(f))!=eof && len < M && !feof(f)) { 
      addchar(S, L, c);
      /* printf("%d ", j++); */
      if (g != NULL) {
	fputc(c, g);
	len++;
      }
    }

  if (len == M || c == eof || feof(f)) {
    fprintf(stderr, "Error -- string %s not found in input file.\n",
	    S);
  }

  free(S);
}

/* Scan string f until you hit the string s or the end of the string,
   copying the text to g in the process. This function returns a pointer
   to the next character, or NULL if it's reached the end of the string. */
char *squeueupsN(char *f, char *g, char *s, int M) {
  char c;
  int i, slen, flen;

  slen = strlen(s);
  flen = strlen(f);
  
  for (i=0; i < flen && i < M && strncmp(s, &(f[i]), slen); i++)
      if (g != NULL) g[i] = f[i];

  if (!strncmp(s, &(f[i]), slen)) {
      if (g != NULL) g[i+slen] = '\0';
      return(&(f[i+slen]));
  } else if (i == strlen(f)) {
      if (g != NULL) g[i] = '\0';
      return(NULL);
  } else {
      if (g != NULL) g[i+1] = '\0';
      return(&(f[i+1]));
  }
}

int iswhitespace(char c) {
    return(c == ' ' || c == '\t' || c == '\n');
}

// Magnitude of the cross product
double crossProductMag(double x1, double y1, double z1, double x2, double y2, double z2) {
    double a, b, c;
    a = y1 * z2 - z1 * y2;
    b = z1 * x2 - x1 * z2;
    c = x1 * y2 - y1 * x2;
    
    return(sqrt(a*a + b*b + c*c));
}

// Area of a planar quadrilateral; points must be ordred either clockwise or ccw.
double quadArea(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
    return(triangleArea(x1, y1, x2, y2, x3, y3)
	   +triangleArea(x3, y3, x4, y4, x1, y1));
}

double triangleArea(double x1, double y1, double x2, double y2, double x3, double y3) {
    return(0.5 * crossProductMag(x2-x1, y2-y1, 0.0, x3-x1, y3-y1, 0.0));
}

// Perform interpolation to point (x, y) based on values fi at (xi, yi).
// To do this maps the quad onto a unit square with coordinates (l, m).
// Fails for degenerate quads. 
// Assumes points listed in ccw or cw order
// See https://www.particleincell.com/2012/quad-interpolation/ for derivation, which makes reference
// to Birdsall, C. K., and Langdon, A.B., “Plasma Physics Via Computer Simulations“, Institute of Physics Publishing, 2000
//    Hughes, T. J. R., “The Finite Element Method“, Dover Publications, 2000
double fourPtInterp(double x1, double y1, double f1, 
		    double x2, double y2, double f2, 
		    double x3, double y3, double f3, 
		    double x4, double y4, double f4, 
		    double x,  double y) {
    double w1, w2, w3, w4; // quadrant areas
    double a[5], b[5];
    double aa, bb, cc, l, m, det;
    
    a[1] =  x1;
    a[2] = -x1 + x2;
    a[3] = -x1 + x3;
    a[4] =  x1 - x2 + x3 - x4;
    
    b[1] =  y1;
    b[2] = -y1 + y2;
    b[3] = -y1 + y3;
    b[4] =  y1 - y2 + y3 - y4;
    
    aa = a[4] * b[3] - a[3] * b[4];
    bb = a[4] * b[1] - a[1] * b[4] + a[2] * b[3] - a[3] * b[2] + x * b[4] - y * a[4];
    cc = a[2] * b[1] - a[1] * b[2] + x * b[2] - y * a[2];
 
    // compute m = (-b+sqrt(b^2-4ac))/(2a)
    det = sqrt(bb * bb - 4 * aa * cc);
    m   = (-bb + det) / (2 * aa);
 
    // compute l
    l   = (x - a[1] - a[3] * m) / (a[2] + a[4] * m);

    w1  = l * m;
    w2  = (1.0 - l) * m;
    w3  = (1.0 - l) * (1.0 - m);
    w4  = l * (1.0 - m);

    // return(w1 * f3 + w2 * f4 + w3 * f1 + w4 * f2);

    return(0.25 * (f1 + f2 + f3 + f4));
}

// Determine if a 2D point is within a triangle               
//     This function returns true if (rx, ry) is in the triangle formed by the points 
//     (x1, y1), ..., (x3, y3).
bool inTriangle(double xx1, double yy1, double xx2, double yy2, double xx3, double yy3, 
		double rxx, double ryy, bool verbose) {
    double sp, sq, pxx, pyy, qxx, qyy, det, sxx, syy;
    bool   in_triangle;
         
    pxx = xx2 - xx1;
    pyy = yy2 - yy1;
    qxx = xx3 - xx1;
    qyy = yy3 - yy1;
    sxx = rxx - xx1;
    syy = ryy - yy1;

    det = pxx * qyy - qxx * pyy;

    if (fabs(det)/ sqrt(pxx * pxx + pyy * pyy + qxx * qxx + qyy * qyy) < 1.0e-12) {
        // The triangle is flat; so is the point colinear with the triangle?
	if (verbose) {
	    fprintf(stderr, "Flat triangle detected in inTriangle()\n");
	    fprintf(stderr, "%12.el4 %12.el4 %12.el4 %12.el4 %12.el4 %12.el4 %12.el4 \n", 
		    pxx, pyy, qxx, qyy, sxx, syy, det);
	}
	if (fabs(syy * pxx - sxx * pyy) < 1.0e-10 * sxx * pxx)
	    in_triangle = true;
	else
	    in_triangle = false;
    } else {

        // Find the coordinates (sp, sq) of the vector (rxx, ryy)-(xx1, yy1) in the
        // basis of (xx2, yy2)-(xx1, yy1) and (xx3, yy3)-(xx1, yy1):

	sp = ( sxx * qyy - syy * qxx) / det;
	sq = (-sxx * pyy + pxx * syy) / det;
	if (sp < 0.0 || sq < 0.0 || (sp + sq - 1.0) > 1.0e-10) 
	    in_triangle = false;
	else
	    in_triangle = true;

	if (verbose) {
	    fprintf(stderr,  "--> %12.4le %12.4le %12.4le %12.4le %12.4le %12.4le %12.4le %12.4le %12.4le\n", 
		    sp, sq, pxx, pyy, qxx, qyy, sxx, syy, det);
	    fprintf(stderr, "zone f=point\n");
            fprintf(stderr, "%12.4le %12.4le \n", rxx, ryy);
	    fprintf(stderr, "zone f=point\n");
	    fprintf(stderr, "%12.4le %12.4le \n", xx1, yy1);
	    fprintf(stderr, "%12.4le %12.4le \n", xx2, yy2);
	    fprintf(stderr, "%12.4le %12.4le \n", xx3, yy3);
	    fprintf(stderr, "%12.4le %12.4le \n", xx1, yy1);
	    fprintf(stderr, "zone f=point\n");
	    fprintf(stderr, "%12.4le %12.4le \n", sp, sq);
	}
    }
    return(in_triangle);
}
				      
bool inQuadrilateral(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, 
		     double x, double y, bool verbose) {
	 return(inTriangle(x1, y1, x2, y2, x3, y3, x, y, verbose) ||
		inTriangle(x1, y1, x3, y3, x4, y4, x, y, verbose));
}

// From Numerical Recipes 3rd ed., this generates the *renomralized* associated Legendre functions--see
// Eq. 6.7.8, p. 294. This is related to the spherical harmonics as follows: 
// 
//   Ylm(theta, phi) = plgndr(l, m, x) * exp(im * phi)
//
// where x = cos(theta). This differs from the associated Legengre functions in that it has a factor of
// sqrt((2l+1)(l-m)!/(4pi(l+m)!)) for normalization.
double plegendre(const int l, const int m, const double x) {
    int i,ll;
    double fact,oldfact,pll,pmm,pmmp1,omx2;
    if (m < 0 || m > l || fabs(x) > 1.0) {
	fprintf(stderr, "Bad arguments in routine plegendre()\n");
	fprintf(stderr, "l = %d, m = %d, x = %g\n", l, m, x);
	exit(1);
    }
    pmm=1.0;
    if (m > 0) {
	omx2=(1.0-x)*(1.0+x);
	fact=1.0;
	for (i=1;i<=m;i++) {
	    pmm *= omx2*fact/(fact+1.0);
	    fact += 2.0;
	}
    }
    pmm=sqrt((2*m+1)*pmm/(4.0*M_PI));
    if (m & 1)
	pmm=-pmm;
    if (l == m)
	return pmm;
    else {
	pmmp1=x*sqrt(2.0*m+3.0)*pmm;
	if (l == (m+1))
	    return pmmp1;
	else {
	    oldfact=sqrt(2.0*m+3.0);
	    for (ll=m+2;ll<=l;ll++) {
		fact=sqrt((4.0*ll*ll-1.0)/(ll*ll-m*m));
		pll=(x*pmmp1-pmm/oldfact)*fact;
		oldfact=fact;
		pmm=pmmp1;
		pmmp1=pll;
	    }
	    return pll;
	}
    }
}

// This routine from Numerical Recipes in C, generates the associated Legrendre functions, which are related
// to the spherical harmonics via 
//
//   Ylm(theta, phi) = sqrt((2l+1)(l-m)!/(4pi(l+m)!)) plgndr(l,m,cos(theta) exp(i m phi)
//
double plgndr(int l, int m, double x)
{
    double fact,pll,pmm,pmmp1,somx2;
    int i,ll;
    
    if (m < 0 || m > l || fabs(x) > 1.0) {
	fprintf(stderr, "Bad arguments in routine PLGNDR");
	exit(1);
    }
    pmm=1.0;
    if (m > 0) {
	somx2=sqrt((1.0-x)*(1.0+x));
	fact=1.0;
	for (i=1;i<=m;i++) {
	    pmm *= -fact*somx2;
	    fact += 2.0;
	}
    }
    if (l == m)
	return pmm;
    else {
	pmmp1=x*(2*m+1)*pmm;
	if (l == (m+1))
	    return pmmp1;
	else {
	    for (ll=(m+2);ll<=l;ll++) {
		pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
		pmm=pmmp1;
		pmmp1=pll;
	    }
	    return pll;
	}
    }
}

// The following is the adaptive integrator taken from Numerical Recipes 3rd ed. p. 292.
// NR3 uses templates, which are not sufficiently universal for me to use, so I have 
// modified the code so it only works with type double, and uses a structure, adaptParms,
// to hold the member variables which were previously part of the Adapt class.
double integrate(double (*func)(double), const double a, const double b, double tol) {
	double m,h,fa,fb,i1,i2,is,erri1,erri2,r,y[13];
	const double alpha=sqrt(2.0/3.0);
	const double beta=1.0/sqrt(5.0);
	const double x1=0.942882415695480;
	const double x2=0.641853342345781;
	const double x3=0.236383199662150;
	const double x[12]={0,-x1,-alpha,-x2,
				  -beta,-x3,0.0,x3,
				  beta,x2,alpha,x1};

	struct adaptParms p;

	// Populate the structure containing the integrator parameters:
	p.terminate = true;
	p.out_of_tolerance = false;
	p.TOL = tol;

	m = 0.5*(a+b);
	h = 0.5*(b-a);
	fa = y[0]=(*func)(a);
	fb = y[12]=(*func)(b);
	for (int i=1; i < 12; i++)
	    y[i] = (*func)(m+x[i] * h);
	i2 = (h/6.0)*(y[0]+y[12]+5.0*(y[4]+y[8]));
	i1 = (h/1470.0)*(77.0*(y[0]+y[12])+432.0*(y[2]+y[10]) + 
			 625.0*(y[4]+y[8])+672.0*y[6]);
	is = h*(0.0158271919734802*(y[0]+y[12])+0.0942738402188500*
		(y[1]+y[11])+0.155071987336585*(y[2]+y[10])+
		0.188821573960182*(y[3]+y[9])+0.199773405226859*
		(y[4]+y[8])+0.224926465333340*(y[5]+y[7])+
		0.242611071901408*y[6]);
	erri1 = fabs(i1-is);
	erri2 = fabs(i2-is);
	r = (erri2 != 0.0) ? erri1/erri2 : 1.0;
	p.toler = (r > 0.0 && r < 1.0) ? p.TOL/r : p.TOL;
	if (is == 0.0)
	    is = b - a;
	is = fabs(is);

	return adaptlob(func,a,b,fa,fb,is, p);
}

double adaptlob(double (*func)(double), const double a, const double b, const double fa,
		const double fb, const double is, struct adaptParms &p) {
	double m,h,mll,ml,mr,mrr,fmll,fml,fm,fmrr,fmr,i1,i2;
	const double alpha=sqrt(2.0/3.0);
	const double beta=1.0/sqrt(5.0);
	const double x1=0.942882415695480;
	const double x2=0.641853342345781;
	const double x3=0.236383199662150;
	const double x[12]={0,-x1,-alpha,-x2,-beta,-x3,0.0,x3,beta,x2,alpha,x1};

	m = 0.5*(a+b);
	h = 0.5*(b-a);
	mll = m-alpha*h;
	ml = m-beta*h;
	mr = m+beta*h;
	mrr = m+alpha*h;
	fmll = (*func)(mll);
	fml = (*func)(ml);
	fm = (*func)(m);
	fmr = (*func)(mr);
	fmrr = (*func)(mrr);
	i2 = h/6.0*(fa+fb+5.0*(fml+fmr));
	i1 = h/1470.0*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);
	if (fabs(i1-i2) <= p.toler*is || mll <= a || b <= mrr) {
	    if ((mll <= a || b <= mrr) && p.terminate) {
		p.out_of_tolerance=true;
		p.terminate=false;
	    }
	    return i1;
	}
	else
	    return adaptlob(func,a,mll,fa,fmll,is,p)+
		adaptlob(func,mll,ml,fmll,fml,is,p)+
		adaptlob(func,ml,m,fml,fm,is,p)+
		adaptlob(func,m,mr,fm,fmr,is,p)+
		adaptlob(func,mr,mrr,fmr,fmrr,is,p)+
		adaptlob(func,mrr,b,fmrr,fb,is,p);
}

double integrate2args(double (*func)(double, double), const double a, const double b, double tol, double phi) {
	double m,h,fa,fb,i1,i2,is,erri1,erri2,r,y[13];
	const double alpha=sqrt(2.0/3.0);
	const double beta=1.0/sqrt(5.0);
	const double x1=0.942882415695480;
	const double x2=0.641853342345781;
	const double x3=0.236383199662150;
	const double x[12]={0,-x1,-alpha,-x2,
				  -beta,-x3,0.0,x3,
				  beta,x2,alpha,x1};

	struct adaptParms2args p;

	// Populate the structure containing the integrator parameters:
	p.terminate        = true;
	p.out_of_tolerance = false;
	p.TOL              = tol;
	p.phi              = phi;

	m = 0.5*(a+b);
	h = 0.5*(b-a);
	fa = y[0]=(*func)(a, phi);
	fb = y[12]=(*func)(b, phi);
	for (int i=1; i < 12; i++)
	    y[i] = (*func)(m+x[i] * h, phi);
	i2 = (h/6.0)*(y[0]+y[12]+5.0*(y[4]+y[8]));
	i1 = (h/1470.0)*(77.0*(y[0]+y[12])+432.0*(y[2]+y[10]) + 
			 625.0*(y[4]+y[8])+672.0*y[6]);
	is = h*(0.0158271919734802*(y[0]+y[12])+0.0942738402188500*
		(y[1]+y[11])+0.155071987336585*(y[2]+y[10])+
		0.188821573960182*(y[3]+y[9])+0.199773405226859*
		(y[4]+y[8])+0.224926465333340*(y[5]+y[7])+
		0.242611071901408*y[6]);
	erri1 = fabs(i1-is);
	erri2 = fabs(i2-is);
	r = (erri2 != 0.0) ? erri1/erri2 : 1.0;
	p.toler = (r > 0.0 && r < 1.0) ? p.TOL/r : p.TOL;
	if (is == 0.0)
	    is = b - a;
	is = fabs(is);

	return adaptlob2args(func,a,b,fa,fb,is, p);
}

double adaptlob2args(double (*func)(double, double), const double a, const double b, const double fa,
		const double fb, const double is, struct adaptParms2args &p) {
	double m,h,mll,ml,mr,mrr,fmll,fml,fm,fmrr,fmr,i1,i2;
	const double alpha=sqrt(2.0/3.0);
	const double beta=1.0/sqrt(5.0);
	const double x1=0.942882415695480;
	const double x2=0.641853342345781;
	const double x3=0.236383199662150;
	const double x[12]={0,-x1,-alpha,-x2,-beta,-x3,0.0,x3,beta,x2,alpha,x1};

	m = 0.5*(a+b);
	h = 0.5*(b-a);
	mll = m-alpha*h;
	ml = m-beta*h;
	mr = m+beta*h;
	mrr = m+alpha*h;
	fmll = (*func)(mll, p.phi);
	fml = (*func)(ml, p.phi);
	fm = (*func)(m, p.phi);
	fmr = (*func)(mr, p.phi);
	fmrr = (*func)(mrr, p.phi);
	i2 = h/6.0*(fa+fb+5.0*(fml+fmr));
	i1 = h/1470.0*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);
	if (fabs(i1-i2) <= p.toler*is || mll <= a || b <= mrr) {
	    if ((mll <= a || b <= mrr) && p.terminate) {
		p.out_of_tolerance=true;
		p.terminate=false;
	    }
	    return i1;
	}
	else
	    return(
		adaptlob2args(func,a,mll,fa,fmll,is,p) +
		adaptlob2args(func,mll,ml,fmll,fml,is,p) +
		adaptlob2args(func,ml,m,fml,fm,is,p) +
		adaptlob2args(func,m,mr,fm,fmr,is,p) +
		adaptlob2args(func,mr,mrr,fmr,fmrr,is,p) +
		adaptlob2args(func,mrr,b,fmrr,fb,is,p)
		);
}

