// Library functions:

// System and file operations:
char    *getUserName(char **env);

// Streams and files:
FILE    *nc(FILE *f);
FILE    *safeopen(char *s, char *f);
FILE    *safeopentb(char *s, char *f, char *tb);
int     fileExists(char *f);
int     countlines(char *name);
char    *newExt(char *f, char *n, int N);
char    *getline(FILE *fp, char *l);
char    *getlineEOF(FILE *fp, char *l, bool &eof);
bool    passLine(FILE *fin, FILE *fout);
void    waitForFileCreation(char *s);

// Parsing streams and string:
int     queueupFileStatus(FILE *f, FILE *g, char *s);
int     queueuplc(FILE *f, FILE *g, char *s);
void    queueup(FILE *f, FILE *g, char *s);
void    queueup_unsgn(FILE *f, FILE *g, char *s);
void    queueupN(FILE *f, FILE *g, char *s, int M);
void    queueupsN(FILE *f, char *g, char *s, int M);
char    *squeueupsN(char *f, char *g, char *s, int M);
char    *getword(char *w, char *s);
double  fetch(char *filename, char *string, bool &found);

// String and characters:
char    *reverse(char *s);
int     iswhitespace(char c);

// Mathemtical functions:
double  maxAbsList(double *v, int N);
double  windowFnSymm(double x, double xmax);
double  windowFn(double x, double xmin, double xmax);
double  windowShift(double x, double x0, double x1);
double  logint2D(int nx, int ny, double *x, double *y, double **f, 
		double x0, double y0);
double  linint2D(int nx, int ny, double *x, double *y, double **f, 
		double x0, double y0);
double  findypos(double *x, double *y, int max, double x0);
double  findy(double *x, double *y, int max, double x0);
double  sign(double x);
double  dmax(int, double , ...);
double  dmin(int, double , ...);
double  max(double a, double b);
double  min(double a, double b);
double  logint(double y1, double y2, double x1, double x2, double x);
double  linint(double y1, double y2, double x1, double x2, double x);
double  myrand(void);
double  myRandSeed(int seed);
double  distsq(double x, double y, double z);
double  dist(double x, double y, double z);
double  plgndr(int l, int m, double x); // from NR
double  plegendre(const int l, const int m, const double x);
double  crossProductMag(double x1, double y1, double z1, double x2, double y2, double z2);
double  triangleArea(double x1, double y1, double x2, double y2, double x3, double y3);
double  quadArea(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
double  fourPtInterp(double x1, double y1, double f1, 
		    double x2, double y2, double f2, 
		    double x3, double y3, double f3, 
		    double x4, double y4, double f4, 
		    double x, double y);

// The following three are for an adaptation of Numerican Recipes' adaptive integrator:
struct  adaptParms {
    double TOL,toler;
    bool terminate,out_of_tolerance;
};
double  integrate(double (*func)(double), const double a, const double b, double tol);
double  adaptlob(double (*func)(double), const double a, const double b, const double fa,
		const double fb, const double is, struct adaptParms &p);
// The following three are for an adaptation of Numerican Recipes' adaptive integrator:
struct  adaptParms2args {
    double TOL,toler;
    bool terminate,out_of_tolerance;
    double phi;
};
double  integrate2args(double (*func)(double, double), const double a, const double b, double tol, double phi);
double  adaptlob2args(double (*func)(double, double), const double a, const double b, const double fa,
		const double fb, const double is, struct adaptParms2args &p);

bool    inQuadrilateral(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, 
		     double x, double y, bool verbose);
bool    inTriangle(double xx1, double yy1, double xx2, double yy2, double xx3, double yy3, 
		double rxx, double ryy, bool verbose);
void    aitoff(double theta, double phi, double *x, double *y);
void    rotPhi(double &x, double &y, double &z, double phi);
void    rotTheta(double &x, double &y, double &z, double theta);
void    unitVector(double &x, double &y, double &z, double theta, double phi);
void    rotyaxis(double &x, double &y, double &z, double theta);


// Plotting:
void tecplot3d(char *filename, char *title, 
	       char *xname, char *yname, char *zname, char *fname, 
	       double *x, double *y, double *z, double ***q, 
	       int xN, int yN, int zN);
void tecplot2d(char *filename, char *title, 
	       char *xname, char *yname, char *fname, 
	       double *x, double *y, double **q, 
	       int xN, int yN);

// Memory allocation:
int     *ialloc(int n);
int     **i2alloc(int n, int m);
int     ***i3alloc(int n, int m, int o);
double  *dalloc(int n);
double  **d2alloc(int n, int m);
double  ***d3alloc(int n, int m, int o);
double  ****d4alloc(int n, int m, int o, int l);
double  *****d5alloc(int n, int m, int o, int l, int q);
void    d2free(double **p, int n);
void    d3free(double ***p, int n, int m);
void    d4free(double ****p, int n, int m, int o);
void    d5free(double *****p, int n, int m, int o, int q);
void    testptr(char *ptr);
char    *stralloc(int n);

// Utilities:
void    usage(int argc, char **argv, int nargs, char *usage_string);
void    printerr(char *s);

// Related to port-config files and spot shapes:
static const int MAXRINGS       = 200;    // max # of rings in the emod and pointing files
static const int MAXPORTS       = 1000;  // max # of points in pointing config file
static const int MAXPCFFILESIZE = 100000; // max input file size (char)

void   read_pcf(char *fname, struct pointingFile &p, bool SpotMask);
void   write_pcf(char *fname, struct pointingFile &p);
void   getEmods(char *filename, pointingFile &fp);

struct pointingFile {
    char header[MAXPCFFILESIZE];

    // int r1[MAXPORTS], r2[MAXPORTS], r3[MAXPORTS];
    int rn[MAXRINGS][MAXPORTS]; // ports in ea ring. NOT zero-referenced!
    int nports[MAXRINGS]; // number of ports in each ring
    int block_sizes[MAXRINGS + 11];
    int nrings;
    int ringsByPort[MAXPORTS]; // the 0-referenced ring for ea. port (and the port #s are 0-reffed)

    int totalNumPorts;

    // This is a block describing the spot shapes and secondary ellipse properties:
    double p[MAXRINGS][7]; 
    double theta0[MAXPORTS], phi0[MAXPORTS]; // original port angles (deg)
    double theta1[MAXPORTS], phi1[MAXPORTS]; // reported port angles (deg)

    // The offsets (in "units" of the target radius) are a calculated quantity:
    double xcent[MAXPORTS], ycent[MAXPORTS];   // The new beam centers corresponding to the repointings
    double spotIntegral[MAXPORTS]; // This is the integral of the generalized spot over the 
                                   // far-field coordinates (computed when file is read in)
    double dLambdaIR[MAXPORTS];    // Detuning amount (A)

    bool   spotMask; // True if spot masking is incorporated
    double RSM;      // spot-masking radius relative to target radius, read from sim-input file
    double xoffSM;   // offset of spot-mask in "x" direction
    double SGSM;     // super-Gaussian order of the spot masking function

    double emods[MAXRINGS];

    bool radiusInterpPresent; // True if the single digit prior to the spot properties block is present 
                              // saying how the radii are to be interpreted 
    int  radiusInterp;
};

enum PCF_SPOT_PARM { PCF_SGexp, PCF_R0, PCF_eta, PCF_f2nd, PCF_SG2nd, PCF_xoff, PCF_eta2nd };

/* NIF port angles, courtesy Reuben, who got them from Stan */
const double theta_nif[] = {
    23.50, 23.50, 23.50, 23.50, 
    44.50, 44.50, 44.50, 44.50, 
    44.50, 44.50, 44.50, 44.50, 
    77.50, 77.50, 77.50, 77.50, 
    77.50, 77.50, 77.50, 77.50, 
    77.50, 77.50, 77.50, 77.50, 
    102.50, 102.50, 102.50, 102.50, 
    102.50, 102.50, 102.50, 102.50, 
    102.50, 102.50, 102.50, 102.50, 
    135.50, 135.50, 135.50, 135.50, 
    135.50, 135.50, 135.50, 135.50, 
    156.50, 156.50, 156.50, 156.50};

const double phi_nif[] = {
    78.750, 168.750,  258.750, 348.750, 
    16.293,  62.455,  106.295, 152.455, 
    196.293, 242.455, 286.295, 332.455, 
    24.375,  54.375,  84.375,  114.375, 
    144.375, 174.375, 204.375, 234.375, 
    264.375, 294.375, 324.375, 354.375, 
    5.625,   35.625,  65.625,  95.625, 
    125.625, 155.625, 185.625, 215.625, 
    245.625, 275.625, 305.625, 335.625, 
    27.545,  73.705,  117.545, 163.707, 
    207.545, 253.705, 297.545, 343.707, 
    11.250,  101.250, 191.250, 281.250};

double generalSpotShape(double x, double y, 
			double SGexp, double SGexp2nd, double SGSM, 
			double eta, double eta2nd, 
			double R0, double RSM, 
			double xcent, double ycent, 
			double xoff, double xoffSM, 
			double f2nd);

double spotIntegral(double SGexp, double SGexp2nd, double SGSM, 
		    double eta, double eta2nd, 
		    double R0, double RSM, 
		    double xcent, double ycent, 
		    double xoff, double xoffSM, 
		    double f2nd);

// Related to spherical projection plots:
struct sphericalProjection {
    int N, M;
    double **theta, **phi, **f; // theta and phi are in rad
};

void     calcSphericalProjection(struct pointingFile &pf, struct sphericalProjection &S, bool kruer, bool useEmods);
void     calcRepointingOffsets(struct pointingFile &p);
void     calcPCFSpotIntegrals(struct pointingFile &p);
void     plotSphericalProjection(struct sphericalProjection &S, char *s, bool thetaOnly, bool append);
void     sphericalProjectionFree(struct sphericalProjection &S);
void     allocSphericalProjection(int N, int M, struct sphericalProjection &S);
void     normalizeSphericalProjection(struct sphericalProjection &S);
double   beam_theta(double theta, double phi, double thetap, double phip);
double   averageSphericalProjection(struct sphericalProjection &S);
double   sigmarmsSphericalProjection(struct sphericalProjection &S, bool thetaOnly);
double   illuminUnif(struct pointingFile &pf, struct sphericalProjection &S, bool thetaOnly);
double   detuningMetric(pointingFile &pf, struct sphericalProjection &S, double dLambdaMaxIR, bool thetaOnly);

struct planeProjection {
    int N, M;
    double **x, **y, **z, **f;
};

void allocPlaneProjection(int N, int M, struct planeProjection &P);

// Related to tecplot files:

/* This data structure contains tecplot data read from an ASCII file as written by Draco, e.g.
   It assumes the zones are all of one type and shape. */
struct tecplotData {
    bool allocated;  // true if the data structures have been allocated already

    char *filename;  // tecplot file name
    int ndim;        // dimensionality (2 or 3)
    int I, J, K;     // array dimensions
    int nvars;       // number of variables
    int nzones;      // number of zones
    char *title;     // overall title
    char **varNames; // variable names
    double    ***d1; // the data cube
    double   ****d2; // the data cube
    double  *****d3; // the data cube
    double   *times; // zone times, if provided
    double      *r0; // the radius values for a 2D data set for j=0; not provided, must be computed
    double  *theta0; // the polar angle values for a 2D data set for i=0; not provided, must be computed
    double    *phi0; // The phi coordinate for 3D data sets; not provided, must be computed
    bool   timeData; // true if the zones have a "solution time" record
    bool      point; // true if input file is in point format, not block; all output files are point
};

void convertSphProj2Tecplot(struct sphericalProjection &S, struct tecplotData &o);
void writeTecplot3DIslice(char *s, struct tecplotData &d, int i0);
void tecplotHalf2FullSphere(struct tecplotData &d, struct tecplotData &e);
void findTecplotPolarCoords(struct tecplotData &d);
void readTecplot(char *s, struct tecplotData &d);
void readTecplotTitle(FILE *fp, struct tecplotData &d);
void readTecplotVarNames(FILE *fp, struct tecplotData &d);
void readTecplotDimensions(FILE *fp, struct tecplotData &d);
void countTecplotZones(FILE *fp, struct tecplotData &d);
void readTecplotData(FILE *fp, struct tecplotData &d);
void readTecplotDatapackingFormat(FILE *fp, struct tecplotData &d);
void writeTecplot(char *s, struct tecplotData &d);
void expandTo3D(struct tecplotData &d, int K);
void interpTecplotToRegularMesh(struct tecplotData &d, struct tecplotData &e, double Rmin, double Rmax, int I, int J); 
void freeTecplot(struct tecplotData &d);
int  findTecplotVarIndex(char *s, struct tecplotData &d);


