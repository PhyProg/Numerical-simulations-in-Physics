#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
//#include <zheev.h>

/* run this program using the console pauser or add your own getch, system("pause") or input loop */

/*

	This program is calculating phonon spectar for arbitrary nonpolar semiconducting material. Your task is to put appropriate
dynamical matrices on appropriate place in memory, or modify ***input function so it can read your form of input data. 
	Complex number problem is solved by defining structure cplx, containing real and imaginary parts Some functions, such as 
summing, multiplying complex numbers, calculating inner products in real and complex vector spaces, calculating exp(ic), whece 
c is real number, etc. All of this operations will be used in calculating force constants as well as dynamical matrices for 
arbitrary q vector. They are implemented as void because of optimisation. 

*/

// Constants, such as number of q vectors, number of atoms in unit cell.

#define N 512			// number of q vectors in reciprocal lattice;									//const int N = 512; 
#define n 2				// number of atoms in unit cell;												//const int n = 2;  
//#define nR 34			// number of neighbour cells force constants that will be calculated,		
						// it's now a completely random number and will be modified optionally;			//const int nR = 16;
#define pi 3.14159265359 	// pi
#define x 0				// Indexes of x,y and z components
#define y 1
#define z 2

#define A 10.187 / 18897161646.3
#define B (2.0000000000000000000*pi*18897161646.3)/10.187
#define au (2.1798741*1.e-18)/((5.2917721092*1.e-11)*(5.2917721092*1.e-11)*(1.8218779*1.e-30))
#define qu ((1.602176565 * 1.602176565)*1.e-38)
#define V (A*A*A)
#define e0 (8.854187817*1.e-12)

// Variables used ih whole program

struct cplx
{
	double Re, Im; // Real and imaginary part of complex number, respectively.
} ;

struct vect
{
	double X,Y,Z;	// Components of vector
};

struct cplx Dq[3*n][3*n]; 		// Dynamical matrix that will be diagonalized
struct cplx fR[3*n][3*n][700];	// Force constants for interaction between 0-th and nR-th cell.
double eigen_Dq[3*n];			// Eigenvalues of Dq
double modes[3*n];				// Modes for arbitrary q
double qa[3];					// Arbitrary q vector
double qv[N][3];
double R[500][3];				// Position vectors of other cells
struct vect a_[3];				// Direct lattice vectors
struct vect b_[3];				// Reciprocal lattice vectors 
int N1 = 8,N2 = 8,N3 = 8;
int nR = 500;
double Q_1[3][3], Q_2[3][3];	// Effective charge tensor for Al and As respectively
double diel[3][3];				// Dielectric tensor for AlAs
double N_Dq[3*n][3*n];			// Non-analytic part of dynamical matrix

// List of functions

// Main procedures for this program 

void input_fc(struct cplx t[3*n][3*n][N]);	// Reading input data from files
void input_q (double t[N][3]);				// Reading q vectors from file
void generate_R ();							// Generates lattice vectors
void calculating_force_constants ();		// Calculating force constants from dynamical matrices
void generate_Dq ();						// Generating dynamical matrices for arbitrary q knowing force constants
void generate_N (struct vect *ne);			// Generating neighbour
void init_lattice ();
void diagonalize_Dq ();						// Finding eigenvalues for dynamical matrix
void write_eigen_values ();					// Writting eigen values 
void write_w (int k);						// Writting modes
void write_force_constants ();				// Writting force constants to separate file
void test ();								// Testing some functions 
void sort (double *a, int num);
void q_GX ();
void q_XW ();
void q_WK ();
void q_KG ();
void q_GL ();
void input_Na();							// Reading macroscopic charge effect factors from file (Q and e(w))
void generate_NDq (double *q_vect);			// Calculating non-analytical part of dynamical matrix based on given values for
											// Born effective charge tensor, dielectric tensor and q-vector
void modify_Dq (struct cplx D[3*n][3*n][N], int l);

// Other functions used in code

struct cplx s (struct cplx a, struct cplx b);						// Sum of complex numbers
struct cplx d (struct cplx a, struct cplx b);						// Diffrence between complex numbers
struct cplx m (struct cplx a, struct cplx b);						// Multiplication of complex numbers
struct cplx conjugate (struct cplx a);								// Conjugate complex number
struct cplx e_i (double a);											// Calculating exp(ix), where x is real number
struct cplx inner_C (struct cplx *t, struct cplx *s, int num);		// Standard Inner product in complex num-dimensional vector space
double inner_R (double *t, double *s, int num);						// Standard Inner product in real num-dimensional vector space
double inner_R3 (struct vect t, struct vect s);						// Ineer product in R^3
double Real (struct cplx t);											// Real part of complex number
double Im (struct cplx t);											// Imaginary part of complex number
struct vect s_v (struct vect s, struct vect t);						// Sum of vectors
struct vect d_v (struct vect s, struct vect t);						// Diffrence between vectors
struct vect m_v (struct vect s, double t);							// Mutiplication of vectors and scalar
double norn_v (struct vect s);										// Norm of vector

// Next few lines contain usual algebraic functions in field of complex numbers, (C,+,*)

struct cplx s (struct cplx a, struct cplx b)
{
	a.Re += b.Re, a.Im += b.Im;
	return a;
}

struct cplx d (struct cplx a, struct cplx b)
{
	a.Re -= b.Re, a.Im -= b.Im;
	return a;
}

struct cplx m (struct cplx a, struct cplx b)
{
	struct cplx p;
	p.Re = a.Re * b.Re - a.Im * b.Im;
	p.Im = a.Re * b.Im + a.Im * b.Re;
	return p;
}

struct cplx conjugate (struct cplx a)
{
	a.Im *= -1;
	return a;
}

struct cplx e_i (double a)
{
	struct cplx p;
	p.Re = cos (a), p.Im = sin (a);
	return p;
}

double Real (struct cplx a)
{
	return a.Re;
}

double Im (struct cplx a)
{
	return a.Im;
}

double inner_R (double *a, double *b, int num)
{
	int i;
	double res = 0;
	for (i = 0; i < num; i++) res += *(a+i) * *(b+i);
	return res;
}

struct cplx inner_C (struct cplx *a, struct cplx *b, int num)
{
	struct cplx p;
	int i;
	p.Re = p.Im = 0;
	for (i = 0; i < num; i++)
	{
		struct cplx temp =  m(a[i],conjugate(b[i]));
		p.Re += temp.Re, p.Im += temp.Im;
	}  
	return p;
}

double inner_R3 (struct vect a, struct vect b)
{
	return a.X*b.X + a.Y*b.Y + a.Z*b.Z;
}

struct vect s_v (struct vect a, struct vect b)
{
	a.X += b.X;
	a.Y += b.Y;
	a.Z += b.Z;
	return a;
}

struct vect m_v (struct vect a, double b)
{
	struct vect ret;
	ret.X = a.X * b;
	ret.Y = a.Y * b;
	ret.Z = a.Z * b;
	return ret;
}

struct vect d_v (struct vect a, struct vect b)
{
	a.X -= b.X;
	a.Y -= b.Y;
	a.Z -= b.Z;
	return a;
}

double norm_v (struct vect a)
{
	return sqrt(inner_R3(a,a));
}

// Now functions for actual code

// This functon is calculating force constants. After declaring variables, it calls another functions to read them and then calculates
// force constants obtained using Furier transform of dynamical matrices elements.

void sort (double *a, int num)
{
	int i,j;
	for (i = 0; i < num-1; i++)
	for (j = i; j < num; j++)
	if (a[i] > a[j])
	{
		double temp;
		temp = a[j];
		a[j] = a[i];
		a[i] = temp;
	}
}

void calculating_force_constants ()
{
	// Variables for calculation dynamical matrices for arbitrary q vector, 
	// such as dynamical matrices for q in first Brilouin zone, q vectors
	// in first Brilouin zone, force constants.

	struct cplx Dyn[3*n][3*n][N]; 		// Dynamical matrices, each for N-th q vector in first Brilouin zone. 
	double q[N][3];			    		// q vectors in first Brilouin zone. ***UNITS

	// Other variables used in this function
	
	int i,j; 			 
	
	// Input parts, calling input_fc() for all data given in separate files and input_q() for file dyn0
	
	input_q (q);
	input_fc(Dyn);
	input_Na ();
	
	// Calculating force constants for q vectors in first Briloun zone, using Fourier transform.
	// Too many loops, unfortunately
	
	for (i = 0; i < 3*n; i++)
	for (j = 0; j < 3*n; j++)
	{
		int k,l;																	// Other variables
		for (k = 0; k < nR; k++)	
		{
			struct cplx temp;
			double jkloi = 0;
			fR[i][j][k].Re = 0; fR[i][j][k].Im = 0; 									// Initialisation
			for (l = 0; l < N; l++)													// Loop. For every pair (j,k) we sum all pairs (j,k) from N dynamical matrices times phase factor
			{
				qv[l][x] = q[l][x];
				qv[l][y] = q[l][y];
				qv[l][z] = q[l][z];
				
				//int i1, i2;

				/*for (i1 = 0; i1 < 3*n; i1++)
				{
				    for (i2 = 0; i2 < 3*n; i2++)
				    printf ("%le %le\t", Dyn[i1][i2][l].Re, Dyn[i1][i2][l].Im);
				    printf ("\n");
				}
				*/
				generate_NDq(q[l]);
				/*int i1, i2;
				for (i1 = 0; i1 < 3*n; i1++)
				{
				    for (i2 = 0; i2 < 3*n; i2++)
				    printf ("%le\t", N_Dq[i1][i2]);
				    printf ("\n");
				}*/
				
				modify_Dq (Dyn,l);		
				
				//int i1, i2;
/*
				for (i1 = 0; i1 < 3*n; i1++)
				{
				    for (i2 = 0; i2 < 3*n; i2++)
				    printf ("%le %le\t", Dyn[i1][i2][l].Re, Dyn[i1][i2][l].Im);
				    printf ("\n");
				}*/
						
				temp = m(Dyn[i][j][l],e_i(-1 * inner_R(q[l],R[k],3)));
				//printf ("%le %le\n", temp.Im, temp.Re);
				//printf ("D %le %le q %le %le %le\n", Dyn[i][j][l].Im, Dyn[i][j][l].Re, q[l][0],q[l][1],q[l][2]);
				
				fR[i][j][k].Re += temp.Re;
				fR[i][j][k].Im += temp.Im;			// Calculation	
				printf ("%le %le\n", temp.Re, fR[i][j][k].Im);
			}	
			//fR[i][j][k].Re /= 512; 
			//fR[i][j][k].Im /= 512;									// Also some kind of calculation
	//	printf ("(%i,%i) R = %le %le %le\n",i,j,R[k][0]*R[k][0]+R[k][1]*R[k][1]+R[k][2]*R[k][2], fR[i][j][k].Re, fR[i][j][k].Im);
	//	printf ("%le %le\n", fR[i][j][k].Re, fR[i][j][k].Im);
		}
		//printf ("%i %i\n", i,j);
	}
	//printf ("izasao\n");
}

void input_Na ()
{
	FILE *f;																// pointer to a file
	int i = 0, j;													  			// i for loop, pos for position of last character in location
	char location[] = ""; // Location of dynamical matrices on disk
	//char location[] = "";
	char line[80];															// Input buffer
	
	f = fopen (location, "r");												// Opening file
	
	while (fgets(line,80,f) && i < 50)
	{
		int k;
		if (i > 30 && i < 34)
		{
			k = i % 31;
			sscanf (line, "%le %le %le", &diel[k][0], &diel[k][1], &diel[k][2]);
		}
		if (i > 37 && i < 41)
		{
			k = i % 38;
			sscanf (line, "%le %le %le", &Q_1[k][0], &Q_1[k][1], &Q_1[k][2]);
		}
		if (i > 41 && i < 45)
		{
			k = i % 42;
			sscanf (line, "%le %le %le", &Q_2[k][0], &Q_2[k][1], &Q_2[k][2]);
		}
		i++;
	}
	/*
	printf ("diel\n");
	
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf ("%le\t", diel[i][j]);
		}
		printf ("\n");
	}
	
	printf ("q1\n");
	
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf ("%le\t", Q_1[i][j]);
		}
		printf ("\n");
	}
	
	printf ("q2\n");
	
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf ("%le\t", Q_2[i][j]);
		}
		printf ("\n");
	}
	*/
}

// Reading q vectors from file dyn0

void input_q (double q[N][3])
{
	FILE *f;																// pointer to a file
	int i = 0;													  			// i for loop, pos for position of last character in location
	char location[] = ""; // Location of dynamical matrices on disk
	//char location[] = "";
	char line[80];															// Input buffer
	
	f = fopen (location, "r");												// Opening file
	
	fgets(line,80,f);														// Two lines with information which are not necessairly for the program
	sscanf(line, "%i %i %i", &N1, &N2, &N3);								// Reading N1,N2,N3
	fgets(line,80,f);
	
	// Reading q vectors from file
	
	//printf ("q vectors\n\n");
	
	while (fgets(line,80,f) && (i<512)) 
	{
		double qx, qy, qz;
		sscanf (line, "%le %le %le", &qx, &qy, &qz);
		q[i][x] = qx*B;
		q[i][y] = qy*B;
		q[i][z] = qz*B;
		//printf ("%le %le %le \n", q[i][x], q[i][y], q[i][z]);
		i++;
	}
}

// Reading dynamical matrices from files

void input_fc(struct cplx Dyn[3*n][3*n][N])
{
	FILE *f;																  // pointer to a file
	int i, pos = 24;													  // i for loop, pos for position of last character in location, 18,51 + 6
	char location[] = ""; 								  // Location of dynamical matrices in disk
	//char location[] = "";
	for (i = 0; i < N; i++)
	{
		int j,k = -1;
		char line[80];
		f = fopen (location, "r");
		
		while (fgets (line, 80, f))
		{
			k++;
			int a;
	
			if (k > 12 && k < 16)
			{
				a = k%13;
				sscanf (line, "%le %le %le %le %le %le", &Dyn[a][x][i].Re, &Dyn[a][x][i].Im, &Dyn[a][y][i].Re, 
														 &Dyn[a][y][i].Im, &Dyn[a][z][i].Re, &Dyn[a][z][i].Im);
			}
			if (k > 16 && k < 20)
			{
				a = k%17;
				sscanf (line, "%le %le %le %le %le %le", &Dyn[a][x+3][i].Re, &Dyn[a][x+3][i].Im, &Dyn[a][y+3][i].Re, 
														 &Dyn[a][y+3][i].Im, &Dyn[a][z+3][i].Re, &Dyn[a][z+3][i].Im);
			}
			if (k > 20 && k < 24)
			{
				a = k%21 + 3;
				sscanf (line, "%le %le %le %le %le %le", &Dyn[a][x][i].Re, &Dyn[a][x][i].Im, &Dyn[a][y][i].Re, 
														 &Dyn[a][y][i].Im, &Dyn[a][z][i].Re, &Dyn[a][z][i].Im);
			}
			if (k > 24 && k < 28)
			{
				a = k%25 + 3;
				sscanf (line, "%le %le %le %le %le %le", &Dyn[a][x+3][i].Re, &Dyn[a][x+3][i].Im, &Dyn[a][y+3][i].Re, 
														 &Dyn[a][y+3][i].Im, &Dyn[a][z+3][i].Re, &Dyn[a][z+3][i].Im);
			}
		}
		/*
		if (i == 150)
		{
			printf ("\n%c%c%c\n", location[pos-2],location[pos-1],location[pos]);
		
			int j;
			for (j = 0; j<3*n; j++)
			{
				for (k = 0; k<3*n; k++)
				{
					printf ("%le   %le\t", Dyn[j][k][i].Re, Dyn[j][k][i].Im);
				}
				printf ("\n");
			}
		}
		 */
		
		fclose (f);
		
		//Change of location
	//	printf ("\n%c%c%c", location[pos-2],location[pos-1],location[pos]);
		
		if (location[pos] == 57)
		location [pos-1] == 57 ? location [pos-2]++, location[pos-1] = location[pos] = 48 : location [pos-1]++, location[pos] = 48;
		else
		location [pos]++;
	}
	
//}
}

void generate_NDq (double *q)
{
	int i,j;
	double qeq = 0, sum1 = 0, sum2 = 0;
	if (q[x] == 0. && q[y] == 0. && q[z] == 0.) 
	{
	    q[x] = 0.0000000000001*B;
	    q[y] = 0.0000000000001*B; 
	    q[z] = 0.0000000000001*B;
	}
	for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	qeq += q[i] * diel[i][j] * q[j];
	qeq *= V;
	qeq *= e0;
	//printf ("qeq %le\n", qeq);
	for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	{
		int k;
		sum1 = 0, sum2 = 0;
		for (k = 0; k < 3; k++)
		{
			sum1 += Q_1[k][i%3] * q[k];
			sum2 += Q_1[k][j%3] * q[k];
			//printf (" sume %le %le\n", sum1, sum2);
		}
		//printf (" sume %le %le\n", sum1, sum2);
		N_Dq[i][j] = qu * sum1 * sum2 / qeq;
		//printf ("%i %i %le\n", i,j, N_Dq[i][j]);
	}
	for (i = 0; i < 3; i++)
	for (j = 3; j < 6; j++)
	{
		int k;
		sum1 = 0, sum2 = 0;
		for (k = 0; k < 3; k++)
		{
			sum1 += Q_1[k][i%3] * q[k];
			sum2 += Q_2[k][j%3] * q[k];
		}
		//printf (" sume %le %le\n", sum1, sum2);
		N_Dq[i][j] = qu * sum1 * sum2 / qeq;
		//printf ("%i %i %le\n", i,j, N_Dq[i][j]);
	}
	for (i = 3; i < 6; i++)
	for (j = 0; j < 3; j++)
	{
		int k;
		sum1 = 0, sum2 = 0;
		for (k = 0; k < 3; k++)
		{
			sum1 += Q_1[k][i%3] * q[k];
			sum2 += Q_2[k][j%3] * q[k];
		}
		//printf (" sume %le %le\n", sum1, sum2);
		N_Dq[i][j] = qu * sum1 * sum2 / qeq;
		//printf ("%i %i %le\n", i,j, N_Dq[i][j]);
	}
	for (i = 3; i < 6; i++)
	for (j = 3; j < 6; j++)
	{
		int k;
		sum1 = 0, sum2 = 0;
		for (k = 0; k < 3; k++)
		{
			sum1 += Q_2[k][i%3] * q[k];
			sum2 += Q_2[k][j%3] * q[k];
		}
		//printf (" sume %le %le\n", sum1, sum2);
		N_Dq[i][j] = qu * sum1 * sum2 / qeq;
		//printf ("%i %i %le\n", i,j, N_Dq[i][j]);
	}
	
	/*for (i = 0; i < 3*n; i++)
	{
		for (j = 0; j < 3*n; j++)
		{
			printf ("%le\t", N_Dq[i][j]);
		}
		printf ("\n");
	}
*/

    int i1, i2;
	for (i1 = 0; i1 < 3*n; i1++)
	{
	   for (i2 = 0; i2 < 3*n; i2++)
	   printf ("%le\t", N_Dq[i1][i2]);
	   printf ("\n");
	}
	
}

void modify_Dq (struct cplx Dyn[3*n][3*n][N], int l) 
{
	int i,j;
	for (i = 0; i < 3*n; i++)
	for (j = 0; j < 3*n; j++)
	Dyn[i][j][l].Re -= N_Dq[i][j];
}

// Generating dynamical matrices for arbitrary q vector

void generate_Dq ()
{
	int i, j;
	for (i = 0; i < 3*n; i++)
	for (j = 0; j < 3*n; j++)
	{
		Dq [i][j].Re = Dq [i][j].Im = 0;
		int k;
		struct cplx temp;
		for (k = 0; k < nR; k++)
		{
			temp = m(fR [i][j][k],e_i (inner_R(qa,R[k],3)));
			Dq [i][j].Re += temp.Re, Dq [i][j].Im += temp.Im;
		}	
	}
	
	/*
	printf ("1\n");
	for (i = 0; i < 3*n; i++)
	{
		for (j = 0; j < 3*n; j++)
		{
			printf ("%le %le   ", Dq[i][j].Re, Dq[i][j].Im);
		}
		printf ("\n");
	}
	printf ("\n\n\n");
	printf ("2\n");
	*/
	
	generate_NDq (qa);
	for (i = 0; i < 3*n; i++)
	for (j = 0; j < 3*n; j++)
	Dq[i][j].Re += N_Dq[i][j];
	
	/*
	for (i = 0; i < 3*n; i++)
	{
		for (j = 0; j < 3*n; j++)
		{
			printf ("%le %le   ", Dq[i][j].Re, Dq[i][j].Im);
		}
		printf ("\n");
	}
	printf ("\n\n\n");
	*/
}

// Writting eigen values of dynnamical matrix for arbitrary q vectors

void write_eigen_values ()
{
	int i;
	for (i = 0; i < 3*n; i++)
		printf ("%f\n", eigen_Dq[i]);
}

// Writting modes, w(q)

void write_w (int k)
{
	int i,j;
	FILE *f;
	char location[] = "";
	//char location[] = "";
	f = fopen (location, "a");
	for (i = 0; i < 3*n; i++)
	{
		fprintf (f,"%le\t", modes[i]);
	}
	//fprintf (f, "%le\t %le\t %le \n", qa[x], qa[y], qa[z]);
	fprintf (f, "\n");
	if (k == 9) fprintf (f, "\n");
	fclose (f);
}

void diagonalize_Dq ()  
{
	double rwork[3*3*n - 2];
	double complex c_H [36], work[100];
	int i, j, n1 = n, lda = n, lwork = 100, info = 1;
	
	for (i = 0; i < 3*n; i++)
	for (j = 0; j < 3*n; j++)
		c_H [i*3*n + j] = Dq[i][j].Re + Dq[i][j].Im * I;
		
//	while (info != 0)
//	zheev_("Vectors", "Upper", &n1, c_H, &lda, eigen_Dq, work, &lwork, rwork, &info);
	
	for (i = 0; i < 3*n; i++)
	modes[i] = sqrt (eigen_Dq[i] * (eigen_Dq[i] > 0 ? 1 : -1));
	sort (modes, 6);
}

// Generating R vectors over which interpolation is done

void generate_R()
{
	struct vect neighbour[12];
	generate_N (neighbour);
	int i,j,k,l, count = 0;
	//for (i = 0; i<12;i++)
	//printf ("%le %le %le \n", neighbour[i].X, neighbour[i].Y, neighbour[i].Z);
	for (i = -8; i < 8; i++)
	for (j = -8; j < 8; j++)
	for (k = -8; k < 8; k++)
	{
		struct vect temp;
		temp = s_v (s_v(m_v(a_[0],i) , m_v(a_[1],j)) , m_v(a_[2],k));
		int m = 0;
		for (l = 0; l < 12; l++)
		{
			if (norm_v(temp) <= norm_v(d_v(temp, neighbour[l]))) m++;
			//printf ("i = %i j = %i k = %i l = %i  m = %i\n",i,j,k,l,m);
			if (m == 12 && l == 11 && i+j+k)
			{
				R[count][x] = (i*a_[0].X + j*a_[1].X + k*a_[2].X);
				R[count][y] = (i*a_[0].Y + j*a_[1].Y + k*a_[2].Y);
				R[count][z] = (i*a_[0].Z + j*a_[1].Z + k*a_[2].Z);
			//	printf ("%le %le %le\n", R[count][x], R[count][y], R[count][z]);
				count++;
			}			
		}
	}
	nR = count;
	//printf ("%i\n", nR);
}

// Help procedure that generates neighbourhood of 0th cell

void generate_N (struct vect *ne)
{
	int i;
	ne[0] = m_v (a_[0],N1);
	ne[1] = m_v (a_[1],N2);
	ne[2] = m_v (a_[2],N3);
	ne[3] = d_v (ne[0],ne[1]);
	ne[4] = d_v (ne[1],ne[2]);
	ne[5] = d_v (ne[2],ne[0]);
	for (i = 6; i < 12; i++)
	ne[i] = m_v (ne[i%6], -1);
}

// Moving from (0,0,0) to (1/2,0,1/2)

void q_GX ()
{
	int i;
	qa[x] = qa[y] = qa[z] = 0;
	for (i = 0; i < 10; i++)
	{
		qa[x] -= 0.1000 * B;
		qa[y] += 0.0000 * B;
		qa[z] += 0.0000 * B;
		//printf ("\n%le %le %le \n", qa[x],qa[y],qa[z]);
		generate_Dq ();
    	diagonalize_Dq ();
		write_w (i);
	}
}

// Moving from (1/2,0,1/2) to (1/2,1/4,3/4)

void q_XW ()
{
	int i;
	for (i = 0; i < 10; i++)
	{
		qa[x] += 0.0000 * B;
		qa[y] += 0.0500 * B;
		qa[z] += 0.0000 * B;
		//printf ("%le %le %le \n", qa[x],qa[y],qa[z]);
		generate_Dq ();
    	diagonalize_Dq ();
		write_w (i);
	}
}

void q_WK ()
{
	int i;
	for (i = 0; i < 10; i++)
	{
		qa[x] += 0.0625 * B;
		qa[y] -= 0.0125 * B;
		qa[z] += 0.0000 * B;
		//printf ("%le %le %le \n", qa[x],qa[y],qa[z]);
		generate_Dq ();
    	diagonalize_Dq ();
		write_w (i);
	}
}

void q_KG ()
{
	int i;
	for (i = 0; i < 10; i++)
	{
		qa[x] += 0.0375 * B;
		qa[y] -= 0.0375 * B;
		qa[z] += 0.0000 * B;
		//printf ("%le %le %le \n", qa[x],qa[y],qa[z]);
		generate_Dq ();
    	diagonalize_Dq ();
		write_w (i);
	}
}

// Moving from (0,0,0) to (1/2,1/2,1/2)

void q_GL ()
{
	int i;
	qa[x] = qa[y] = qa[z] = 0;
	for (i = 0; i < 10; i++)
	{
		qa[x] -= 0.0500 * B;
		qa[y] += 0.0500 * B;
		qa[z] += 0.0500 * B;
		//printf ("%le %le %le \n", qa[x],qa[y],qa[z]);
		generate_Dq ();
    	diagonalize_Dq ();
		write_w (i);
	}
}

void init_lattice ()
{
	a_[0].X = -1 * A/2;
	a_[0].Y = 0;
	a_[0].Z = A/2;
	
	a_[1].X = 0;
	a_[1].Y = A/2;
	a_[1].Z = A/2;
	
	a_[2].X = -1 * A/2;
	a_[2].Y = A/2;
	a_[2].Z = 0;
	
	b_[0].X = -1 * B;
	b_[0].Y = -1 * B;
	b_[0].Z = B;
	
	b_[1].X = B;
	b_[1].Y = B;
	b_[1].Z = B;
	
	b_[2].X = -1 * B;
	b_[2].Y = B;
	b_[2].Z = -1 * B;
}

void q_test ()
{
	int i;
	FILE *f;
	    char location[] = "";
	    //char location[] = "";
	    f = fopen (location, "a");
	for (i = 0; i < N; i++)
	{
		qa[x] = qv[i][x];
		qa[y] = qv[i][y];
		qa[z] = qv[i][z];
		generate_Dq ();
		fprintf (f, "\n\n%i\n\n", i);
		int j,k;
		for (j = 0; j < 3*n; j++)
		{
		    for (k = 0; k < 3*n; k++)
    	    fprintf (f,"%le %le\t", Dq[j][k].Re, Dq[j][k].Im);
    	    fprintf (f,"\n");
    	}
    	printf ("%i \n", i);
    	//diagonalize_Dq ();
		//write_w (7869);
	}
}

int main(int argc, char *argv[]) {
	int counter; 

	init_lattice ();
	generate_R ();
	calculating_force_constants ();
	
	// moving q vectors through first Briloun zone
	
//	q_GX ();
//	q_XW ();
//	q_WK ();
//	q_KG ();
//	q_GL ();
  q_test ();
	
	return 0;
}
