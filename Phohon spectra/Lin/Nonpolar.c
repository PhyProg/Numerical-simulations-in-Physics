#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
#define nR 16			// number of neighbour cells force constants that will be calculated,		
						// it's now a completely random number and will be modified optionally;			//const int nR = 16;
#define pi 3.141592 	// pi
#define x 0				// Indexes of x,y and z components
#define y 1
#define z 2

#define A 1
#define B 1

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
struct cplx fR[3*n][3*n][nR];	// Force constants for interaction between 0-th and nR-th cell.
double eigen_Dq[3*n];			// Eigenvalues of Dq
double modes[3*n];				// Modes for arbitrary q
double qa[3];					// Arbitrary q vector
double R[nR][3];				// Position vectors of other cells
struct vect a_[3];				// Direct lattice vectors
struct vect b_[3];				// Reciprocal lattice vectors 
int N1,N2,N3;

// List of functions

// Main procedures for this program 

void input_fc(struct cplx t[3*n][3*n][N]);	// Reading input data from files
void input_q (double t[N][3]);				// Reading q vectors from file
void generate_R ();							// Generates lattice vectors
void calculating_force_constants ();		// Calculating force constants from dynamical matrices
void generate_Dq ();						// Generating dynamical matrices for arbitrary q knowing force constants
void generate_N (struct vect *ne);			// Generating neighbour
void diagonalize_Dq ();						// Finding eigenvalues for dynamical matrix
void write_eigen_values ();					// Writting eigen values 
void write_w ();							// Writting modes
void write_force_constants ();				// Writting force constants to separate file
void test ();								// Testing some functions 
void q_GX ();
void q_XW ();
void q_WK ();
void q_KG ();
void q_GL ();

// Other functions used in code

struct cplx s (struct cplx a, struct cplx b);						// Sum of complex numbers
struct cplx d (struct cplx a, struct cplx b);						// Diffrence between complex numbers
struct cplx m (struct cplx a, struct cplx b);						// Multiplication of complex numbers
struct cplx conjugate (struct cplx a);								// Conjugate complex number
struct cplx e_i (double a);											// Calculating exp(ix), where x is real number
struct cplx inner_C (struct cplx *t, struct cplx *s, int num);		// Standard Inner product in complex num-dimensional vector space
double inner_R (double *t, double *s, int num);						// Standard Inner product in real num-dimensional vector space
double inner_R3 (struct vect t, struct vect s);						// Ineer product in R^3
double Re (struct cplx t);											// Real part of complex number
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

double Re (struct cplx a)
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
	for (i = 0; i < num; i++) res += *(a+i) * *(b+1); 
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
	a.X *= b;
	a.Y *= b;
	a.Z *= b;
	return a;
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

void calculating_force_constants ()
{
	// Variables for calculation dynamical matrices for arbitrary q vector, 
	// such as dynamical matrices for q in first Brilouin zone, q vectors
	// in first Brilouin zone, force constants.

	struct cplx Dyn[3*n][3*n][N]; 		// Dynamical matrices, each for N-th q vector in first Brilouin zone. 
	double q[N][3];			    		// q vectors in first Brilouin zone. ***UNITS
	double R[nR][3];					// Positions of other cells relative to 0th cell

	// Other variables used in this function
	
	int i,j; 			 
	
	// Input parts, calling input_fc() for all data given in separate files and input_q() for file dyn0
	
	input_q (q);
	input_fc(Dyn);
	
	// Calculating force constants for q vectors in first Briloun zone, using Fourier transform.
	// Too many loops, unfortunately
	
	for (i = 0; i < 3*n; i++)
	for (j = 0; j < 3*n; j++)
	{
		int k,l;																	// Other variables
		for (k = 0; k < nR; k++)	
		{
			struct cplx temp;
			fR[i][j][k].Re = fR[i][j][k].Im = 0; 									// Initialisation
			for (l = 0; l < N; l++)													// Loop. For every pair (j,k) we should sum all pairs (j,k) from N dynamical matrices times phase factor
			{
				
				temp = m(Dyn[i][j][l],e_i(-1 * inner_R(q[l],R[k],3)));
				fR [i][j][k].Re += temp.Re, fR[i][j][k].Im += temp.Im ;				// Calculation	
			//	printf ("%i %i %i %i\n", i,j,k,l);
			}	
		fR[i][j][k].Re /= N, fR[i][j][k].Im /= N;									// Also some kind of calculation
		}
	}
}

// Reading q vectors from file dyn0

void input_q (double q[N][3])
{
	FILE *f;																// pointer to a file
	int i = 0;													  			// i for loop, pos for position of last character in location
	char location[] = ""; // Location of dynamical matrices in disk
	char line[80];															// Input buffer
	
	f = fopen (location, "r");												// Opening file
	
	fgets(line,80,f);														// Two lines with information which are not necessairly for the program
	sscanf(line, "%i %i %i", &N1, &N2, &N3);
	fgets(line,80,f);
	
	// Reading q vectors from file
	
	while (fgets(line,80,f) && (i<512)) 
	{
		sscanf (line, "%le %le %le", q[i]+x,q[i]+y,q[i]+z);
		q[i][x] *= B;
		q[i][y] *= B;
		q[i][z] *= B;
		i++;
	}
}

// Reading dynamical matrices from files

void input_fc(struct cplx Dyn[3*n][3*n][N])
{
	FILE *f;																  // pointer to a file
	int i = 0, pos = 18;													  // i for loop, pos for position of last character in location
	char location[] = ""; // Location of dynamical matrices in disk

	for (i = 1; i++ < N;)
	{
		int j,k = 0;
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
			//printf ("\n");
		}
		
		fclose (f);
		
		//Change of location
		
		if (location[pos] == 57)
		location [pos-1] == 57 ? location [pos-2]++, location[pos-1] = location[pos] = 48 : location [pos-1]++, location[pos] = 48;
		else
		location [pos]++;
		//printf ("%s \n", location); //Check
	}
	
//}
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
}

// Writting eigen values of dynnamical matrix for arbitrary q vectors

void write_eigen_values ()
{
	int i;
	for (i = 0; i < 3*n; i++)
		printf ("%f\n", eigen_Dq[i]);
}

// Writting modes, w(q)

void write_w ()
{
	int i,j;
	FILE *f;
	char location[] = "";
	f = fopen (location, "a");
	for (i = 0; i < 3*n; i++)
	{
		fprintf (f,"%le\t", modes[i]);
	}
	fprintf (f, "%le\t %le\t %le \n", qa[x], qa[y], qa[z]);	
}

void diagonalize_Dq ()  
{
	doublecomplex c_H [3*n * 3*n], work[n];
	integer n1 = n, lda, lwork[n], info;
	double rwork;
	int i, j;
	
	for (i = 0; i < 3*n; i++)
	for (j = 0; j < 3*n; j++)
	{
		c_H[3*i + j].re = Dq[i][j].Re;
		c_H[3*i + j].im = Dq[i][j].Im;
	}
	while (!info)
	zheev_("N", "U", &n1, c_H, &lda, modes, work, lwork, &rwork, &info);
	
	for (i = 0; i < 3*n; i++)
	modes[i] = sqrt (eigen_Dq[i] * (eigen_Dq[i] > 0 ? 1 : -1));
}

// Generating R vectors over which interpolation is done

void generate_R()
{
	struct vect neighbour[12];
	generate_N (neighbour);
	int i,j,k,l, count = 0;
	for (i = 0; i < 5; i++)
	for (j = 0; j < 5; j++)
	for (k = 0; k < 5; k++)
	{
		struct vect temp;
		temp = s_v (s_v(m_v(a_[0],i) , m_v(a_[1],j)) , m_v(a_[2],k));
		for (l = 0; l < 12; l++)
		{
			int m = 0;
			if (norm_v(temp) < norm_v(d_v(temp, neighbour[l]))) m++;
			if (m == 12)
			{
				R[count][x] = (i*a_[0].X + j*a_[1].X + k*a_[2].X) * A;
				R[count][y] = (i*a_[0].Y + j*a_[1].Y + k*a_[2].Y) * A;
				R[count][z] = (i*a_[0].Z + j*a_[1].Z + k*a_[2].Z) * A;
				count++;
			}			
		}
	}
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
	for (i = 0; i++ < 10;)
	{
		qa[x] += 0.0500 * B;
		qa[y] += 0.0000 * B;
		qa[z] += 0.0500 * B;
		generate_Dq ();
    	diagonalize_Dq ();
		write_w ();
	}
}

// Moving from (1/2,0,1/2) to (1/2,1/4,3/4)

void q_XW ()
{
	int i;
	for (i = 0; i++ < 10;)
	{
		qa[x] += 0.0000 * B;
		qa[y] += 0.0075 * B;
		qa[z] += 0.0075 * B;
		generate_Dq ();
    	diagonalize_Dq ();
		write_w ();
	}
}

void q_WK ()
{
	int i;
	for (i = 0; i++ < 10;)
	{
		
		generate_Dq ();
    	diagonalize_Dq ();
		write_w ();
	}
}

void q_KG ()
{
	int i;
	for (i = 0; i++ < 10;)
	{
		
		generate_Dq ();
    	diagonalize_Dq ();
		write_w ();
	}
}

// Moving from (0,0,0) to (1/2,1/2,1/2)

void q_GL ()
{
	int i;
	qa[x] = qa[y] = qa[z] = 0;
	for (i = 0; i++ < 10;)
	{
		qa[x] += 0.0500 * B;
		qa[y] += 0.0500 * B;
		qa[z] += 0.0500 * B;
		generate_Dq ();
    	diagonalize_Dq ();
		write_w ();
	}
}

int main(int argc, char *argv[]) {
	int counter; 
	
	generate_R ();
	calculating_force_constants ();
	
	// moving q vectors through first Briloun zone
	
	q_GX ();
	q_XW ();
	q_WK ();
	q_KG ();
	q_GL ();
	
	return 0;
}
