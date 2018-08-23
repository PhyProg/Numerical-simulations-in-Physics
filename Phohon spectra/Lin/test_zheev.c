#include <stdio.h>
#include <stdlib.h>
#include <zheev.h>
#include <complex.h>
//#include <blaswrap.h>
#include <f2c.h>

/* run this program using the console pauser or add your own getch, system("pause") or input loop */

struct _dcomplex { double re, im; };
//typedef struct _dcomplex doublecomplex;


double complex herm[9] = 
{
   9.14 +  0.00*I, -4.37 - 9.22*I, -1.98 - 1.72*I, -8.96 - 9.50*I,
  -4.37 +  9.22*I, -3.35 + 0.00*I,  2.25 - 9.51*I,  2.57 + 2.40*I,
  -1.98 +  1.72*I,  2.25 + 9.51*I, -4.82 + 0.00*I, -3.24 + 2.04*I,
  -8.96 +  9.50*I,  2.57 - 2.40*I, -3.24 - 2.04*I,  8.44 + 0.00*I
};

doublereal eigen[3], rwork[3];
doublecomplex work[3];
int i, n = 3, lda = 3, lwork, info;

int main(int argc, char *argv[]) {
	
	//for (i = 0; i < 9; i++)
	//printf ("%le %le\n", herm[i].re, herm[i].im);
	while (!info)
	zheev("Vectors", "Upper", &n, herm, &lda, eigen, work, &lwork, rwork, &info);
	
	for (i = 0; i < 3; i++)
	printf ("%le ", eigen[i]);
	
	return 0;
}
