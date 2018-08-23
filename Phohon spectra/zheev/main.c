#include <stdio.h>
#include <stdlib.h>
#include <zheev.h>

/* run this program using the console pauser or add your own getch, system("pause") or input loop */

struct _dcomplex { double re, im; };
//typedef struct _dcomplex doublecomplex;


doublecomplex herm[9] = 
{
	{1,1}, {1,1}, {1,1},
	{1,1}, {1,1}, {1,1},
	{1,1}, {1,1}, {1,1}
};

doublereal eigen[3], rwork[3];
doublecomplex work[3];

int main(int argc, char *argv[]) {
	integer i, n = 3, lda = 3, lwork, info;
	//for (i = 0; i < 9; i++)
	//printf ("%le %le\n", herm[i].re, herm[i].im);
	while (!info)
	zheev_("N", "U", &n, herm, &lda, eigen, work, &lwork, rwork, &info);
	
	for (i = 0; i < 3; i++)
	printf ("%le ", eigen[i]);
	
	return 0;
}
