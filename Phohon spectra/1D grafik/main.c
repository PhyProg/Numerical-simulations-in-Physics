#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.141592
#define mr m1*m2/(m1+m2)
/* run this program using the console pauser or add your own getch, system("pause") or input loop */

double m1,m2,k,q,wa, wo,a;

int i;
int main(int argc, char *argv[]) {
	m1 = 1;
	m2 = 5;
	k  = 15;
	a  = 1;
	q  = -pi/a;
	FILE *f1 = fopen ("", "r");
	FILE *f = fopen ("", "w");
//	while (q <= pi/a)
//	{
		//wa = sqrt ((k/mr)*(1-sqrt(1-mr/(m1*m2) * sin(q*a/2)*sin(q*a/2))));
		//wo = sqrt ((k/mr)*(1+sqrt(1-mr/(m1*m2) * sin(q*a/2)*sin(q*a/2))));
		for (i = 0; i < 62; i++)
		{
				fscanf (f1, "%le %le", &wa, &wo);
					wa = sqrt (wa);
		wo = sqrt (wo);
		fprintf (f, "%le \t %le\n", wa,wo);
		}
	

		q += 0.05*pi/a;
	//}
	fclose (f);
	fclose (f1);
	return 0;
}
