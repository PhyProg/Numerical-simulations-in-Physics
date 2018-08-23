#include <iostream>
#include <cstdio>
#include <cmath>

/* run this program using the console pauser or add your own getch, system("pause") or input loop */

double D[10][10];
double l[10], D1[10][10], D2[10][10];
double k = 0.01;

void diag ()
{
	for (int i = 0; i<10; i++)
	{
		while (Det(D1)*Det(D2) > 0)
		{
			l[i] += k;
			for (int j = 0; j<10; j++);
			{
				D1[j][j] = D[j][j] - l[i];
				D2[j][j] = D[j][j] - 2*l[i];
			}
			if (l[i] > 100.0) return;
		}
		l[i+1]=l[i]+k;
	}
	return;
}

int main(int argc, char** argv) {

	for (int i = 0; i<10; i++)
	for (int j = 0; j<10; j++)
	{
		D[i][j] = ((i + j)>10 ? (i+j) : (i+j-10));
		D1[i][j] = D2[i][j] = D[i][j];
	}


	for (int i = 0; i++<10;)
	{
		for (int j = 0; j++<10;)
		printf ("%i ",D[i][j]);
		printf ("\n");
		l[i] = -100.0;
	}

	diag();



	return 0;
}
