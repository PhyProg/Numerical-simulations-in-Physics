#include <iostream>
#include <cstdio>
#include <cmath>

#define N 31
#define n_c 1
/* run this program using the console pauser or add your own getch, system("pause") or input loop */

double b1, b2, n_O2[N], n_N2[N], n_N[N];
double T_n[N], T_e, Beta[21][N], h[N];
double C_i[4], D_i[4], Alfa[21][3][N];
double N_e[10000][N], dN_e[10000][N], N_e1[10][N], b[10000], h1[10000], t_h[10000], t_m[10000], t_s[10000], test, Te[10000][N];
double q1[10000][N], bet[10000][N], alf[10000][N], betN[10000][N], alfN[10000][N];
int DT = 1, I = 0;
FILE *f_IN, *f_OUT, *f_OUTa, *f_INn, *f_OUTn;

//Konstante

double B1c = 1.0e-43, B2c1 = 1.1617e-39, B2c2 = 3.4665e-42, B2c3 = 3.2825e-45;
double B2a1,B2a2,B2a3;

double beta(double T_1, double T_2, double n)
{
	//T_1 temperatura elektrona, T_2 temperatura N, n koncentracija O2
	//printf ("%f %f %f\n", T_1,T_2,n);
	B2a1 = 7.8193e2 - 3.2964*T_2;
	B2a2 = -1.9159e2+ 3.7646*T_2 - 4.5446e-3 * T_2*T_2;
	B2a3 = -7.6834e1+ 1.2277e-2*T_2 - 7.6427e-3 * T_2*T_2 + 1.7856e-5 * T_2*T_2*T_2;
	b2 = (B2c1 - B2c2*T_2 + B2c3*T_2*T_2) * pow(1/T_1, 0.65) * exp (-B2a1/T_1 - pow (B2a2/T_1,2) - pow (B2a3/T_1,3)) * n*n; //(22)
	//printf ("%.50f \n %.50f\n", b1,b2);
	return b2+b1;
}

double q (double T)
{
	double Q;

	return Q;
}

double C (double T)
{
	double c;

	return c;
}

double D (double T,double Ne, int i)
{
	double d;
	d = q(T)-beta(T, T_n[i],n_O2[i])*N-C(T)*Ne*Ne*sqrt(300)/sqrt(T);
	return d;
}

double abs (double a, double b)
{
	return a-b > 0 ? a-b : b-a;
}

void ucitavanje ()
{
	char line[38], line1[154];
    f_IN = fopen ("", "r");
    f_INn = fopen ("", "r");
    f_OUT = fopen ("", "w");
    f_OUTa = fopen ("", "w");
    f_OUTn = fopen ("", "w");
    for (int i = 0; i < N; i++)
	{
		fgets (line, 38, f_IN);
		//printf ("Linija : %s", line);
		sscanf (line,"%le %le %le %le\n",h+i, n_N2+i, n_O2+i, T_n+i);
		//printf ("Vrednosti : %f %f %f %f\n",h[i], n_N2[i], n_O2[i], T_n[i]);
	}
	/*for (int i = 0; i < N; i++)
	{
		fgets (line, 38, f_INn);
		//printf ("Linija : %s", line);
		sscanf (line,"%le %le\n", &N_e[0][i], &N_e[1][i]);
		//printf ("Vrednosti : %f %f %f %f\n",h[i], n_N2[i], n_O2[i], T_n[i]);
	}*/
	C_i[0] = 1.8e-13;
	C_i[1] = 1.6e-13;
	C_i[2] = 4.5e-13;
	C_i[3] = (5+20*n_c)*(1e-13);
	D_i[0] = -0.39;
	D_i[1] = -0.55;
	D_i[2] = -0.83;
	D_i[3] = -0.5;

	fgets (line1, 154, f_INn);

	//printf ("%s\n", line1);

	while (fgets (line1, 154, f_INn))
	{
		sscanf (line1, "%le %le %le %le %le %le %le %le %le", t_h+I, t_m+I, t_s+I, &test, &test, &test, &test, b+I, h1+I);

		//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);

		for (int i = 0; i < N; i++)
		{
			N_e[I][i] = 1.43e13 * exp (-b[I]*h1[I]) * exp ((b[I]-0.15)*h[i]);
			if ((t_h[I] == 12.000000) && (t_m[I] == 00.000000) && (t_s[I] ==  0.0000000))
			{
				//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);
				N_e1[0][i] = N_e[I][i];
			}
			if ((t_h[I] == 12.000000) && (t_m[I] == 10.000000) && (t_s[I] ==  0.0000000))
			{
				//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);
				N_e1[1][i] = N_e[I][i];
			}
			if ((t_h[I] == 12.000000) && (t_m[I] == 20.000000) && (t_s[I] ==  0.0000000))
			{
				//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);
				N_e1[2][i] = N_e[I][i];
			}
			if ((t_h[I] == 12.000000) && (t_m[I] == 30.000000) && (t_s[I] ==  0.0000000))
			{
				//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);
				N_e1[3][i] = N_e[I][i];
			}
			if ((t_h[I] == 12.000000) && (t_m[I] == 35.000000) && (t_s[I] ==  0.0000000))
			{
				//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);
				N_e1[4][i] = N_e[I][i];
			}
			if ((t_h[I] == 12.000000) && (t_m[I] == 40.000000) && (t_s[I] ==  0.0000000))
			{
				//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);
				N_e1[5][i] = N_e[I][i];
			}
			if ((t_h[I] == 12.000000) && (t_m[I] == 44.000000) && (t_s[I] ==  0.0000000))
			{
				//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);
				N_e1[6][i] = N_e[I][i];
			}
			if ((t_h[I] == 12.000000) && (t_m[I] == 47.000000) && (t_s[I] ==  0.0000000))
			{
				//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);
				N_e1[7][i] = N_e[I][i];
			}
			if ((t_h[I] == 12.000000) && (t_m[I] == 48.000000) && (t_s[I] ==  30.000000))
			{
				//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);
				N_e1[8][i] = N_e[I][i];
			}
			if ((t_h[I] == 12.000000) && (t_m[I] == 50.000000) && (t_s[I] ==  0.0000000))
			{
				//printf ("%f %f %f\n", t_h[I], t_m[I], t_s[I]);
				N_e1[9][i] = N_e[I][i];
			}
		}
		I++;
	}

}

void racunanje_beta ()
{
	for (DT = 0; DT < 21; DT++)
    {
    	//fprintf(f_OUT, "DT = %i\n\n\t h\t Tn \t Dbeta \n", DT);
    	fprintf(f_OUT, "\n\nDT = %i\n\n\t h \t Beta \n", DT);
    	//for (int i = 0; i<N; i++)
    	//fprintf(f_OUT,"%f \n", T_n[i]);
    	for (int i = 0; i<N; i++)
    	{
    		b1 = B1c * n_O2[i] * n_N2[i];
    		Beta[DT][i] = beta(T_n[i] + DT, T_n[i], n_O2[i]);// - beta(T_n[i],T_n[i],n_O2[i]);
    	//	printf ("DBeta = %f\t B = %f \t B = %f \n", DBeta, beta(T_n[i] + DT, T_n[i], n_O2[i]), beta(T_n[i],T_n[i],n_O2[i]));
    	//    fprintf(f_OUT,"\t %2.f %f %.20f \n", h[i] ,T_n[i], DBeta);
    	//    fprintf(f_OUT,"\t%f %.20f \n",T_n[i], DBeta);
		    fprintf(f_OUT,"%f \t %.20f \n", h[i], Beta[DT][i]);
		}
	}
	for (DT = 0; DT < 21; DT++)
	{
		fprintf(f_OUT, "\n\nDT = %i\n\n\t Tn \t Dbeta \n", DT);
		for (int i = 0; i<N; i++)
		{
			double DBeta = Beta[DT][i]-Beta[0][i];
			fprintf(f_OUT,"%f \t %.20f \n", h[i], DBeta/Beta[0][i]);
		}
	}
}

void racunanje_alfa ()
{
	for(int i = 0; i<3; i++)
	{
		fprintf(f_OUTa, "\n\nElement = %i\n\n", i);
		for(DT = 0; DT<21; DT++)
		{
			fprintf(f_OUTa, "\n\n DT = %i\n\n", DT);
			fprintf(f_OUTa, "h \t Alfa\n");
			for(int j = 0; j<N; j++)
			{
				T_e = T_n[j]+DT;
				Alfa[DT][i][j] = C_i[i] * pow(T_e/300, D_i[i]);
				//printf ("%le \n", Alfa[i][j]);
				fprintf(f_OUTa, "%le \t %le\n", h[j] , Alfa[DT][i][j]);
			}
//			fprintf(f_OUTa, "Te\n");
/*			for(int j = 0; j<N; j++)
			{
				fprintf(f_OUTa, "%le\n", T_n[j]+DT);
			}*/
		}
		for(DT = 0; DT<21; DT++)
		{
			fprintf(f_OUTa, "\n\n DT = %i\n\n", DT);
			fprintf(f_OUTa, "h \t DAlfa\n");
			for(int j = 0; j<N; j++)
			{
				double DAlfa = Alfa[DT][i][j]-Alfa[0][i][j];
				fprintf(f_OUTa,"%f \t %.20f \n", h[j], DAlfa/Alfa[0][i][j]);
			}
		}
	}

	//deo za alfa efektivno

	/*

	double ae;

	for (int i=0; i<N; i++)
	{
    	ae =
    }

	*/
}

void racunanje_N ()
{
	/*fprintf (f_OUTn,"test\n\n\n");
	for (int i = 0; i < N; i++)
	{
		fprintf (f_OUTn,"%f %f\n", N_e1[][i], N_e1[j][i]);
	}
	*/
	for (int j = 0; j < 10; j++)
	{
		fprintf (f_OUTn, "\n\n %i\n\n",j);
		for (int i = 0; i < N; i++)
		{
			fprintf (f_OUTn, "\n\nDN/N \n\n");
			fprintf (f_OUTn,"%f\n", (N_e1[j][i]-N_e1[9][i])/N_e1[9][i]);
			fprintf (f_OUTn, "\n\n DN2/N2 \n\n");
			fprintf (f_OUTn,"%f\n", (N_e1[j][i]*N_e1[j][i]-N_e1[9][i]*N_e1[9][i])/(N_e1[9][i]*N_e1[9][i]));
		}
	}
	/*
	for (int i = 0; i < N; i++)
	{
		dN_e[0][i] = N_e[1][i]-N_e[0][i];
		dN_e[I-1][i] = N_e[I][i]-N_e[I-1][i];
	}

	for (int j = 1; j < I-1; j++)
	{
		for (int i = 0; i < N; i++)
		{
			dN_e[j][i] = (N_e[j+1][i]-N_e[j-1][i])/2;
		}
	}
	*/
}

void racunanje_T ()
{
	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<9000; j++)
		{
			for (int k = 0; k<20; k++)
			{
				Te[j][i] = abs (D(T_n[i]+k, N_e[j][i],i),dN_e[j][i]) < abs (D(T_n[i]+k+1, N_e[j][i],i),dN_e[j][i]) ? T_n[i]+k : T_n[i]+k+1;
			}
		}
	}
}

void racunanje_O ()
{
	//q
	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<9000; j++)
		{
			q1[j][i] = q(Te[j][i]);
			bet[j][i] = beta(Te[j][i], T_n[i], n_O2[i]);
			alf[j][i] = C(Te[j][i]) * sqrt (300/Te[j][i]);
			betN[j][i] = bet[j][i] * N_e[j][i];
			alfN[j][i] = alf[j][i] * N_e[j][i] * N_e[j][i];
		}
	}
}

int main(int argc, char** argv) {

    ucitavanje();
    racunanje_beta();
    racunanje_alfa ();
    racunanje_N ();
    racunanje_T ();
    racunanje_O ();

	return 0;
}
