/*
 *  idealBEC.c
 *  
 *
 *  Created by yuhey IWATA.
 *  Copyright 2011 Yuhey. All rights reserved.
 *
 */

#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_zeta.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sum.h>
#include<gsl/gsl_deriv.h>

//Bose-Einstein functionとchemical potentialの定義
double BEF(double, double);
double mu(double);
double Specific_heat(double, double, double);
double Energy(double, void *);
double Number(double);


//t:Temperature
//energy : Energy
//number : Condensation number
//sf : Specific feat

int main( )
{
	double t, energy, number, sf, sf2;
	int i, N;
	FILE *fp;
	fp = fopen("BEC_ideal2.dat", "w");
	
	//partition
	N = 50;

	fprintf(fp, "#  Temp  ChemiPo  Energy  Num\n");
	
	for(i=0; i<=2*N; i++)	
	{
		t = 1.0*i/N;
		
		
		energy = Energy(t, 0);
		number = Number(t);
		sf = Specific_heat(t,2.0/N, sf2);

		
		printf("%20.16f  %20.16f  %20.16f  %20.16f  %20.16f \n", t, mu(t), energy, number, sf);
		fprintf(fp ,"%20.16f  %20.16f  %20.16f  %20.16f  %20.16f\n", t, mu(t), energy, number, sf);

		sf2 = sf;
	}

	fclose(fp);

	return 0;
	
}


double Number(double t)
{
	if(t < 1.0)
	{
		return 1.0 - pow(t, 3.0/2.0);
	}
	else
	{
		return 0.0;
	}
}

double Specific_heat(double t, double h, double sf2)
{
/*
	gsl_function F;
	double result, abserr;
	F.function = &Energy;
	F.params = 0;
	
	gsl_deriv_backward(&F, t, 1.0e-5, &result, &abserr);
	
//	printf("%f %f %f \n", result, result - 99.0/4.0 * gsl_sf_zeta(11.0/2.0)/gsl_sf_zeta(9.0/2.0) * pow(t, 9.0/2.0) , t);
	return result;	
	*/
	if(t <= 1) {
		return 15.0/4.0 * gsl_sf_zeta(5.0/2.0)/gsl_sf_zeta(3.0/2.0) * pow(t,3.0/2.0);
	}
	if(t > 1)
	{
	//	printf("%f\n", t+15.0*h);
		return sf2/2.0 + ((Energy(t+2.0*h,0) - Energy(t,0))/(2.0*h)
				+ (Energy(t+3.0*h,0) - Energy(t,0))/(3.0*h)
				+ (Energy(t+4.0*h,0) - Energy(t,0))/(4.0*h) 
				+ (Energy(t+6.0*h,0) - Energy(t,0))/(6.0*h) 
				+ (Energy(t+7.0*h,0) - Energy(t,0))/(7.0*h) 
				+ (Energy(t+9.0*h,0) - Energy(t,0))/(9.0*h)
				+ (Energy(t+10.0*h,0) - Energy(t,0))/(10.0*h) 
				+ (Energy(t+12.0*h,0) - Energy(t,0))/(12.0*h)
				+ (Energy(t+13.0*h,0) - Energy(t,0))/(13.0*h)
				+ (Energy(t+15.0*h,0) - Energy(t,0))/(15.0*h))/20.0;
	}
}

double Energy(double t, void * params)
{
	if(t < 1.0)
	{
		return 3.0/2.0 * pow(t,5.0/2.0) * gsl_sf_zeta(5.0/2.0)/gsl_sf_zeta(3.0/2.0);
	}
	else if(t > 1.0)
	{
		return 3.0/2.0 * pow(t,5.0/2.0) * BEF(5.0/2.0, t)/gsl_sf_zeta(3.0/2.0);
	}
	else
	{
		return 3.0/2.0 * gsl_sf_zeta(5.0/2.0)/gsl_sf_zeta(3.0/2.0);
	}
	
}
	
double BEF(double p, double t)
{
	double g;
	int n, N=10;
	double r[N], sum_accel, err, sum = 0.0;
	double np1;

	
	if(t < 1.0)
	{
		g = gsl_sf_zeta(p);
		printf("error!\n");
	}
	else if(t > 1.0)
	{
//		exp(mu(t)/t)^n/n^(p)need mu function
		gsl_sum_levin_u_workspace * w = gsl_sum_levin_u_alloc(N);
		for (n = 0; n < N; n++)
		{
			np1 = n + 1.0;
			r[n] = pow(exp( mu(t)/t ) ,np1) / pow(np1, p);
			sum += r[n];
		}		
		gsl_sum_levin_u_accel(r, N, w, &sum_accel, &err);
		g = sum_accel;
		gsl_sum_levin_u_free(w);
	}
	
	return g;
}




//高温領域でのmuの温度依存性を調べる関数
double mu(double t)
{
	int j, N;
	double mu1, mu2, g2=-5.0, g, g_exact;
	int N2 = 10;
	double r[N2], sum_accel, err, sum = 0.0;
	int n;
	double np1;
	
	g_exact = gsl_sf_zeta(3.0/2.0)/pow(t, 3.0/2.0);
	
	//partition
	N = 1000;
	
	for(j=0; j<=6.5*N; j++)
	{
		mu1 = 0.0 - 1.0*j/N;

		gsl_sum_levin_u_workspace * w = gsl_sum_levin_u_alloc(N2);
		for (n = 0; n < N2; n++)
		{
			np1 = n + 1.0;
			r[n] = pow(exp( mu1/t ) ,np1) / pow(np1, 3.0/2.0);
			sum += r[n];
		}		
		gsl_sum_levin_u_accel(r, N2, w, &sum_accel, &err);
		g = sum_accel;
		gsl_sum_levin_u_free(w);
		if( fabs( g -  g_exact ) < fabs( g2 - g_exact ) )
		{
			//printf("%f %f %20.16f %20.16f %20.16f\n",t, mu1, g, g2, g_exact);
		}else
		{
			//printf("%f %f %20.16f %20.16f %20.16f\n",t, mu1, g, g2, g_exact);
			break;
		}
		g2 = g;
		mu2 = mu1;

	} 	
	return mu2;
}

