/*
 *  simpleBEC.c
 *  
 *
 *  Created by yuhey IWATA on 11/07/07.
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
//number : Condensate particle number
//sh : Specific heat

int main( )
{
	double t=0.0, h, energy, number, sh, sh2=4.5;
	int i;
	FILE *fp;
	fp = fopen("BEC7_2.dat", "w");
	
	//step
	h = 0.01;

	fprintf(fp, "#  Temp  ChemiPo  Energy  Num\n");
	
	while (t <= 2.01)
	{
		energy = Energy(t, 0);
		number = Number(t);
		sh = Specific_heat(t, h, sh2);

		
		printf("%20.16f  %20.16f  %20.16f  %20.16f  %20.16f \n", t, mu(t), energy, number, sh);
		fprintf(fp ,"%20.16f  %20.16f  %20.16f  %20.16f  %20.16f\n", t, mu(t), energy, number, sh);

		sh2 = sh;
		if (t>=0.9 && t<=1.01) {
			t += h*(1.01 - pow(cos(t-1.0),10));
		}
		else
		{
			t += h;
		}
	}

	fclose(fp);

	return 0;
	
}


double Number(double t)
{
	if(t < 1.0)
	{
		return 1.0 - pow(t, 9.0/2.0);
	}
	else
	{
		return 0.0;
	}
}

double Specific_heat(double t, double h, double sh2)
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
		return 99.0/4.0 * gsl_sf_zeta(11.0/2.0)/gsl_sf_zeta(9.0/2.0) * pow(t,9.0/2.0);
	}
	if(t > 1)
	{
	//	printf("%f\n", t+15.0*h);
		return sh2/6.0 + 5.0*((Energy(t+1.0*h,0) - Energy(t,0))/(1.0*h)
					   +  (Energy(t+3.0*h,0) - Energy(t,0))/(3.0*h))/12.0 ;
	}
}

double Energy(double t, void * params)
{
	if(t < 1.0)
	{
		return 9.0/2.0 * pow(t,11.0/2.0) * gsl_sf_zeta(11.0/2.0)/gsl_sf_zeta(9.0/2.0);
	}
	else if(t > 1.0)
	{
		return 9.0/2.0 * pow(t,11.0/2.0) * BEF(11.0/2.0, t)/gsl_sf_zeta(9.0/2.0);
	}
	else
	{
		return 9.0/2.0 * gsl_sf_zeta(11.0/2.0)/gsl_sf_zeta(9.0/2.0);
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
	
	g_exact = gsl_sf_zeta(9.0/2.0)/pow(t, 9.0/2.0);
	
	//partition
	N = 10000;
	
	for(j=0; j<=15.0*N; j++)
	{
		mu1 = 0.0 - 1.0*j/N;

		gsl_sum_levin_u_workspace * w = gsl_sum_levin_u_alloc(N2);
		for (n = 0; n < N2; n++)
		{
			np1 = n + 1.0;
			r[n] = pow(exp( mu1/t ) ,np1) / pow(np1, 9.0/2.0);
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

