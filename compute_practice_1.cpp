// compute_practice.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double g_u[502][502];
double g_l[502][502];
double g_mat_a[502], g_mat_b, g_mat_c;

template<class T>inline T max(const T& a, const T& b)
{
	return a < b ? b : a;	
}
template<class T>inline T min(const T& a, const T& b)
{
	return a <= b ? a : b;
}

void init()
{
	int i;
	g_mat_b = 0.16;
	g_mat_c = -0.064;
	for (i = 1; i <= 501; i++)
		g_mat_a[i] = (1.64 - 0.024*i)*sin(0.2*i) - 0.64*exp(0.1 / i);
}

inline double getValue(int i, int j)
{
	if (i == j)
	{
		return g_mat_a[i];
	}
	else if (abs(i - j) == 1)
	{
		return g_mat_b;
	}
	else if (abs(i - j) == 2)
	{
		return g_mat_c;
	}
	else
		return 0;
}

double powerMethod(double offset)
{
	int i;
	double u[502], y[502];
	double beta_k = 0, beta_k_1 = 0, n = 0;
	for (i = 1; i <= 501; i++)
	{
		u[i] = 1;
		y[i] = 0;
	}

	for (int k = 1; k <= 10000; k++)
	{
		n = 0;
		for (i = 1; i <= 501; i++)
			n += u[i] * u[i];
		n = sqrt(n);
		for (i = 1; i <= 501; i++)
			y[i] = u[i] / n;
		for (i = 1; i <= 501; i++)
		{
			u[i] = 0;
			for (int j = 1; j <= 501; j++)
				u[i] += ((i == j) ? (getValue(i, j) - offset) : getValue(i, j))*y[j];
		}
		beta_k_1 = beta_k;
		beta_k = 0;
		for (i = 1; i <= 501; i++)
			beta_k += y[i] * u[i];
		if (k > 1 && fabs((beta_k_1 - beta_k) / (beta_k)) <= 1e-12)
			break;
	}
	return (beta_k + offset);
}

void LUDivision(double offset)
{
	int i, k, j, t;
	double sum;
	for (k = 1; k <= 501; k++)
	{
		for (j = 1; j <= 501; j++)
		{
			g_u[k][j] = 0; 
			g_l[k][j] = 0;
			if (k == j)
				g_l[k][j] = 1;
		}
	}
	for (k = 1; k <= 501; k++)
	{
		for (j = k; j <= min(k + 2, 501); j++)
		{
			sum = 0;
			for (t = max(1, max(k - 2, j - 2)); t <= (k - 1); t++)
				sum += g_l[k][t] * g_u[t][j];
			g_u[k][j] = ((k == j) ? (getValue(k, j) - offset) : getValue(k, j)) - sum;
		}
		if (k == 501)
			continue;
		for (i = k + 1; i <= min(k + 2, 501); i++)
		{
			sum = 0;
			for (t = max(1, max(i - 2, k - 2)); t <= (k - 1); t++)
				sum += g_l[i][t] * g_u[t][k];
			g_l[i][k] = (((i == k) ? (getValue(i, k) - offset) : getValue(i, k)) - sum) / g_u[k][k];
		}
	}
}

void solve(double x[], double b[])
{
	int i, t;
	double y[502];
	double sum;
	y[1] = b[1];
	for (i = 2; i <= 501; i++)
	{
		sum = 0;
		for (t = max(1, i - 2); t < i; t++)
			sum += g_l[i][t] * y[t];
		y[i] = b[i] - sum;
	}
	x[501] = y[501] / g_u[501][501];
	for (i = 500; i >= 1; i--)
	{
		sum = 0;
		for (t = i + 1; t <= min(i + 2, 501); t++)
			sum += g_u[i][t] * x[t];
		x[i] = (y[i] - sum) / g_u[i][i];
	}
}

double inversePowerMethod(double offset)
{
	int i;
	double u[502], y[502];
	double beta_k = 0, beta_k_1 = 0, n = 0;
	LUDivision(offset);
	for (i = 1; i <= 501; i++)
	{
		u[i] = 1;
		y[i] = 0;
	}
	for (int k = 1; k <= 10000; k++)
	{
		n = 0;
		for (i = 1; i <= 501; i++)
			n += u[i] * u[i];
		n = sqrt(n);
		for (i = 1; i <= 501; i++)
			y[i] = u[i] / n;
		solve(u, y);
		beta_k_1 = beta_k;
		beta_k = 0;
		for (i = 1; i <= 501; i++)
			beta_k += y[i] * u[i];
		beta_k = 1 / beta_k;
		if (k > 1 && fabs((beta_k_1 - beta_k) / (beta_k)) <= 1e-12)
			break;
	}
	return (beta_k + offset);
}

int main()
{
	int i, k;
	double lambda_temp_1, lambda_temp_2, lambda_1, lambda_501, lambda_s, mig_u[40], det;
	double lambdai[40];
	init();
	
	lambda_temp_1 = powerMethod(0);
	lambda_temp_2 = powerMethod(lambda_temp_1);

	lambda_1 = min(lambda_temp_1, lambda_temp_2);
	lambda_501 = max(lambda_temp_1, lambda_temp_2);	
	lambda_s = inversePowerMethod(0);

	det = 1;
	for (i = 1; i <= 501; i++)
		det *= g_u[i][i];
	for (k = 1; k <= 39; k++)
	{
		mig_u[k] = lambda_1 + k*(lambda_501 - lambda_1) / 40;	
		lambdai[k] = inversePowerMethod(mig_u[k]);
	}

	printf("-----------Result-------------\n");
	printf("--question 1\n");
	printf("λ1=%1.11e\n", lambda_1);
	printf("λ501=%1.11e\n", lambda_501);
	printf("λs=%1.11e\n", lambda_s);

	printf("\n--question 2\n");
	for (k = 1; k <= 39; k++)	
		printf("λi%d=%1.11e \n", k, lambdai[k]);

	printf("\n--question 3\n");
	printf("cond(A)=%1.11e\n", fabs(lambda_temp_1 / lambda_s));
	printf("detA=%1.11e \n", det);
	
	system("pause");
	return 0;
}
