#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
#include "cpg_workspace.h"
#include "cpg_solve.h"
#include "func.h"


static c_int print_V = 1;
static c_int i;
static c_int K = 30;
static c_float r_scale = 20.0;
static c_int dim_x = 12;
static c_int dim_u = 4;
static c_float robot_radius = 0.025;
static c_int nb_iterations = 50;
static c_int x_end = 12;		// dim_x;
static c_int A_bar_end = 156;	// dim_x * (1 + dim_x);
static c_int B_bar_end = 204;	// dim_x * (1 + dim_x + dim_u);
static c_int C_bar_end = 252;	// dim_x * (1 + dim_x + dim_u + dim_u);
static c_int z_bar_end = 264;	// dim_x * (1 + dim_x + dim_u + dim_u + 1);
static c_int sigma = 8;		// t_f_guess
static c_float dt = 0.275862069; // 1. / (K - 1) * sigma;
static c_int nb_obstacles = 1;
static c_int nb_obstacles_s = 2;
static c_int nb_obstacles_c = 1;
static c_float p[3] = {0., 0., 0.1}; // p[nb_obstacles*3], 3 because 2-Dimensions and radius
static c_float p_s[8] = {0.25, 0.2, 0.1, 0.15, -0.25, -0.2, 0.15, 0.15}; // p[nb_obstacles_s*4], 4 because 3-Dimensions and radius
static c_float p_c[4] = {-0.375, 0.25, 0.3, 0.1}; // p[nb_obstacles_c*4], 4 because 3-Dimensions and radius
static c_float x_init[12] = {-0.4, -0.4, 0.025, 0., 0., 0., 0., 0., 0., 0., 0., 0.};
static c_float x_final[12] = {0.125, 0.4, 0.175, 0., 0., 0., 0., 0., 0., 0., 0., 0.,};
static c_int x_size = 360; // dim_x * K;
static c_int u_size = 120; // dim_u * K;
static const c_int nt = 150;
static c_float rho_0 = 0.0;
static c_float rho_1 = 0.25;
static c_float rho_2 = 0.9;
static c_float alpha = 2.0;
static c_float beta = 3.2;
static c_float weight_nu = 100000.0;
c_float tr_radius = 5.;

c_float get_linear_cost()
{
	c_float cost = 0.0;
	for (c_int j = 0; j < K; j++)
	{
		cost += CPG_Result.prim->var1[j];
		// cost += CPG_Result.prim->s_prime_0[j];
		cost += CPG_Result.prim->var2[j];
		// cost += CPG_Result.prim->s_prime_1[j];
		cost += CPG_Result.prim->var3[j];
		cost += CPG_Result.prim->var4[j];
	}	return cost;
}

void f_f(c_float x[dim_x], c_float u[dim_u], c_float *result)
{
	f(u[0], u[1], u[2], u[3], x[6], x[9], x[8], x[11], x[7], x[10], x[3], x[4], x[5], result);
}

void A_f(c_float x[dim_x], c_float u[dim_u], c_float *result)
{
	A(u[0], x[6], x[9], x[8], x[11], x[7], x[10], result);
}

void B_f(c_float x[dim_x], c_float u[dim_u], c_float *result)
{
	B(x[6], x[8], x[7], result);
}

void mat_mul(c_int dim_m1, c_int dim_m2, c_float (*m1)[dim_m1], c_float m2[dim_m1 * dim_m2], c_float res[dim_x * dim_m2])
{
	c_float sum;
	for (int i = 0; i < dim_x; ++i)
	{
		for (int j = 0; j < dim_m2; ++j)
		{
			sum = 0.0;
			for (int k = 0; k < dim_m1; ++k)
			{
				sum += m1[i][k] * m2[k * dim_m2 + j];
			}
			res[i * dim_m2 + j] = sum;
		}
	}
}

/* INPUT:  A - array of pointers to rows of a square matrix having dimension dim_x
 *         Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *         The permutation matrix is not stored as a matrix, but in an integer vector P of size dim_x+1 containing column indexes where the permutation matrix has "1". The last element P[dim_x]=S+dim_x, where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S    
*/
c_int LUPDecompose(c_float (*A)[dim_x], c_int *P) {

	c_int i, j, k, imax; 
	c_float maxA, absA;
	c_float tmp[dim_x];
	c_float Tol = 0.000001;

	for (i = 0; i <= dim_x; i++)
	{
		P[i] = i; //Unit permutation matrix, P[dim_x] initialized with dim_x
	}
	for (i = 0; i < dim_x; i++)
	{
		maxA = 0.0;
		imax = i;

		for (k = i; k < dim_x; k++)
		{
			if ((absA = fabs(A[k][i])) > maxA)
			{ 
				maxA = absA;
				imax = k;
			}
		}

		if (maxA < Tol)
		{
			return 0; //failure, matrix is degenerate
		}

		if (imax != i)
		{
			//pivoting P
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			//pivoting rows of A
			for (c_int l = 0; l < dim_x; l++)
			{
				tmp[l] = A[i][l];
			}
			for (c_int l = 0; l < dim_x; l++)
			{
				A[i][l] = A[imax][l];
			}
			for (c_int l = 0; l < dim_x; l++)
			{
				A[imax][l] = tmp[l];
			}

			//counting pivots starting from dim_x (for determinant)
			P[dim_x]++;
		}

		for (j = i + 1; j < dim_x; j++)
		{
			A[j][i] /= A[i][i];
			for (k = i + 1; k < dim_x; k++)
			{
				A[j][k] -= A[j][i] * A[i][k];
			}
		}
	}

	return 1;  //decomposition done 
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; dim_x - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
void LUPSolve(c_float **A, c_int *P, c_float *b, c_float *x) {

	for (c_int i = 0; i < dim_x; i++) {
		x[i] = b[P[i]];

		for (c_int k = 0; k < i; k++)
			x[i] -= A[i][k] * x[k];
	}

	for (c_int i = dim_x - 1; i >= 0; i--) {
		for (c_int k = i + 1; k < dim_x; k++)
			x[i] -= A[i][k] * x[k];

		x[i] /= A[i][i];
	}
}

/* INPUT: A,P filled in LUPDecompose; dim_x - dimension
 * OUTPUT: IA is the inverse of the initial matrix
 */
void LUPInvert(c_float (*A)[dim_x], c_float (*IA)[dim_x]) {
	c_int dec = 0;
	c_int P[dim_x+1];
	dec = LUPDecompose(A, P);
	if(dec == 0){
		printf("ERREUR\n");
		c_int err = 0;
		scanf("%d\n",&err);
		return;
	}
	for (c_int j = 0; j < dim_x; j++)
	{
		for (c_int i = 0; i < dim_x; i++)
		{
			IA[i][j] = P[i] == j ? 1.0 : 0.0;

			for (c_int k = 0; k < i; k++)
			{
				IA[i][j] -= A[i][k] * IA[k][j];
			}
		}

		for (c_int i = dim_x - 1; i >= 0; i--)
		{
			for (c_int k = i + 1; k < dim_x; k++)
			{
				IA[i][j] -= A[i][k] * IA[k][j];
			}

			IA[i][j] /= A[i][i];
		}
	}
	return;
}

/* INPUT: A,P filled in LUPDecompose; dim_x - dimension. 
 * OUTPUT: Function returns the determinant of the initial matrix
 */
c_float LUPDeterminant(c_float (*A)[dim_x]) {

	c_int dec = 0;
	c_int P[dim_x+1];
	dec = LUPDecompose(A, P);
	if(dec == 0){
		printf("ERREUR\n");
		c_int err = 0;
		scanf("%d\n",&err);
		return;
	}

	c_float det = A[0][0];

	for (c_int i = 1; i < dim_x; i++)
		det *= A[i][i];

	return (P[dim_x] - dim_x) % 2 == 0 ? det : -det;
}

c_float determinant(c_float a[dim_x][dim_x], c_float k)
{
	c_int i, j, m, n, c;
	c_float s = 1, det = 0, b[dim_x][dim_x];

	if (k == 1)
	{
		return (a[0][0]);
	}
	else
	{
		det = 0;
		for (c = 0; c < k; c++)
		{
			m = 0;
			n = 0;
			for (i = 0; i < k; i++)
			{
				for (j = 0; j < k; j++)
				{
					b[i][j] = 0;
					if (i != 0 && j != c)
					{
						b[m][n] = a[i][j];
						if (n < (k - 2))
							n++;
						else
						{
							n = 0;
							m++;
						}
					}
				}
			}
			det = det + s * (a[0][c] * determinant(b, k - 1));
			s = -1 * s;
		}
	}

	return det;
}

void transpose(c_float num[dim_x][dim_x], c_float fac[dim_x][dim_x], c_float r, c_float (*res)[dim_x])
{
	c_int i, j;
	c_float b[dim_x][dim_x], d;

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < r; j++)
		{
			b[i][j] = fac[j][i];
		}
	}

	d = determinant(num, r);
	for (i = 0; i < r; i++)
	{
		for (j = 0; j < r; j++)
		{
			res[i][j] = b[i][j] / d;
		}
	}
	return;
}

void mat_inv(c_float num[dim_x][dim_x], c_float (*res)[dim_x])
{
	c_int p, q, m, n, i, j;
	c_float b[dim_x][dim_x], fac[dim_x][dim_x];
	for (q = 0; q < dim_x; q++)
	{
		for (p = 0; p < dim_x; p++)
		{
			m = 0;
			n = 0;
			for (i = 0; i < dim_x; i++)
			{
				for (j = 0; j < dim_x; j++)
				{
					if (i != q && j != p)
					{
						b[m][n] = num[i][j];
						if (n < (dim_x - 2))
							n++;
						else
						{
							n = 0;
							m++;
						}
					}
				}
			}
			fac[q][p] = pow(-1, q + p) * determinant(b, dim_x - 1);
		}
	}
	transpose(num, fac, dim_x, res);
	return;
}

void update_params(c_float weight_nu, c_float tr_radius, c_float *tmp, c_float *A_bar, c_float *B_bar, c_float *C_bar, c_float *z_bar, c_float *X_last, c_float *U_last)
{
	c_int count = 0;

	cpg_update_weight_nu(weight_nu);
	count++;
	for (c_int j = 0; j < (x_size * (nb_obstacles + nb_obstacles_s + nb_obstacles_c)); j++)
	{
		(cpg_update_tmp(j, tmp[j]));
	}
	count += x_size * (nb_obstacles + nb_obstacles_s + nb_obstacles_c);

	for (c_int j = 0; j < (dim_x * dim_x * (K - 1)); j++)
	{
		(cpg_update_A_bar(j, A_bar[j]));
	}
	count += dim_x * dim_x * (K - 1);

	for (c_int j = 0; j < (dim_x * dim_u * (K - 1)); j++)
	{
		(cpg_update_B_bar(j, B_bar[j]));
	}
	count += dim_x * dim_u * (K - 1);

	for (c_int j = 0; j < (dim_x * dim_u * (K - 1)); j++)
	{
		(cpg_update_C_bar(j, C_bar[j]));
	}
	count += dim_x * dim_u * (K - 1);

	for (c_int j = 0; j < (dim_x * (K - 1)); j++)
	{
		(cpg_update_z_bar(j, z_bar[j]));
	}
	count += dim_x * (K - 1);

	for (c_int j = 0; j < x_size; j++)
	{
		(cpg_update_X_last(j, X_last[j]));
	}
	count += x_size;

	for (c_int j = 0; j < u_size; j++)
	{
		(cpg_update_U_last(j, U_last[j]));
	}
	count += u_size;

	cpg_update_tr_radius(tr_radius);
}

c_float norm_1(c_int dim_vector, c_float vector[dim_vector])
{
	c_float sum = 0.0;
	for (c_int j = 0; j < dim_vector; j++)
	{
		if (vector[j] > 0.0)
		{
			sum += vector[j];
		}
		else
		{
			sum -= vector[j];
		}
	}
	return sum;
}

c_float norm_1_3D(c_float x, c_float y, c_float z, c_float a, c_float b, c_float c)
{
	return fabs(x-a) + fabs(y-b) + fabs(z-c);
}

c_float norm_2(c_float x, c_float y, c_float a, c_float b)
{
	c_float x1 = (x - a) * (x - a);
	c_float x2 = (y - b) * (y - b);
	return sqrtf(x1 + x2);
}

c_float norm_2_3D(c_float x, c_float y, c_float z, c_float a, c_float b, c_float c)
{
	c_float x1 = (x - a) * (x - a);
	c_float x2 = (y - b) * (y - b);
	c_float x3 = (z - c) * (z - c);
	return sqrtf(x1 + x2 + x3);
}

void _ode_dVdt(c_float V[dim_x * (1 + dim_x + dim_u + dim_u + 1)], c_float t, c_float *U, c_int k, c_int ind, c_float dVdt[dim_x * (1 + dim_x + dim_u + dim_u + 1)])
{
	c_float x[dim_x];
	c_float u[dim_u];
	c_float V_A_bar[dim_x][dim_x];
	c_float Phi_A_xi[dim_x][dim_x];
	c_float tmp_Phi_A_xi[dim_x][dim_x];
	c_float A_subs[dim_x * dim_x];
	c_float B_subs[dim_x * dim_u];
	c_float f_subs[dim_x];
	c_float z_t[dim_x];
	c_float result_tmp_x[dim_x * dim_x];
	c_float result_tmp_u[dim_x * dim_u];
	c_float result_tmp_f[dim_x];
	c_float result_tmp_f_2[dim_x];

	for (c_int j = 0; j < dim_x; j++)
	{
		x[j] = V[j];
	}

	for (c_int j = 0; j < dim_u; j++)
	{
		u[j] = U[k * dim_u + j] + (t / dt) * (U[(k + 1) * dim_u + j] - U[k * dim_u + j]);
	}

	// printf("matrix A\n");

	for (c_int j = 0; j < dim_x; j++)
	{
		for (c_int l = 0; l < dim_x; l++)
		{
			V_A_bar[j][l] = V[dim_x + j * dim_x + l];
			tmp_Phi_A_xi[j][l] = V[dim_x + j * dim_x + l];
			// printf("%f ",V[dim_x + j * dim_x + l]);
		}
		// printf("\n");
	}
	// mat_inv(V_A_bar, Phi_A_xi);
	LUPInvert(tmp_Phi_A_xi, Phi_A_xi);
	// printf("matrix IA\n");
	// for (c_int j = 0; j < dim_x; j++)
	// {
	// 	for (c_int l = 0; l < dim_x; l++)
	// 	{
	// 		printf("%f ",Phi_A_xi[j][l]);
	// 	}
	// 	printf("\n");
	// }

	// char go ='a';
	// printf("wait : ");
	// scanf("%c", &go);
	// while ((getchar()) != '\n');

	A_f(x, u, A_subs);
	B_f(x, u, B_subs);
	f_f(x, u, f_subs);

	c_int count = 0;
	for (c_int j = 0; j < dim_x; j++)
	{
		dVdt[count + j] = f_subs[j];
	}
	count += dim_x;
	c_float A_subs_matrix[dim_x][dim_x];
	c_float V_A_bar_flatten[dim_x * dim_x];
	for (c_int j = 0; j < dim_x; j++)
	{
		for (c_int l = 0; l < dim_x; l++)
		{
			A_subs_matrix[j][l] = A_subs[j * dim_x + l];
			V_A_bar_flatten[j * dim_x + l] = V_A_bar[j][l];
		}
	}
	mat_mul(dim_x, dim_x, A_subs_matrix, V_A_bar_flatten, result_tmp_x);
	for (c_int j = 0; j < dim_x; j++)
	{
		for (c_int l = 0; l < dim_x; l++)
		{
			dVdt[count + j * dim_x + l] = result_tmp_x[j * dim_x + l];
		}
	}
	count += dim_x * dim_x;

	mat_mul(dim_x, dim_u, Phi_A_xi, B_subs, result_tmp_u); // result * ((dt - t) / dt);
	for (c_int j = 0; j < dim_x; j++)
	{
		for (c_int l = 0; l < dim_u; l++)
		{
			dVdt[count + j * dim_u + l] = result_tmp_u[j * dim_u + l] * ((dt - t) / dt);
		}
	}
	count += dim_x * dim_u;
	for (c_int j = 0; j < dim_x; j++)
	{
		for (c_int l = 0; l < dim_u; l++)
		{
			dVdt[count + j * dim_u + l] = result_tmp_u[j * dim_u + l] * (t / dt);
		}
	}
	count += dim_x * dim_u;
	c_float B_subs_matrix[dim_x][dim_u];
	for (c_int j = 0; j < dim_x; j++)
	{
		for (c_int l = 0; l < dim_u; l++)
		{
			B_subs_matrix[j][l] = B_subs[j * dim_u + l];
		}
	}

	mat_mul(dim_x, 1, A_subs_matrix, x, result_tmp_f);
	mat_mul(dim_u, 1, B_subs_matrix, u, result_tmp_f_2);
	for (c_int j = 0; j < dim_x; j++)
	{
		z_t[j] = f_subs[j] - result_tmp_f[j] - result_tmp_f_2[j];
	}

	mat_mul(dim_x, 1, Phi_A_xi, z_t, result_tmp_f);
	for (c_int j = 0; j < dim_x; j++)
	{
		dVdt[count + j] = result_tmp_f[j];
	}
	count += dim_x * dim_x;

	return;
}

void euler(void (*f)(c_float *, c_float, c_float *, c_int, c_int, c_float *), c_float t0, c_float tf, c_float *y0, c_float *U, c_int k, c_int dim_V, c_float V[dim_V]) // c_float V[dim_x * (1 + dim_x + dim_u + dim_u + 1)])
{
	// c_float y[nt + 1][dim_V]; //[dim_x * (1 + dim_x + dim_u + dim_u + 1)]; // [nt+1][V0.size]
	// c_float y[nt][dim_V];
	c_float y2[2][dim_V];
	c_int flag = 0; // if flag=0 result goes in y2[1] else result goes in y2[0]
	c_float dVdt[dim_V];	  //[dim_x * (1 + dim_x + dim_u + dim_u + 1)];
	c_float h = (tf - t0) / (nt-1);
	c_float t[nt];
	for (c_int j = 0; j < nt - 1; j++)
	{
		t[j] = t0 + j * h; // linspace
	}
	t[nt-1] = tf;
	for (c_int j = 0; j < dim_V; j++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); j++)
	{
		// y[0][j] = y0[j];
		y2[0][j] = y0[j];
	}

	for (c_int j = 1; j < nt; j++)
	{
		// y[j+1] = y[j] + h * f(y[j], t[j], U, k);

		// f(y[j-1], t[j], U, k, j, dVdt);

		if (flag == 0)
		{
			f(y2[0], t[j], U, k, j, dVdt);
			// flag = 1;
		}
		else
		{
			f(y2[1], t[j], U, k, j, dVdt);
			// flag = 0;
		}

		if (j == nt - 1)
		{
			if (flag == 0)
			{
				for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); l++)
				{
					V[l] = y2[0][l] + h * dVdt[l];
				}
			}
			else
			{
				for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); l++)
				{
					V[l] = y2[1][l] + h * dVdt[l];
				}
			}

			// for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); l++)
			// {
			// 	V[l] = y[j-1][l] + h * dVdt[l];
			// }
		}
		else
		{
			if (flag == 0)
			{
				for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); k++)
				{
					y2[1][l] = y2[0][l] + h * dVdt[l];
					flag = 1;
				}
			}
			else
			{
				for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); k++)
				{
					y2[0][l] = y2[1][l] + h * dVdt[l];
					flag = 0;
				}
			}
			
			// for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); k++)
			// {
			// 	y[j][l] = y[j-1][l] + h * dVdt[l];
			// }
		}
	}
	return;
}

void calculate_discretization(c_float *X, c_float *U, c_float *A_bar, c_float *B_bar, c_float *C_bar, c_float *z_bar, c_float *V0)
{
	/* c_float V0[dim_x * (1 + dim_x + dim_u + dim_u + 1)];
	for (c_int j = 0; j < dim_x * (1 + dim_x + dim_u + dim_u + 1); j++)
	{
		V0[j] = 0.0;
	}
	c_int index = dim_x;
	for (c_int j = 0; j < dim_x; j++)
	{
		for (c_int k = 0; k < dim_x; k++)
		{
			if (k == j)
			{
				V0[index] = 1.0;
			}
			else
			{
				V0[index] = 0.0;
			}
			index++;
		}
	} */
	for (c_int k = 0; k < K - 1; k++)
	{
		// printf("V0 : \n");
		for (c_int j = 0; j < dim_x; j++)
		{
			V0[j] = X[k * dim_x + j];
			// printf("%f ",V0[j]);
		}
		// printf("\n");

		c_float V[dim_x * (1 + dim_x + dim_u + dim_u + 1)];

		euler(_ode_dVdt, 0, dt, V0, U, k, (dim_x * (1 + dim_x + dim_u + dim_u + 1)), V);

		c_float Phi[dim_x][dim_x];
		// printf("Phi :\n");
		for (c_int j = 0; j < dim_x; j++)
		{
			for (c_int l = 0; l < dim_x; l++)
			{
				Phi[j][l] = V[dim_x + j * dim_x + l];
				// printf("%f ",Phi[j][l]);
				A_bar[l * dim_x + j + (k * dim_x * dim_x)] = V[dim_x + j * dim_x + l];
			}
			// printf("\n");
		}
		/* for (c_int j = 0; j < dim_x * dim_x; j++)
		{
			A_bar[j + (k * dim_x * dim_x)] = V[j + dim_x];
		} */
		c_float matrix2[dim_x * dim_u];
		// printf("matrix2 :\n");
		c_float result[dim_x * dim_u];
		for (c_int j = 0; j < dim_x; j++)
		{
			for (c_int l = 0; l < dim_u; l++)
			{
				matrix2[j * dim_u + l] = V[A_bar_end + j * dim_u + l];
				// printf("%f ",matrix2[j * dim_u + l]);
			}
			// printf("\n");
		}
		mat_mul(dim_x, dim_u, Phi, matrix2, result);
		// for (c_int j = 0; j < dim_x * dim_u; j++)
		// {
		// 	B_bar[j + k * dim_x * dim_u] = result[j]; 
		// }
		// printf("Result :\n");
		for (c_int j = 0; j < dim_x; j++)
		{
			for (c_int l = 0; l < dim_u; l++)
			{
				B_bar[l * dim_x + j + k * dim_x * dim_u] = result[j*dim_u+l];
				// printf("%f ",result[j*dim_u+l]);
				matrix2[j * dim_u + l] = V[B_bar_end + j * dim_u + l];
			}
			// printf("\n");
		}
		// char go ='a';
		// printf("B_bar done : ");
		// scanf("%c", &go);
		// while ((getchar()) != '\n');
		mat_mul(dim_x, dim_u, Phi, matrix2, result);
		// for (c_int j = 0; j < dim_x * dim_u; j++)
		// {
		// 	C_bar[j + k * dim_x * dim_u] = result[j]; 
		// }
		for (c_int j = 0; j < dim_x; j++)
		{
			for (c_int l = 0; l < dim_u; l++)
			{
				C_bar[l * dim_x + j + k * dim_x * dim_u] = result[j*dim_u+l];
				matrix2[j * dim_u + l] = V[B_bar_end + j * dim_u + l];
			}
		}
		c_float matrix_2[dim_x];
		c_float result2[dim_x];
		for (c_int j = 0; j < dim_x; j++)
		{
			matrix_2[j] = V[C_bar_end + j];
		}
		mat_mul(dim_x, 1, Phi, matrix_2, result2);
		for (c_int j = 0; j < dim_x; j++)
		{
			z_bar[j + k * dim_x] = result2[j];
		}
	}

	return;
}

void _dx(c_float x[dim_x], c_float t, c_float *U, c_int k, c_int ind, c_float *f_subs)
{
	c_float u[dim_u];
	// u = u_t0 + (t / dt) * (u_t1 - u_t0);
	for (c_int j = 0; j < dim_u; j++)
	{
		u[j] = U[k * dim_u + j] + (t / dt) * (U[(k + 1) * dim_u + j] - U[k * dim_u + j]);
	}
	f(u[0], u[1], u[2], u[3], x[6], x[9], x[8], x[11], x[7], x[10], x[3], x[4], x[5], f_subs);

	return;
}

void integrate_nonlinear_piecewise(c_float *new_X, c_float *new_U, c_float *X_nl)
{
	c_float result_euler[dim_x];
	c_float V0[dim_x];
	for (c_int j = 0; j < dim_x; j++)
	{
		X_nl[j] = new_X[j];
	}

	for (c_int k = 0; k < K - 1; k++)
	{
		for (c_int j = 0; j < dim_x; j++)
		{
			V0[j] = new_X[k * dim_x + j];
		}
		euler(_dx, 0, dt, V0, new_U, k, dim_x, result_euler);
		for (c_int j = 0; j < dim_x; j++)
		{
			X_nl[(k + 1) * dim_x + j] = result_euler[j];
		}
	}
	return;
}

c_float get_nonlinear_cost(c_float *X) // TODO modify to take into account every type of obstacles
{
	c_float cost = 0.0;
	c_float dist_to_obstacle;
	for (c_int j = 0; j < nb_obstacles; j++)
	{
		for (c_int k = 0; k < K; k++)
		{
			dist_to_obstacle = norm_2(X[0 + k * dim_x], X[1 + k * dim_x], p[j * 3 + 0], p[j * 3 + 1]);
			if (dist_to_obstacle < (p[j * 3 + 2] + robot_radius))
			{
				cost += p[j * 3 + 2] + robot_radius - dist_to_obstacle;
			}
		}
	}
	for (c_int j = 0; j < nb_obstacles_s; j++)
	{
		for (c_int k = 0; k < K; k++)
		{
			dist_to_obstacle = norm_2_3D(X[0 + k * dim_x], X[1 + k * dim_x], X[2 + k * dim_x], p_s[j * 4 + 0], p_s[j * 4 + 1], p_s[j * 4 + 2]);
			if (dist_to_obstacle < (p_s[j * 4 + 3] + robot_radius))
			{
				cost += p_s[j * 4 + 3] + robot_radius - dist_to_obstacle;
			}
		}
	}
	for (c_int j = 0; j < nb_obstacles_c; j++)
	{
		for (c_int k = 0; k < K; k++)
		{
			dist_to_obstacle = norm_1_3D(X[0 + k * dim_x], X[1 + k * dim_x], X[2 + k * dim_x], p_c[j * 4 + 0], p_c[j * 4 + 1], p_c[j * 4 + 2]);
			if (dist_to_obstacle < (p_c[j * 4 + 3] + robot_radius))
			{
				cost += p_c[j * 4 + 3] + robot_radius - dist_to_obstacle;
			}
		}
	}
	return cost;
}

void dump(c_float *new_X, c_float *new_U, c_float *tmp_obstacle, c_float *A_bar, c_float *B_bar, c_float *C_bar, c_float *z_bar)
{
	FILE *fptrX;
	fptrX = fopen("../../../X_output_c.txt", "w");
	for (i = 0; i < x_size; i++)
	{
		fprintf(fptrX, "%e\n", new_X[i]);
	}
	fclose(fptrX);
	FILE *fptrU;
	fptrU = fopen("../../../U_output_c.txt", "w");
	for (i = 0; i < u_size; i++)
	{
		fprintf(fptrU, "%e\n", new_U[i]);
	}
	fclose(fptrU);
	FILE *fptrA;
	fptrA = fopen("../../../A_output_c.txt", "w");
	for (i = 0; i < (dim_x * dim_x * (K-1)); i++)
	{
		fprintf(fptrA, "%e\n", A_bar[i]);
	}
	fclose(fptrA);
	FILE *fptrB;
	fptrB = fopen("../../../B_output_c.txt", "w");
	for (i = 0; i < (dim_x * dim_u * (K-1)); i++)
	{
		fprintf(fptrB, "%e\n", B_bar[i]);
	}
	fclose(fptrB);
	FILE *fptrC;
	fptrC = fopen("../../../C_output_c.txt", "w");
	for (i = 0; i < (dim_x * dim_u * (K-1)); i++)
	{
		fprintf(fptrC, "%e\n", C_bar[i]);
	}
	fclose(fptrC);
	FILE *fptrz;
	fptrz = fopen("../../../z_output_c.txt", "w");
	for (i = 0; i < (dim_x * (K-1)); i++)
	{
		fprintf(fptrz, "%e\n", z_bar[i]);
	}
	fclose(fptrz);
	FILE *fptrO;
	fptrO = fopen("../../../O_output_c.txt", "w");
	for (i = 0; i < (dim_x * K * (nb_obstacles + nb_obstacles_s + nb_obstacles_c)); i++)
	{
		fprintf(fptrO, "%e\n", tmp_obstacle[i]);
	}
	fclose(fptrO);
	FILE *fptrR;
	fptrR = fopen("../../../R_output_c.txt", "w");
	fprintf(fptrR, "%e\n", tr_radius);
	fclose(fptrR);
	FILE *fptrNu;
	fptrNu = fopen("../../../Nu_output_c.txt", "w");
	fprintf(fptrNu, "%e\n", weight_nu);
	fclose(fptrNu);
}

void read(c_float *new_X, c_float *new_U, c_float *tmp_obstacle, c_float *A_bar, c_float *B_bar, c_float *C_bar, c_float *z_bar, c_float *nu)
{
	c_float num;
	FILE *fptrX;
	if ((fptrX = fopen("../../../X_output.txt","r")) == NULL){
		printf("Error! opening file");
		// Program exits if the file pointer returns NULL.
		exit(1);
	}
	for (i = 0; i < x_size; i++)
	{
		fscanf(fptrX,"%lf\n", &num);
		new_X[i] = num;
	}
	fclose(fptrX);
	FILE *fptrU;
	if ((fptrU = fopen("../../../U_output.txt","r")) == NULL){
		printf("Error! opening file");
		// Program exits if the file pointer returns NULL.
		exit(1);
	}
	for (i = 0; i < u_size; i++)
	{
		fscanf(fptrU,"%lf\n", &num);
		new_U[i] = num;
	}
	fclose(fptrU);
	FILE *fptrR;
	if ((fptrR = fopen("../../../R_output.txt","r")) == NULL){
		printf("Error! opening file");
		// Program exits if the file pointer returns NULL.
		exit(1);
	}
	fscanf(fptrR,"%lf\n", &num);
	tr_radius = num;
	fclose(fptrR);
	FILE *fptrNu;
	if ((fptrNu = fopen("../../../Nu_output.txt","r")) == NULL){
		printf("Error! opening file");
		// Program exits if the file pointer returns NULL.
		exit(1);
	}
	fscanf(fptrNu,"%lf\n", &num);
	weight_nu = num;
	fclose(fptrNu);
	FILE *fptrO;
	if ((fptrO = fopen("../../../O_output.txt","r")) == NULL){
		printf("Error! opening file");
		// Program exits if the file pointer returns NULL.
		exit(1);
	}
	for (i = 0; i < (dim_x * K * (nb_obstacles + nb_obstacles_s + nb_obstacles_c)); i++)
	{
		fscanf(fptrO,"%lf\n", &num);
		tmp_obstacle[i] = num;
	}
	fclose(fptrO);
	FILE *fptrA;
	if ((fptrA = fopen("../../../A_output.txt","r")) == NULL){
		printf("Error! opening file");
		// Program exits if the file pointer returns NULL.
		exit(1);
	}
	for (i = 0; i < (dim_x * dim_x * (K-1)); i++)
	{
		fscanf(fptrA,"%lf\n", &num);
		A_bar[i] = num;
	}
	fclose(fptrA);
	FILE *fptrB;
	if ((fptrB = fopen("../../../B_output.txt","r")) == NULL){
		printf("Error! opening file");
		// Program exits if the file pointer returns NULL.
		exit(1);
	}
	for (i = 0; i < (dim_x * dim_u * (K-1)); i++)
	{
		fscanf(fptrB,"%lf\n", &num);
		B_bar[i] = num;
	}
	fclose(fptrB);
	FILE *fptrC;
	if ((fptrC = fopen("../../../C_output.txt","r")) == NULL){
		printf("Error! opening file");
		// Program exits if the file pointer returns NULL.
		exit(1);
	}
	for (i = 0; i < (dim_x * dim_u * (K-1)); i++)
	{
		fscanf(fptrC,"%lf\n", &num);
		C_bar[i] = num;
	}
	fclose(fptrC);
	FILE *fptrz;
	if ((fptrz = fopen("../../../z_output.txt","r")) == NULL){
		printf("Error! opening file");
		// Program exits if the file pointer returns NULL.
		exit(1);
	}
	for (i = 0; i < (dim_x * (K-1)); i++)
	{
		fscanf(fptrz,"%lf\n", &num);
		z_bar[i] = num;
	}
	fclose(fptrz);
	FILE *fptrnu;
	if ((fptrnu = fopen("../../../nu_output.txt","r")) == NULL){
		printf("Error! opening file");
		// Program exits if the file pointer returns NULL.
		exit(1);
	}
	for (i = 0; i < (dim_x * (K-1)); i++)
	{
		fscanf(fptrnu,"%lf\n", &num);
		nu[i] = num;
	}
	fclose(fptrnu);
}

int main(int argc, char *argv[])
{
	c_float X_nl[x_size];
	c_float X[x_size];
	c_float new_X[x_size];
	c_float U[u_size];
	c_float new_U[u_size];
	c_float tmp_obstacle[dim_x * K * (nb_obstacles + nb_obstacles_s + nb_obstacles_c)];
	c_float A_bar[dim_x * dim_x * (K - 1)];
	c_float B_bar[dim_x * dim_u * (K - 1)];
	c_float C_bar[dim_x * dim_u * (K - 1)];
	c_float z_bar[dim_x * (K - 1)];
	c_float nu[dim_x * (K - 1)];
	c_float V0[dim_x * (1 + dim_x + dim_u + dim_u + 1)];

	c_int first = 1;
	c_int converged = 0;

	c_float last_nonlinear_cost = 0.0;

	c_float actual_change;
	c_float predicted_change;

	for (c_int j = 0; j < dim_x * (1 + dim_x + dim_u + dim_u + 1); j++)
	{
		V0[j] = 0.0;
	}
	c_int index = dim_x;
	for (c_int j = 0; j < dim_x; j++)
	{
		for (c_int k = 0; k < dim_x; k++)
		{
			if (k == j)
			{
				V0[index] = 1.0;
			}
			else
			{
				V0[index] = 0.0;
			}
			index++;
		}
	}

	for (c_int j = 0; j < u_size; j++)
	{
		U[j] = 0.0;
	}
	double alpha1;
	double alpha2;
	for (c_int j = 0; j < K; j++)
	{
		alpha1 = (K - j) * 1.0 / K;
		alpha2 = j * 1.0 / K;
		for (c_int l = 0; l < dim_x; l++)
		{
			X[j * dim_x + l] = x_init[l] * alpha1 + x_final[l] * alpha2;
		}
	}
	c_int it;
	for (it = 0; it < nb_iterations; it++)
	{
		printf("\n\nIteration %d\n", it);
		// printf("X : \n");
		// for (c_int j = 0; j < dim_x; j++)
		// {
		// 	for (c_int l = 0; l < K; l++)
		// 	{
		// 		printf("%f ",X[l * dim_x + j]);
		// 	}
		// 	printf("\n");
		// }
		// printf("U : \n");
		// for (c_int j = 0; j < K; j++)
		// {
		// 	for (c_int l = 0; l < dim_u; l++)
		// 	{
		// 		printf("%f ",U[j * dim_u + l]);
		// 	}
		// 	printf("\n");
		// }

		calculate_discretization(X, U, A_bar, B_bar, C_bar, z_bar, V0);
		c_int nb_total = 0;
		for (c_int j = 0; j < nb_obstacles; j++) // TODO modify to take into account every type of obstacles
		{
			for (c_int k = 0; k < K; k++)
			{
				tmp_obstacle[0 + (nb_total+j) * K * dim_x + k * dim_x] = (X[0 + k * dim_x] - p[j * 3 + 0]) * 1.0 / (0.000001 + norm_2(X[0 + k * dim_x], X[1 + k * dim_x], p[j * 3 + 0], p[j * 3 + 1]));
				
				tmp_obstacle[1 + (nb_total+j) * K * dim_x + k * dim_x] = (X[1 + k * dim_x] - p[j * 3 + 1]) * 1.0 / (0.000001 + norm_2(X[0 + k * dim_x], X[1 + k * dim_x], p[j * 3 + 0], p[j * 3 + 1])); // need to flatten it
			}
		}
		nb_total += nb_obstacles;
		for (c_int j = 0; j < nb_obstacles_s; j++)
		{
			for (c_int k = 0; k < K; k++)
			{
				tmp_obstacle[0 + (nb_total+j) * K * dim_x + k * dim_x] = (X[0 + k * dim_x] - p_s[j * 4 + 0]) * 1.0 / (0.000001 + norm_2_3D(X[0 + k * dim_x], X[1 + k * dim_x], X[2 + k * dim_x], p_s[j * 4 + 0], p_s[j * 4 + 1], p_s[j * 4 + 2]));
				
				tmp_obstacle[1 + (nb_total+j) * K * dim_x + k * dim_x] = (X[1 + k * dim_x] - p_s[j * 4 + 1]) * 1.0 / (0.000001 + norm_2_3D(X[0 + k * dim_x], X[1 + k * dim_x], X[2 + k * dim_x], p_s[j * 4 + 0], p_s[j * 4 + 1], p_s[j * 4 + 2]));

				tmp_obstacle[2 + (nb_total+j) * K * dim_x + k * dim_x] = (X[2 + k * dim_x] - p_s[j * 4 + 2]) * 1.0 / (0.000001 + norm_2_3D(X[0 + k * dim_x], X[1 + k * dim_x], X[2 + k * dim_x], p_s[j * 4 + 0], p_s[j * 4 + 1], p_s[j * 4 + 2]));
			}
		}
		nb_total += nb_obstacles_s;
		for (c_int j = 0; j < nb_obstacles_c; j++)
		{
			for (c_int k = 0; k < K; k++)
			{
				tmp_obstacle[0 + (nb_total+j) * K * dim_x + k * dim_x] = (X[0 + k * dim_x] - p_c[j * 4 + 0]) * 1.0 / (0.000001 + norm_1_3D(X[0 + k * dim_x], X[1 + k * dim_x], X[2 + k * dim_x], p_c[j * 4 + 0], p_c[j * 4 + 1], p_c[j * 4 + 2]));
				
				tmp_obstacle[1 + (nb_total+j) * K * dim_x + k * dim_x] = (X[1 + k * dim_x] - p_c[j * 4 + 1]) * 1.0 / (0.000001 + norm_1_3D(X[0 + k * dim_x], X[1 + k * dim_x], X[2 + k * dim_x], p_c[j * 4 + 0], p_c[j * 4 + 1], p_c[j * 4 + 2]));

				tmp_obstacle[2 + (nb_total+j) * K * dim_x + k * dim_x] = (X[2 + k * dim_x] - p_c[j * 4 + 2]) * 1.0 / (0.000001 + norm_1_3D(X[0 + k * dim_x], X[1 + k * dim_x], X[2 + k * dim_x], p_c[j * 4 + 0], p_c[j * 4 + 1], p_c[j * 4 + 2]));
			}
		}
		nb_total += nb_obstacles_c;

		// dump(X, U, tmp_obstacle, A_bar, B_bar, C_bar, z_bar);
		// printf("dump done\n");
		// char go ='a';
		// printf("wait before read : ");
		// scanf("%c", &go);
		// while ((getchar()) != '\n');
		// read(X, U, tmp_obstacle, A_bar, B_bar, C_bar, z_bar, nu);
		
		// dump(X, U, tmp_obstacle, A_bar, B_bar, C_bar, z_bar);
		// printf("dump done\n");
		// printf("wait before read : ");
		// scanf("%c", &go);
		// while ((getchar()) != '\n');

		update_params(weight_nu, tr_radius, tmp_obstacle, A_bar, B_bar, C_bar, z_bar, X, U);

		while (1)
		{
			printf("Solving...\n");

			cpg_solve();

			for (i = 0; i < x_size; i++)
			{
				new_X[i] = CPG_Result.prim->X[i];
			}
			for (i = 0; i < u_size; i++)
			{
				new_U[i] = CPG_Result.prim->U[i];
			}

			// dump(new_X, new_U, tmp_obstacle, A_bar, B_bar, C_bar, z_bar);
			// printf("dump done\n");
			// char go ='a';
			// printf("wait before read : ");
			// scanf("%c", &go);
			// while ((getchar()) != '\n');
			// read(new_X, new_U, tmp_obstacle, A_bar, B_bar, C_bar, z_bar, nu);

			integrate_nonlinear_piecewise(new_X, new_U, X_nl);

			c_float linear_cost_dynamics = norm_1(dim_x * (K - 1), CPG_Result.prim->nu);
			printf("linear_cost_dynamics = %.9f\n", linear_cost_dynamics);

			c_float dynamics[x_size];
			for (c_int l = 0; l < x_size; l++)
			{
				dynamics[l] = new_X[l] - X_nl[l];
			}
			c_float nonlinear_cost_dynamics = norm_1(x_size, dynamics);
			printf("nonlinear_cost_dynamics : %.9f\n", nonlinear_cost_dynamics);

			c_float linear_cost_constraints = get_linear_cost();
			printf("linear_cost_constraints : %.9f\n", linear_cost_constraints);

			c_float nonlinear_cost_constraints = get_nonlinear_cost(new_X);
			printf("nonlinear_cost_constraints : %.9f\n", nonlinear_cost_constraints);

			c_float linear_cost = linear_cost_dynamics + linear_cost_constraints; // J
			printf("linear_cost : %.9f\n", linear_cost);

			c_float nonlinear_cost = nonlinear_cost_dynamics + nonlinear_cost_constraints; // L
			printf("nonlinear_cost : %.9f\n", nonlinear_cost);

			if (last_nonlinear_cost == 0)
			{
				last_nonlinear_cost = nonlinear_cost;
				for (i = 0; i < x_size; i++)
				{
					X[i] = new_X[i];
				}
				for (i = 0; i < u_size; i++)
				{
					U[i] = new_U[i];
				}
				break;
			}

			actual_change = last_nonlinear_cost - nonlinear_cost; // delta_J
			predicted_change = last_nonlinear_cost - linear_cost; // delta_L

			printf("actual_change : %.9f\n", actual_change);
			printf("predicted_change : %.9f\n", predicted_change);
			if (actual_change < 0.000001 && actual_change > -0.000001)
			{
				converged = 1;
				printf("Converged...\n");
				break;
			}
			else
			{
				printf("Not converged...\n");
				c_float rho = actual_change / predicted_change;
				printf("rho = %.9f\n", rho);
				if (rho < rho_0)
				{
					// reject solution
					tr_radius /= alpha;
					printf("Trust region too large. Solving again with radius = %.9f\n", tr_radius);
				}
				else
				{
					for (i = 0; i < x_size; i++)
					{
						X[i] = new_X[i];
					}
					for (i = 0; i < u_size; i++)
					{
						U[i] = new_U[i];
					}

					printf("Solution accepted.\n");

					if (rho < rho_1)
					{
						printf("Decreasing radius.\n");
						tr_radius /= alpha;
					}
					else if (rho >= rho_2)
					{
						printf("Increasing radius.\n");
						tr_radius *= beta;
					}
					last_nonlinear_cost = nonlinear_cost;
					break;
				}
			}
			cpg_update_tr_radius(tr_radius);
		}

		if (converged)
		{
			printf("Converged after %d iterations.\n", it);
			dump(X, U, tmp_obstacle, A_bar, B_bar, C_bar, z_bar);
			break;
		}
	}
	if (!converged)
	{
		printf("Maximum number of iterations reached without convergence.\n");
	}
	return 0;
}
