// constr_unit_circle5.c

#include <stdio.h>
#include <math.h>

#define NMAX 100

#include "optimal_test.h"



double falpha(double alpha)
{
  int i;

 printf("vector_n = %d, alpha = %lf\n", vector_n, alpha);

  for(i=0; i < vector_n; i++)
    gtemp[i] = xnm1[i] - alpha*grad_vector[i]; 

  for(i=0; i < vector_n; i++)
      printf(
    "gtemp[%d] = %lf,  xnm1[%d] = %lf, grad_vector[%d] = %lf\n",
      i, gtemp[i], i, xnm1[i], i,
              grad_vector[i]);

  return  objective_function(gtemp);

} // falpha



double golden(double (*fp)(double), double x1, double x3, double eps)
{
 double x2, fx2,  fx3, x4, fx4;
 double phi = 1.618;
 double phi1 = 2.618;

 x2 = x1 + (x3-x1)/phi1;
 fx2 = (*fp)(x2);
 x4 = x1 + (x3-x1)/phi;
 fx4 = (*fp)(x4);

do {

 if (fx2 > fx4)
 {
   x1 = x2;
   x2 = x4;
   fx2 = fx4;
   x4 = x1 + (x3-x1)/phi;
   fx4 = (*fp)(x4);
 }
 else 
  {
   x3 = x4;
   x4 = x2;
   fx4 = fx2;
   x2 = x1 + (x3-x1)/phi1;
   fx2 = (*fp)(x2);
  } /* else */  
} while ( (x3 - x1) > eps);

 return ( (x1+x3)/2);
} /* golden */




int vector_convergence_test(double arr[], int n, double epsilon)
{
  int i;

  for(i=0; i < n; i++)
    printf("arr[%d] = %lf\n", i, arr[i]);

  for(i=0; i < n; i++)
   if (fabs(arr[i]) > epsilon)
     return 0.0;
  return 1;

} // vector_convergence_test

void copy_vector(double dest[], double source[], int n)
{
  int i;
  
  for(i=0; i < n; i++)
   dest[i] = source [i];  

} // copy_vector


void find_initial_alphas(double (*falpha)(double),
  double *alpha_1, double *alpha_2)
{
  int going_down_flag;
  double falpha1, falpha2, alpha1, alpha2;

  falpha1 = (*falpha)(0.0);

  alpha1 = 0.0;
  alpha2 = 0.0009765625; // 1/1024

  going_down_flag =  1;

  while(going_down_flag == 1)
  {
    falpha2 = (*falpha)(alpha2);
    if(falpha2 >= falpha1)
       going_down_flag = 0;
    else
     { 
       alpha1 = alpha2;
       alpha2 = 2.0*alpha2; 
     } // else   
  } // while

  printf("alpha1 = %lf, alpha2 = %lf\n", alpha1, alpha2);

  *alpha_1 = alpha1;
  *alpha_2 = alpha2;
} // find_initial_alphas


void steepest(double xn[], double x0[], int n, 
FUN_PTR f, GRAD_FUN_PTR grad, double epsilon,
VECTOR_CONVERGENCE_TEST v)
{
    
  double temp, alpha_1, alpha_2, alpha_k;
  int i;
  
  vector_n = n;
  copy_vector(xn, x0, n);
  grad(grad_vector, x0);
  objective_function = f;
  copy_vector(xnm1, xn, n);
  grad(grad_vector, xnm1);

  while(v(grad_vector, n, 0.001) == 0)
  {
      
    find_initial_alphas(falpha,&alpha_1, &alpha_2);
  
    alpha_k = golden(falpha, alpha_1, alpha_2, epsilon);   
    
    printf("alpha_k = %lf\n", alpha_k);
    
    for(i=0; i < n; i++)
      xn[i] = xnm1[i] - alpha_k * grad_vector[i];

    printf("xn:\n");
    for(i=0; i < n; i++)
      printf(" xn[%d] = %lf\n",  i, xn[i]);
    copy_vector(xnm1, xn, n);
    grad(grad_vector, xn);

  } // while 

 
}  // steepest



double original_objective_function(double x[])
{
  double t1, t2, t3, t4;
  
  t1 = x[0] - 1;
  t2 = x[1] - 1;
  t3 = x[2] - 1;
  t4 = x[3] - 1;

  return (t1*t1*t1*t1 + 2*t2*t2*t2*t2+ 3*t3*t3*t3*t3+4*t4*t4*t4*t4);
} // original_objective_function

double h0(double x[])
{
    int i;
    double temp = 0;
    for (i = 0; i < global_n; i++)
        temp = temp + x[i] * x[i];
    temp = temp - 1.0;
    return temp;
} // h0

double f(double x[])
{   
    double term1, term2, term3, term4, temp5;    
    term1 = 4 * pow((x[0] - 1), 3) + 2 * x[0] * x[4];
    term2 = 8 * pow((x[1] - 1), 3) + 2 * x[1] * x[4];
    term3 = 12 * pow((x[2] - 1), 3) + 2 * x[2] * x[4];
    term4 = 16 * pow((x[3] - 1), 3) + 2 * x[3] * x[4];
    return term1 * term1 + term2 * term2 + term3 * term3 + term4 * term4 + h0(x) * h0(x);
} // f




double approx_partial_derivative(double (*obj_f)(double x[]),
    int i, double x[])
{
    double temp1, temp2, xi_orig, result, h;
    double eps_const = 1048576.0;

    xi_orig = x[i];
    h = x[i] / eps_const;

    x[i] = xi_orig + h;

    temp1 = (*obj_f)(x);

    x[i] = xi_orig - h;

    temp2 = (*obj_f)(x);

    result = (temp1 - temp2) / (2 * h);

    x[i] = xi_orig;

    return result;

} // approx_partial_derivative

void g(double grad[], double x[])
{
  int i;
  for(i=0; i < vector_n; i++)
   grad[i] = approx_partial_derivative(f,i, x); 

  for(i=0; i < vector_n; i++)
    printf(" grad[%d] = %lf\n", i, grad[i]);

} // g

double LagrangeFunc(double x[]) {
    return original_objective_function(x) + x[global_n] * h0(x);

}

void compute_delta(double x[])
{
  int i, j;

  for(i=0; i < global_m; i++)
    for(j=0; j < global_n; j++)
       delta_h[i][j] =  approx_partial_derivative(global_h[i],j, x); 
} // compute_delta


double approx_partial_second_order_derivate(int i, int j, double x[]) {
    double xi, dxi, t1, t2;
    double eps_const = 1048576;
    xi = x[i];  dxi = xi / eps_const;
    x[i] = xi + dxi;  t1 = approx_partial_derivative(LagrangeFunc, j, x);
    x[i] = xi - dxi;  t2 = approx_partial_derivative(LagrangeFunc, j, x);
    return (t1 - t2) / (2 * dxi);
}// approx_partial_second_order_derivate

void approx_FLH(double FLH[][NMAX],double x[], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            FLH[i][j] = approx_partial_second_order_derivate(i, j, x);
            if (i != j)  FLH[j][i] = FLH[i][j];
        }
    }
    printf("Hessian F:\n");
    for (int i = 0; i < global_n; i++)    {
        for (int j = 0; j < global_n; j++)
            printf("%10.4lf", FLH[i][j]);
        printf("\n");
    } // for
}


void swap_row(double W[][NMAX], int n, int m1, int m2)
{
    int i;
    double temp;
    for (i = 0; i <= n; i++) {
        temp = W[m1][i];
        W[m1][i] = W[m2][i];
        W[m2][i] = temp;
    } /* for */

} /* swap_rows */

double compute_F_determinant(double F[][NMAX], int n)
{
    int i, j, k, p, itemp, sign;
    double MaxValue, RelativeValue, result;
    double W[NMAX][NMAX];

    sign = 1;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            W[i][j] = F[i][j];

    for (k = 0; k < n; k++)
    {
        p = k;
        MaxValue = fabs(W[k][k]);
        for (i = k + 1; i < n; i++)
            if (fabs(W[i][k]) > MaxValue) {
                p = i;
                MaxValue = fabs(W[i][k]);
            }// if
        if (p != k) {
            swap_row(W, n, k, p);
            sign = 0 - sign;
        } // if
        RelativeValue = W[k][k];

        for (i = k + 1; i < n; i++)
        {
            if (i != k)
            {
                RelativeValue = W[i][k] / W[k][k];;
                W[i][k] = 0.0;
                for (j = k + 1; j <= n; j++)
                    W[i][j] = W[i][j] - RelativeValue * W[k][j];
            } // if
        } // for

    } /* for */

    result = W[0][0];
    for (i = 1; i < n; i++)
        result *= W[i][i];
    result = sign * result;
    return result;
}//compute_F_determinant

void compute_minus_F(double F[][NMAX], int n)
{
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            F[i][j] = -F[i][j];
        } // for

} // compute_minus_F(int n)

int test_non_zero_matrix(double F[][NMAX], int n, double epsilon)
{
    int i, j, flag;
    flag = 0;
    for (i = 0; (i < n) && (flag == 0); i++)
        for (j = 0; (j < n) && (flag == 0); j++)
            if (fabs(F[i][j]) > epsilon)
                flag = 1;

    return flag;

} // test_non_zero_matrix

int positive_definite_test(double F[][NMAX], int n, double epsilon) {
    // result  1 : positive definite    // result  2 : positive semi definite    // result  3 : not positive definite 

    int i, flag, np1;
    double Det_i;
    np1 = n + 1;
    flag = 0;  // Initial
    for (i = 1; (i < np1) && (flag != 3); i++)
    {
        Det_i = compute_F_determinant(F, i);
        printf("Det_i(%d) = %lf\n", i, Det_i);

        if (Det_i < -epsilon)
        {
            flag = 3;  // not positive or semi-positive
        } // if
        else
            if (Det_i > epsilon)
            {
                if (flag == 0)
                    flag = 1;  // possibly positive definite
            } // if
            else // Det_i considered == 0
            {
                if (flag == 1) // up to now 
                         //positive definite
                    flag = 2; // from now on
                              // possibly semi positive definite
            } // else
    } // for

    return flag;

} // positive_definite_test
int definite_test(double F[][NMAX], int n, double epsilon)
{
    // result 0 : zero matrix    // result -1 : negative definite    // result -2 : negative semi definite
    // result  1 : positive definite    // result  2 : positive semi definite    // result  3 : not definite either way

    int flag1, flag2;
    flag1 = test_non_zero_matrix(F, n, epsilon);
    if (flag1 == 0)
        return 0;

    flag1 = positive_definite_test(F, n, epsilon);
    if ((flag1 == 1) || (flag1 == 2))
        return flag1;

    compute_minus_F(F, n);
    flag2 = positive_definite_test(F, n, epsilon);
    compute_minus_F(F, n);

    if ((flag2 == 1) || (flag2 == 2))
        return -flag2;
    else
        return 3;
} // definite_test

static void swap_rows(double W[][NMAX],int m1, int m2){
    int i;
    double temp;
    for (i = 0; i <= 2 * global_n; i++)    {
        temp = W[m1][i];
        W[m1][i] = W[m2][i];
        W[m2][i] = temp;
    } /* for */
} /* swap_rows */

void inv_gen_gaussian() {
    int i, j, k, p, itemp;
    double MaxValue, RelativeValue;
   
    for (i = 0; i < global_m; i++)
        for (j = 0; j < naux; j++)
            W[i][j] = Aaux[i][j];
    for (i = 0; i < global_m; i++)
        for (j = naux; j < 2 * naux; j++)
            W[i][j] = 0.0;

    for (i = 0; i < global_m; i++)
        W[i][naux + i] = 1.0;

    for (k = 0; k < global_m; k++)    {  
        p = k;
        MaxValue = fabs(W[k][k]);
        for (i = k + 1; i < global_m; i++)
            if (fabs(W[i][k]) > MaxValue)     {
                p = i;
                MaxValue = fabs(W[i][k]);
            }// if
        if (p != k)    swap_rows(W, k, p);    

        RelativeValue = W[k][k];
      //  printf("RelativeValue = %6.2lf\n", RelativeValue);

        if (fabs(RelativeValue) < precision_epsilon)        {
            int u, v, itemp;
            double maxvalue = 0.0, temp;
            for (u = k + 1; u < naux; u++)            {
                if (fabs(W[k][u]) > maxvalue)                {
                    v = u;                   
                    maxvalue = fabs(W[k][u]);
                } // if  
            } // for
            itemp = basis[k];
            basis[k] = basis[v];
            basis[v] = itemp;
            for (u = 0; u < global_m; u++)            {
                temp = W[u][k];
                W[u][k] = W[u][v];
                W[u][v] = temp;
            }// for
        } // if

        RelativeValue = W[k][k];
      //  printf("RelativeValue = %6.2lf\n", RelativeValue);
        W[k][k] = 1.0;
        for (j = k + 1; j < 2 * naux; j++)
            W[k][j] = W[k][j] / RelativeValue;

        for (i = 0; i < global_m; i++)        {
            if (i != k)            {
                RelativeValue = W[i][k];
                W[i][k] = 0.0;
                for (j = k + 1; j <= 2 * naux; j++)
                    W[i][j] = W[i][j] - RelativeValue * W[k][j];
            } // if
        } // for      
    } /* for */

    for (i = 0; i < global_m; i++)
        for (j = 0; j < naux; j++)
            B[i][j] = W[i][j + global_m];

} /*  gen_gaussian */


void print_original_system()
{
    int i, j;

    printf("Original System:\n");
    for (i = 0; i < global_m; i++)    {
        for (j = 0; j < global_n; j++)
            printf("%10.3lf", delta_h[i][j]);
        printf("\n");
    } /* for */

    printf("F:\n");
    for (i = 0; i < global_n; i++)    {
        for (j = 0; j < global_n; j++)
            printf("%10.3lf", FLH[i][j]);
        printf("\n");
    } /* for */

} /* print_original_system */

void matrix_multiplication(double Dest[][NMAX],
    double Mat1[][NMAX], double Mat2[][NMAX],
    int n1, int n2, int n3){
    int i, j, k;
    double sum;
    for (i = 0; i < n1; i++)
    {
        for (j = 0; j < n3; j++)
        {
            sum = 0;
            for (k = 0; k < n2; k++)
            {
                sum = sum + Mat1[i][k] * Mat2[k][j];
            } // for
            Dest[i][j] = sum;
        } // for
    } // for1
} //  matrix_multiplication


void compute_E()
{
    int i, j, k, l;
    double sum;

    for (i = 0; i < (global_n - global_m); i++)
    {
        for (k = 0; k < global_n; k++)
            E[i][k] = 0;
       
        E[i][basis[i + global_m]] = 1.0;
        for (k = global_m; k < global_n; k++)
        {
            Eaux[k] = 0;
        } // for
        Eaux[i + global_m] = -1.0;
        sum = 0;
        for (j = 0; j < global_m; j++)
        {
            E[i][basis[j]] = -W[j][i + global_m];
        } // for2
    } // for1
    printf("E:\n");
    for (i = 0; i < global_n - global_m; i++)
    {
        for (j = 0; j < global_n; j++)
            printf("%8.3lf", E[i][j]);
        printf("\n");
    } // for   

    for (i = 0; i < global_n - global_m; i++)
    {
        for (j = 0; j < global_n; j++)
            Et[j][i] = E[i][j];
    } // for   

    printf("Et:\n");
    for (i = 0; i < global_n; i++)
    {
        for (j = 0; j < global_n - global_m; j++)
            printf("%8.3lf", Et[i][j]);
        printf("\n");
    } // for   


} //compute_E

void compute_Reduced_F()
{
    int i, j;
    double sum;  

    print_original_system();
    naux = global_n;
    for (i = 0; i < global_m; i++)
        for (j = 0; j < global_n; j++)
            Aaux[i][j] = delta_h[i][j];

    for (i = 0; i < global_n; i++)
        basis[i] = i;
    basis_size = global_n;

    for (i = 0; i < global_m; i++)
        for (j = 0; j < basis_size; j++)
            Aaux[i][j] = delta_h[i][basis[j]];

    printf("\n\n");
    for (i = 0; i < global_m; i++)
    {
        sum = 0.0;
        for (j = 0; j < global_n; j++)
        {
            sum = sum + y[j] * delta_h[i][j];
          //  printf(" + %6.2lf ", y[j] * delta_h[i][j]);
        } // for
       // printf(" = %6.2lf\n", sum);

    } // for
    inv_gen_gaussian();

    printf("\nBasis:\n");

    for (i = 0; i < basis_size; i++)
        printf(" %d ", basis[i]);
    printf("\n");

    compute_E();
    matrix_multiplication(temp,   E, FLH,
        global_n - global_m, global_n, global_n);

    printf("\ntemp:\n");
    for (i = 0; i < (global_n - global_m); i++)
    {
        for (j = 0; j < global_n; j++)
            printf("   %6.2lf  ", temp[i][j]);
        printf("\n");
    } // for
    printf("\n");

    matrix_multiplication(Reduced_F,
        temp, Et,
        global_n - global_m, global_n, global_n - global_m);

    printf("\nSolution Reduced_F:\n");
    for (i = 0; i < (global_n - global_m); i++)
    {
        for (j = 0; j < (global_n - global_m); j++)
            printf("   %6.3lf  ", Reduced_F[i][j]);
        printf("\n");
    } // for
 
} //compute_Reduced_F()


int main()
{
  double xstar[20], x0[20];
  int i, flag;  

  x0[0] = 0.3;   x0[1] = 0.4;   x0[2] = 0.5;   x0[3] = 0.7;   x0[4] = 1.21;
  global_n = 4;   global_m = 1;
  global_h[0] = h0;

  steepest(xstar, x0, global_n+global_m, 
   f, g, 0.001, vector_convergence_test);

  printf("optimal solution:\n xstar[0] = %lf,  xstar[1] = %lf"
      " xstar[2] = %lf, xstar[3] = %lf\n, xstar[4] = %lf\n",
             xstar[0], xstar[1], xstar[2], xstar[3], xstar[4]);
  printf("optimal value  = %lf\n", f(xstar));

 for(i=0; i < global_m; i++)
    lambda[i] = xstar[i+global_n]; 

  approx_FLH(FLH, xstar, global_n + global_m);
  compute_delta(xstar);
  compute_Reduced_F();

  flag = definite_test(FLH, global_n - global_m, epsilon);
  
  
  if(flag == 1)
   printf("Reduced_R is POSITIVE definite, solution is minimum\n");
  else
  if(flag == 2)
   printf("Reduced_R is POSITIVE SEMI-definite\n");
  else
  if(flag == 3)
   printf("Reduced_R is NON-definite\n");
  else
  if(flag == -1)
   printf("Reduced_R is NEGATIVE-definite, solution is maximum\n");
  else
  if(flag == -2)
   printf("Reduced_R is NEGATIVE SEMI-definite\n");
  else
  if(flag == 0)
   printf("Reduced_R is the zero matrix, both POSITIVE AND NEGATIVE SEMI-definite\n"); 

} // main
