// steepest_golden7.c

#include <stdio.h>
#include <math.h>

#define NMAX 100
#define M_PI 3.1415926

typedef double (*FUN_PTR)(double[]); // pointer to function 

typedef void (*GRAD_FUN_PTR)(double grad[], double x[]); // pointer 
//to function 

typedef int VECTOR_CONVERGENCE_TEST(double arr[], int n, double epsilon);

FUN_PTR objective_function;
double grad_vector[NMAX];
double xnm1[NMAX];
double gtemp[NMAX];
double diff[NMAX];
double Hes_matrix[NMAX][NMAX];
int vector_n;

void subtract_vector(double dest[], double source[], int n)
{
    int i;

    for (i = 0; i < n; i++)
        dest[i] = dest[i] - source[i];

} // subtract_vector

double falpha(double alpha)
{
    int i;

    printf("vector_n = %d, alpha = %lf\n", vector_n, alpha);

    for (i = 0; i < vector_n; i++)
        gtemp[i] = xnm1[i] - alpha * grad_vector[i];

    for (i = 0; i < vector_n; i++)
        printf(
            "gtemp[%d] = %lf,  xnm1[%d] = %lf, grad_vector[%d] = %lf\n",
            i, i, i,
            gtemp[i], xnm1[i], grad_vector[i]);

    return  objective_function(gtemp);

} // falpha



double golden(double (*fp)(double), double x1, double x3, double eps)
{
    double x2, fx2, fx3, x4, fx4;
    double phi = 1.618;
    double phi1 = 2.618;

    x2 = x1 + (x3 - x1) / phi1;
    fx2 = (*fp)(x2);
    x4 = x1 + (x3 - x1) / phi;
    fx4 = (*fp)(x4);

    do {

        if (fx2 > fx4)
        {
            x1 = x2;
            x2 = x4;
            fx2 = fx4;
            x4 = x1 + (x3 - x1) / phi;
            fx4 = (*fp)(x4);
        }
        else
        {
            x3 = x4;
            x4 = x2;
            fx4 = fx2;
            x2 = x1 + (x3 - x1) / phi1;
            fx2 = (*fp)(x2);
        } /* else */
    } while ((x3 - x1) > eps);

    return ((x1 + x3) / 2);
} /* golden */




int vector_convergence_test(double arr[], int n, double epsilon)
{
    int i;

    for (i = 0; i < n; i++)
        printf("arr[%d] = %lf\n", i, arr[i]);

    for (i = 0; i < n; i++)
        if (fabs(arr[i]) > epsilon)
            return 0;
    return 1;

} // vector_convergence_test

void copy_vector(double dest[], double source[], int n)
{
    int i;

    for (i = 0; i < n; i++)
        dest[i] = source[i];

} // copy_vector


void find_initial_alphas(double (*falpha)(double),
    double* alpha_2)
{
    int going_down_flag;
    double falpha1, falpha2, alpha1, alpha2;

    falpha1 = (*falpha)(0.0);

    alpha1 = 0.0;
    alpha2 = 0.0009765625; // 1/1024

    going_down_flag = 1;

    while (going_down_flag == 1)
    {
        falpha2 = (*falpha)(alpha2);
        if (falpha2 >= falpha1)
            going_down_flag = 0;
        else
        {
            alpha1 = alpha2;
            falpha1 = falpha2;
            alpha2 = 2.0 * alpha2;
        } // else   
    } // while

    printf("alpha1 = %lf, alpha2 = %lf\n", alpha1, alpha2);

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
    copy_vector(diff, xn, n);

    while ((v(diff, n, 0.001) == 0) || (v(grad_vector, n, 0.001) == 0))
    {

        find_initial_alphas(falpha, &alpha_2);

        alpha_k = golden(falpha, 0.0, alpha_2, epsilon);

        printf("alpha_k = %lf\n", alpha_k);

        for (i = 0; i < n; i++)
            diff[i] = alpha_k * grad_vector[i];
        for (i = 0; i < n; i++)
            xn[i] = xnm1[i] - diff[i];

        printf("xn:\n");
        for (i = 0; i < n; i++)
            printf(" xn[%d] = %lf\n", i, xn[i]);

        copy_vector(xnm1, xn, n);
        grad(grad_vector, xn);

    } // while 
}  // steepest


double f(double x[])
{
    double temp1, temp2;
    //f(x,y) = -sin(x+2y) - cos(3x+4y)
    temp1 = -sin(x[0] + 2 * x[1]);//-sin(x+2y)
    temp2 = -cos(3 * x[0] + 4 * x[1]);//- cos(3x+4y)
    return temp1 + temp2;

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


void approx_g(double grad[], double x[])
{
    int i, j;

    for (i = 0; i < vector_n; i++)
        grad[i] = approx_partial_derivative(f, i, x);

} // approx_g

void swap_rows(double W[][NMAX], int n, int m1, int m2)
{
    int i;
    double temp;
    for (i = 0; i <= n; i++)    {
        temp = W[m1][i];
        W[m1][i] = W[m2][i];
        W[m2][i] = temp;
    } /* for */

} /* swap_rows */

double compute_F_determinant(double F[][NMAX],int n)
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
            if (fabs(W[i][k]) > MaxValue)    {
                p = i;
                MaxValue = fabs(W[i][k]);
            }// if
        if (p != k)        {
            swap_rows(W,n, k, p);
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

void compute_minus_F(double F[][NMAX],int n)
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

int positive_definite_test(double F[][NMAX],int n, double epsilon){
    // result  1 : positive definite    // result  2 : positive semi definite    // result  3 : not positive definite 

    int i, flag, np1;
    double Det_i;
    np1 = n + 1;
    flag = 0;  // Initial
    for (i = 1; (i < np1) && (flag != 3); i++)
    {
        Det_i = compute_F_determinant(F,i);
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
    flag1 = test_non_zero_matrix(F,n, epsilon);
    if (flag1 == 0)
        return 0;

    flag1 = positive_definite_test(F, n, epsilon);
    if ((flag1 == 1) || (flag1 == 2))
        return flag1;

    compute_minus_F(F,n);
    flag2 = positive_definite_test(F, n, epsilon);
    compute_minus_F(F,n);

    if ((flag2 == 1) || (flag2 == 2))
        return -flag2;
    else
        return 3;
} // definite_test

double approx_partial_second_order_derivate(int i, int j, double x[]) {
    double xi,  dxi, t1, t2;
    double eps_const = 1048576;
    xi = x[i];  dxi = xi / eps_const;
    x[i] = xi + dxi;  t1 = approx_partial_derivative(f, j, x);
    x[i] = xi - dxi;  t2 = approx_partial_derivative(f, j, x);   
    return (t1-t2)/(2*dxi);
}// approx_partial_second_order_derivate


void hes_matrix_fun(double x[], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {           
               Hes_matrix[i][j] = approx_partial_second_order_derivate(i, j, x);
               if (i != j)  Hes_matrix[j][i] = Hes_matrix[i][j];
        }
    }
    printf("Hessian F:\n");
    for (int i = 0; i < n; i++)    {
        for (int j = 0; j < n; j++) {
            printf(" %lf ", Hes_matrix[i][j]);
        }
        printf("\n");
    } // for
}

void check(int f) {
    if (f == 1)
        printf("F is POSITIVE definite\n");
    else
        if (f == 2)
            printf("F is POSITIVE SEMI-definite\n");
        else
            if (f == 3)
                printf("F is NON-definite\n");
            else
                if (f == -1)
                    printf("F is NEGATIVE-definite\n");
                else
                    if (f == -2)
                        printf("F is NEGATIVE SEMI-definite\n");
                    else
                        if (f == 0)
                            printf("F is the zero matrix, both POSITIVE AND NEGATIVE SEMI-definite\n");
}
int main()
{
    double xstar[2], x0[2], value;
    int flag;
    // initail points
    x0[0] = 3.14;  x0[1] = -0.76;

    steepest(xstar, x0, 2, f, approx_g, 0.00001, vector_convergence_test);
    printf("\noptimal solution:\nxstar[0] = %lf,  xstar[1] = %lf",
        xstar[0], xstar[1]);
    printf("\nIn degrees: xstar[0] = %lf,  xstar[1] = %lf\n",
        xstar[0] * 180.0 / M_PI, xstar[1] * 180.0 / M_PI);
    printf("\noptimal value  = %lf\n", f(xstar));
    hes_matrix_fun(xstar,2);
    flag = definite_test(Hes_matrix, 2, 0.00001);
    check(flag);
} // main
