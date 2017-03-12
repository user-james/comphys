#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Define pi as it is used for the test function*/
#define pi 3.14159265358979323846 

FILE *energies;

/* Defining functions to be used*/
double simpsons(double, double, double, double (*f)(double, double), double);
double gaussian(double, double);
double gaussian_sq(double, double);
void test_integration(double);
double H_oscillator(double, double);
double exp_val(double, double);

int main(int argc, char* argv[]){
/*
 * This program calculates the expectation value of a given wave using numerical integration 
 */
    double x0 = 0.5;
    double step = 0.001;
    double expected = 0;

    //test_integration(step);           /* uncomment this line if you want to test the numerical intgration for a given step size */
    
    energies = fopen("energies", "w");
    while(x0 < 1.55){
        expected = exp_val(x0, step);
        fprintf(energies, "%.2f \t %lf\n", x0, expected);
        x0 += 0.05;
    }
    fclose(energies); 


    return 0;
}

double H_oscillator(double x, double x0){
/* Returns the product psi*H(x)*psi at x for a gaussian wave in a harmonic potential 
 *
 * H(x) = [-d/dx^2 + kx^2]/2 with k = 1.3 
 */
    double k = 1.3;
    return ( 1/(x0*x0) - 2*x*x*pow(x0, -4) + 0.5*k*x*x ) * gaussian_sq(x, x0);
}


double exp_val(double x0, double step){
/*
 * Calculates the expected value of the gaussian wave -> exp[-(x/x0)^2] for a given x0
 * and the Hamiltonian H = [-d/dx^2 + kx^2]/2 with k = 1.3
 *
 * This function uses the variance of the Gaussian (x0) to calculate limits that include
 * 5 stdevs either side of the mean. This ensures that whatever x0 we choose the integration
 * is approximately across the entire real line
 */
    double left = -5*x0, right = 5*x0;
    double result = 0;

    result = simpsons(left, right, step, H_oscillator, x0);
    result = result / simpsons(left, right, step, gaussian_sq, x0);


    return result; 

}



double simpsons(double left, double right, double step, double (*f)(double, double), double x0 ){
/*
 * Integrates the function f(x, x0) using the Simpson's Rule for numerical integration between 
 * the limits left and right using a step size = step 
 */
    int N = (right - left) / step;
    int i = 0;
    double result = 0;
    double a = left, b = left;

    for(i=0; i<N; i++){
        a = i*step + left; b = (i+1)*step + left;
        result += (b - a) * ( (*f)(a, x0) + 4*(*f)((a+b)/2, x0) + (*f)(b, x0) ) / 6.0;
    }

    return result;


}

double gaussian_sq(double x, double x0){
/*
 * Returns the squared gaussian function f(x;x0) = exp( -(x/x0)^2 ) with x0 given as a parameter
 */     
    return exp(-2*(x*x)/(x0*x0));
}


double gaussian(double x, double x0){
/*
 * Returns the gaussian function f(x;x0) = exp( -(x/x0)^2 ) with x0 given as a parameter
 */     
    return exp(-(x*x)/(x0*x0));

}

double sinewave(double x, double p){
/*
 * Returns the function f(x) = p*sin(x)
 */
    return p*sin(x);
}


void test_integration(double step){
/*
 * Returns the integrals of several functions and there analytical results
 * This is used to test our numerical intgegration algorithm for a given step size
 * 
 * The input is just the step size and there is no output except results printed to screen
 */
    printf("\nSTEP SIZE = %.2f\n\n\n", step);
    
    double integral = simpsons(0, pi, step, sinewave, 2);
    printf("\tIntegral of 2*Sin(x) on [0, pi]:\n");
    printf("\tAnalytical = 4 \t Numerical = %lf\n\n\n", integral);


    integral = simpsons(-10, 10, step, gaussian, 2);
    printf("\tIntegral of exp(-(x/2)^2) across entire real line:\n");
    printf("\tAnalytical = 3.544908 \t Numerical (on [-10, 10]) = %lf \n\n\n", integral);

    integral = simpsons(-10, 10, step, gaussian, 5); 
    printf("\tIntegral of exp(-(x/5)^2) across entire real line:\n");
    printf("\tAnalytical = 8.862269 \t Numerical (on [-10, 10]) = %lf \n", integral);
    integral = simpsons(-100, 100, step, gaussian, 5);
    printf("\tAnalytical = 8.862269 \t Numerical (on [-100, 100]) = %lf \n\n\n", integral);

    printf("\nNote: Because the gaussian -> exp[-(x/x0)^2] 'spreads out' as you increase x0 if we are integrating across the entire real line\n we must have our integrating range large enough for the given x0\n\n");

}
