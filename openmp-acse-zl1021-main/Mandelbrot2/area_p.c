#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

# define NPOINTS 2000
# define MAXITER 2000


struct complex {
    double real;
    double imag;
};


int main() {

    int numoutside = 0;
    double area, error;
    double start, finish;
    struct complex z, c;

    start = omp_get_wtime();
#pragma omp parallel private(z, c) reduction(+:numoutside) default(none)
    {
    #pragma omp for schedule(static,10)
        for (int i = 0; i < NPOINTS; i++) {
            for (int j = 0; j < NPOINTS; j++) {
                c.real = -2.0 + 2.5 * (double) (i) / (double) (NPOINTS) + 1.0e-7;
                c.imag = 1.125 * (double) (j) / (double) (NPOINTS) + 1.0e-7;
                z = c;
                for (int iter = 0; iter < MAXITER; iter++) {
                    double ztemp = (z.real * z.real) - (z.imag * z.imag) + c.real;
                    z.imag = z.real * z.imag * 2 + c.imag;
                    z.real = ztemp;
                    if ((z.real * z.real + z.imag * z.imag) > 4.0e0) {
                        numoutside++;
                        break;
                    }
                }
            }
        }

    }
    finish = omp_get_wtime();

    area = 2.0 * 2.5 * 1.125 * (double) (NPOINTS * NPOINTS - numoutside) / (double) (NPOINTS * NPOINTS);
    error = area / (double) NPOINTS;

    printf("Area of Mandlebrot set = %12.8f +/- %12.8f\n", area, error);
    printf("Time = %12.8f seconds\n", finish - start);
}
