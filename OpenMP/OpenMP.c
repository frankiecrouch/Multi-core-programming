#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <omp.h> 
#include <math.h>

#define array_size 2000000

//function to perform the final calculation for the pearson correlation
double calculate_pearson(double pearson_array[5]);

int main(int argc, char* argv[]) {

//********DECLARE VARIABLES************//

/*******SERIAL VARIABLES******/
double start_serial, end_serial; // timer variables

//pointers for arrays
double *x;
double *y;
double *xy;
double *xsqr;
double *ysqr;

int i;//iterations

//sum variables
double xsum;
double ysum;
double xysum;
double xsqr_sum;
double ysqr_sum;

xsum = ysum = xysum = xsqr_sum = ysqr_sum = 0; 

double sums_array_serial[5];//array to hold sums for final calculation

double coeff_serial; // pearson coeff from serial code

/*****PARALLEL VARIABLES******/
double start_parallel, end_parallel; // timer variables

//pointers for arrays
double *p_x;
double *p_y;
double *p_xy;
double *p_xsqr;
double *p_ysqr;

double p_xsum;
double p_ysum;
double p_xysum;
double p_xsqr_sum;
double p_ysqr_sum;

p_xsum = p_ysum = p_xysum = p_xsqr_sum = p_ysqr_sum = 0;

double sums_array_parallel[5];//array to hold sums for final calculations

double coeff_parallel; //pearson coeff from parallel code

/********************* SERIAL CODE ************************************/
//start timer
start_serial = omp_get_wtime();
#pragma omp flush(start_serial) // use flush to make sure command is exectued now e.g. make sure we get an accurate time recording

	//allocate memory for the arrays
    x = (double *)malloc(array_size * sizeof(double));
    y = (double *)malloc(array_size * sizeof(double));
    xy = (double *)malloc( array_size * sizeof(double));
    xsqr = (double *)malloc( array_size * sizeof(double));
    ysqr = (double *)malloc( array_size * sizeof(double));

        /* calculate: xsum, ysum, xysum, xsqr_sum, ysqr_sum */
        for (i = 0; i < array_size ; i++) {
                    x[i] = sin(i);
                    y[i] = sin(i+2);

                    xy[i] = x[i] * y[i];
                    xsqr[i] = x[i] * x[i];
                    ysqr[i] = y[i] * y[i];
                    xsum += x[i];
                    ysum += y[i];
                    xysum += xy[i];
                    xsqr_sum += xsqr[i];
                    ysqr_sum += ysqr[i];
        }

    free(x);
    free(y);
    free(xy);
    free(xsqr);
    free(ysqr);

        //assign values to the array ready to do the final calculation
        sums_array_serial[0] = xsum;
        sums_array_serial[1] = ysum;
        sums_array_serial[2] = xysum;
        sums_array_serial[3] = xsqr_sum;
        sums_array_serial[4] = ysqr_sum;

    /* calculate pearson*/
    coeff_serial = calculate_pearson(sums_array_serial);

    
//end timer
end_serial = omp_get_wtime();
#pragma omp flush(end_serial)

    //* pritn schedule type**//
    const char* s = getenv("OMP_SCHEDULE");
    printf("Schedule type and chunk size: %s\n", s);

    //* print the result */
    printf("Serial - Pearson Correlation Coefficient : %f\n", coeff_serial);

    //print run time
    printf("Serial time: %1.2f\n", end_serial-start_serial);


/******************* PARALLEL CODE ************************************/
//start timer
start_parallel = omp_get_wtime(); 
#pragma omp flush(start_parallel)

    //allocate memory for the arrays
    p_x = (double *)malloc(array_size * sizeof(double));
    p_y = (double *)malloc(array_size * sizeof(double));
    p_xy = (double *)malloc(array_size * sizeof(double));
    p_xsqr = (double *)malloc(array_size * sizeof(double));
    p_ysqr = (double *)malloc(array_size * sizeof(double));

    //calculate: xsum, ysum, xysum, xsqr_sum, ysqr_sum in parallel
    #pragma omp parallel for \
        reduction(+:p_xsum, p_ysum, p_xysum, p_xsqr_sum, p_ysqr_sum)\
        schedule(runtime)
            for(i = 0 ; i < array_size ; i++){
                p_x[i] = sin(i);
                p_y[i] = sin(i+2);

                p_xy[i] = p_x[i] * p_y[i];
                p_xsqr[i] = p_x[i] * p_x[i];
                p_ysqr[i] = p_y[i] * p_y[i];

                p_xsum += p_x[i];
                p_ysum += p_y[i];
                p_xysum += p_xy[i];
                p_xsqr_sum += p_xsqr[i];
                p_ysqr_sum += p_ysqr[i];

            }
    
    free(p_x);
    free(p_y);
    free(p_xy);
    free(p_xsqr);
    free(p_ysqr);

        //assign values to the array ready to do the final calculation
        sums_array_parallel[0] += p_xsum;
        sums_array_parallel[1] += p_ysum;
        sums_array_parallel[2] += p_xysum;
        sums_array_parallel[3] += p_xsqr_sum;
        sums_array_parallel[4] += p_ysqr_sum;

        /* calculate pearson */
        coeff_parallel = calculate_pearson(sums_array_parallel);

//end timer
end_parallel = omp_get_wtime();
#pragma omp flush(end_parallel)

    //print result
    printf("Parallel - Pearson Correlation Coefficient: %f\n", coeff_parallel);

    //print run time
    printf("Parallel time: %1.2f\n", end_parallel-start_parallel);

    //print speed up
    printf("Speed up: %1.2f\n \n", (end_serial-start_serial)/(end_parallel-start_parallel));

return 0;
}


//this function takes the results from the calcuate_sums function and calculates the pearson coefficient//
//see report for the formula used
double calculate_pearson(double pearson_array[5]){

    double num; //numerator
    double deno; //denominator

    //calculate the numerator
    num = (pearson_array[2] - (pearson_array[0] * pearson_array[1]/array_size));

    //calculate the denominator
    deno = (pearson_array[3] - (pearson_array[0] * pearson_array[0]/array_size)) * (pearson_array[4] - (pearson_array[1] * pearson_array[1]/array_size));

    //calculate correlation coefficient 
    return num / sqrt(deno);

} 