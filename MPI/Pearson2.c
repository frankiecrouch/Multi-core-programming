#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <mpi.h> 
#include <math.h>

#define array_size 2000000

//function to calculate the sum of x, y , xy, xsqr, ysqr
void calculate_sums(double *x, double *y, double *xy, double *xsqr, double *ysqr, double sums_array[5], int iterations);

//function to performe the final calculation for the pearson correlation
double calculate_pearson(double pearson_array[5]);



int main(void) {

//********DECLARE VARIABLES************//

//pointers for arrays
double *global_x;
double *global_y;
double *xy;
double *xsqr;
double *ysqr;
double *local_x; 
double *local_y; 
double *local_xy;
double *local_xsqr;
double *local_ysqr;

int rem;
int array_size_plus;

double coeff_serial, coeff_parallel; //coefficient variable
int comm_sz, my_rank; // MPI variables
int local_array_size; // number of elements each process receives

//arrays to hold sums
double sums_array_serial[5]; 
double sums_array_parallel[5];
double total_sums_array[5];

double start_serial, end_serial, start_parallel, end_parallel; // timer variables

int i;


//initiate MPI
MPI_Init(NULL, NULL); 

MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);//number of processors

MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);//rank of a process

/**********************************************************************/
/* allocate memory for the x and y arrays and populate with values*/
/* this is performed in serial*/ 

    //frist calculate the remainder when the array size is divided by the number of processes
    rem = (array_size % comm_sz);

    array_size_plus = array_size + (comm_sz - rem);

if (my_rank == 0){

    printf("Number of processes: %d\n", comm_sz);
    printf("Array size: %d\n", array_size);
    
    //if the remainder is not 0 => extend the X ad Y array and assign the extra elements with a value of 0
    if (rem != 0){        

        

        global_x = (double *)malloc(array_size_plus * sizeof(double));
        global_y = (double *)malloc(array_size_plus * sizeof(double));

        for (i = 0; i < (comm_sz - rem); i++){
            global_x[(array_size + i)] = 0;
            global_y[(array_size + i)] = 0;
        }
    }

    //if the remainder is 0 => proceed as in Pearson1.c
    else {
        global_x = (double *)malloc(array_size * sizeof(double));
        global_y = (double *)malloc(array_size * sizeof(double));
    }
        
        /* assign values such that: xi = sin(i) and yi = sin(i+2)*/
        for (i = 0; i < array_size; i++) {
            global_x[i] = sin(i);
            global_y[i] = sin (i+2);
        }
}

/********************* SERIAL CODE ************************************/
//start timer
start_serial = MPI_Wtime();

if (my_rank == 0){
    //create buffers for xy, xsqr, ysqr arrays
    xy = (double *)malloc( array_size * sizeof(double));
    xsqr = (double *)malloc( array_size * sizeof(double));
    ysqr = (double *)malloc( array_size * sizeof(double));

    /* calculate: xsum, ysum, xysum, xsqr_sum, ysqr_sum and put the values into an array of size 5*/
    calculate_sums(global_x, global_y, xy, xsqr, ysqr, sums_array_serial, array_size);

    /* calculate pearson*/
    coeff_serial = calculate_pearson(sums_array_serial);

    /* print the result */
    printf("Serial - Pearson Correlation Coefficient : %f\n", coeff_serial);
    
    //end timer
    end_serial = MPI_Wtime();

    //print run time
    printf("Serial time: %1.2f\n", end_serial-start_serial);

    free(xy);
    free(xsqr);
    free(ysqr);

}

/******************* PARALLEL CODE ************************************/

//start timer
start_parallel = MPI_Wtime();

/**********************************************************************/
/*MPI_Scatter*/
/*Send each process a chunk of the x and y arrays and store the value in local_x and local_y*/

    //calculate the size of the chunk to send to each process e.g. the size of the global_x / no of processes
    rem = (array_size % comm_sz);

    if (rem != 0){
        local_array_size = array_size_plus/comm_sz;
    }

    else{
        local_array_size = array_size/comm_sz;
    }
    
	//create buffer of memory for array x
	local_x = (double *)malloc(local_array_size * sizeof(double));
	// perform MPI_scatter for array x
	MPI_Scatter(global_x, local_array_size, MPI_DOUBLE, local_x, local_array_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//create buffer of memory for array y
	local_y = (double *)malloc(local_array_size * sizeof(double));
	// perform MPI_scatter for array y
	MPI_Scatter(global_y, local_array_size, MPI_DOUBLE, local_y, local_array_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//global_x and global_y can now be freed
if (my_rank ==0) {

    free(global_x);
    free(global_y);
}

/**********************************************************************/
/* all processes calculate: xsum, ysum, xysum, xsqr_sum, ysqr_sum for their chunk of the array*/

    /*first create buffer of memory for array*/
    local_xy = (double *)malloc(local_array_size * sizeof(double));
    local_xsqr = (double *)malloc(local_array_size * sizeof(double));
    local_ysqr = (double *)malloc(local_array_size * sizeof(double));

    /* calculate: xsum, ysum, xysum, xsqr_sum, ysqr_sum and put the values into an array of size 5*/
    calculate_sums(local_x, local_y, local_xy, local_xsqr, local_ysqr, sums_array_parallel, local_array_size);

//local arrays can now be freed
free(local_x);
free(local_y);
free(local_xy);
free(local_xsqr);
free(local_ysqr);

/**********************************************************************/
/* MPI_Reduce*/
/* Get the total sum by adding each of the local sums */
	MPI_Reduce(sums_array_parallel, total_sums_array, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


/***********************************************************************/
/*final calculation of the pearson coeffiecient*/
/*this is done in serial by process 0 only*/

if (my_rank == 0){

    /* calculate pearson using the results of the MPI_Reduce*/
    coeff_parallel = calculate_pearson(total_sums_array);

    //print result
    printf("Parallel - Pearson Correlation Coefficient: %f\n", coeff_parallel);

    //end timer
    end_parallel = MPI_Wtime();

    //print run time
    printf("Parallel time: %1.2f\n", end_parallel-start_parallel);

    //print speed up
    printf("Speed up: %1.2f\n", (end_serial-start_serial)/(end_parallel-start_parallel));
}


MPI_Finalize();
return 0;
}


//this function takes the x and y array and calculates sum of x, y , xy, xsqr, ysqr//
//this results are stored in an array (sums_array)//
void calculate_sums(double *x, double *y, double *xy, double *xsqr, double *ysqr, double sums_array[5], int iterations){

    int i;
    double xsum;
    double ysum;
    double xysum;
    double xsqr_sum;
    double ysqr_sum;

    xsum = ysum = xysum = xsqr_sum = ysqr_sum = 0;

        for (i = 0; i < iterations; i++) {
                xy[i] = x[i] * y[i];
                xsqr[i] = x[i] * x[i];
                ysqr[i] = y[i] * y[i];
                xsum += x[i];
                ysum += y[i];
                xysum += xy[i];
                xsqr_sum += xsqr[i];
                ysqr_sum += ysqr[i];
        }

    sums_array[0] = xsum;
    sums_array[1] = ysum;
    sums_array[2] = xysum;
    sums_array[3] = xsqr_sum;
    sums_array[4] = ysqr_sum;

return;
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