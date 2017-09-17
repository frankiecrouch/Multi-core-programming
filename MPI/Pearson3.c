#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <mpi.h> 
#include <math.h>

#define array_size 2000000

//function to calculate the sum of x, y , xy, xsqr, ysqr
void calculate_sums(double *x, double *y, double *xy, double *xsqr, double *ysqr, double sums_array[5], int *iterations, int my_rank);

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
int *array_size_array;

//varibales for MPI_Scatterv
int *sendcounts;    // number of elements to send to each process
int *displs;        // displacement from 0 (of the X or Y array) e.g. where each segment starts from
int rem;            // remainder when array_size is divided by 
int sum;            // used to calculate the displacement

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

if (my_rank == 0){

    printf("Number of processes: %d\n", comm_sz);
    printf("Array size: %d\n", array_size);

    //create buffer of memory for array x and array y
    global_x = (double *)malloc(array_size * sizeof(double));
    global_y = (double *)malloc(array_size * sizeof(double));
        
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
    array_size_array = (int *)malloc(1 * sizeof(int));

    array_size_array[0] = array_size;

    /* calculate: xsum, ysum, xysum, xsqr_sum, ysqr_sum and put the values into an array of size 5*/
    calculate_sums(global_x, global_y, xy, xsqr, ysqr, sums_array_serial, array_size_array, my_rank);

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
/*MPI_Scatterv*/
/*Send each process a chunk of the x and y arrays and store the value in local_x and local_y*/

    //calculate the size of the chunk to send to each process and the displacement from [0] for each array    
    sendcounts = (int *)malloc( sizeof(int)*comm_sz);
    displs = (int *)malloc( sizeof(int)*comm_sz);

    //calculate the remanider when the array_size is divided by the number of processes
    rem = array_size % comm_sz;
    sum = 0;
    // calculate send counts and displacements
    for (int i = 0; i < comm_sz; i++) {
        sendcounts[i] = array_size/comm_sz;
        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }

        displs[i] = sum;
        sum += sendcounts[i];
    }
    
    //create buffer of memory for array x
    local_x = malloc(sendcounts[my_rank] * sizeof(double));
    // perform MPI_scatterv for array x
    MPI_Scatterv(global_x, sendcounts, displs, MPI_DOUBLE, local_x, sendcounts[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);


    //create buffer of memory for array x
    local_y = malloc(sendcounts[my_rank] * sizeof(double));
    // perform MPI_scatterv for array y
    MPI_Scatterv(global_y, sendcounts, displs, MPI_DOUBLE, local_y, sendcounts[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
   




//global_x and global_y can now be freed
if (my_rank ==0) {

    free(global_x);
    free(global_y);
}

/**********************************************************************/
/* all processes calculate: xsum, ysum, xysum, xsqr_sum, ysqr_sum for their chunk of the array*/

    /*first create buffer of memory for array*/
    local_xy = malloc(sendcounts[my_rank] * sizeof(double));
    local_xsqr = malloc(sendcounts[my_rank] * sizeof(double));
    local_ysqr= malloc(sendcounts[my_rank] * sizeof(double));

    /* calculate: xsum, ysum, xysum, xsqr_sum, ysqr_sum and put the values into an array of size 5*/
    calculate_sums(local_x, local_y, local_xy, local_xsqr, local_ysqr, sums_array_parallel, sendcounts, my_rank);

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
void calculate_sums(double *x, double *y, double *xy, double *xsqr, double *ysqr, double sums_array[5], int *iterations, int my_rank){

    int i;
    double xsum;
    double ysum;
    double xysum;
    double xsqr_sum;
    double ysqr_sum;

    xsum = ysum = xysum = xsqr_sum = ysqr_sum = 0;

   
        for (i = 0; i < iterations[my_rank]; i++) {
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