#include <stdio.h>
#include <stdint.h>
#include <complex.h>
#include <fftw3-mpi.h>


void FWHT( fftw_complex *vec, const uint64_t nQ, const uint64_t dim ){
    uint64_t i, iter;
    uint64_t temp;
    fftw_complex veci;
    
    iter = 1;
    while( iter <= nQ ){
        temp = 1 << (nQ-iter);
        
        i = 0;
        while( i < dim ){
            veci = vec[i];
            vec[i] += vec[i + temp];
            vec[i + temp] = veci - vec[i + temp];
            
            //test = ((i+1)/temp) & 1; //be careful, uint64_t vs int
            i += 1 + ( -((i+1)/temp & 1) & temp );
        }
        
        ++iter;
    }
}


int main(int argc, char **argv){
    ptrdiff_t N0, N1, nQ=4, n0, n1, temp;
    uint64_t dim = 1 << nQ;
    fftw_complex vec[16] = {1,0,1,0,0,1,1,0,0,1,0,1,1,0,0,1};
    fftw_complex res[16];
    fftw_plan planT, planU;
    fftw_complex *data;
    ptrdiff_t alloc_local, local_n0, local_0_start, i, j;
    ptrdiff_t local_n1, local_1_start;
    int id, np;
     
    MPI_Init( &argc, &argv );
    fftw_mpi_init();
     
    MPI_Comm_rank( MPI_COMM_WORLD, &id );
    MPI_Comm_size( MPI_COMM_WORLD, &np );
    n0 = nQ/2;
    n1 = nQ - n0;
    N0 = 1 << n0; //this should be the number of processors?
    //at the very least, np | N1  and  np | N0
    N1 = 1 << n1; //also defines the order of the WHT product

    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_2d_transposed( N0, N1, MPI_COMM_WORLD,
                                                     &local_n0, &local_0_start,
                                                     &local_n1, &local_1_start );
    data = fftw_alloc_complex( alloc_local );

    /*
       NOTE: this should all be doable for odd and even nQ so long as 
             the number of processors divides both strides
    */
    planT = fftw_mpi_plan_many_transpose( N0, N1, 2, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                          (double *)data, (double *)data, MPI_COMM_WORLD, FFTW_ESTIMATE );
    planU = fftw_mpi_plan_many_transpose( N1, N0, 2, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                          (double *)data, (double *)data, MPI_COMM_WORLD, FFTW_ESTIMATE );
     
    /* initialize data to some function my_function(x,y) */
    for( i = 0; i < local_n0; ++i ) for( j = 0; j < N1; ++j )
        data[i*N1 + j] = vec[(local_0_start + i)*N1 + j];

    FWHT( vec, nQ, dim );

    for( j = 0; j < local_n0; ++j){
        FWHT( &data[j*N1], n1, N1 );
    }

    fftw_execute( planT );

    for( j = 0; j < local_n1; ++j){
        FWHT( &data[j*N0], n0, N0 );
    }
    fftw_execute( planU );
    /*
    for( i = 0; i < local_n0; ++i )
        for( j = 0; j < N1; ++j )
            printf( "%d: data[%d] = (%f, %f)\n", id, i*N1+j, creal(data[i*N1+j]), cimag(data[i*N1+j]));
    printf("\n\n");
    */

    //MPI_Gather( data, N1*local_n0, MPI_C_DOUBLE_COMPLEX, res, N1*local_n0, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD );
    //MPI_Gather( data, N1*local_n0, MPI_DOUBLE, res, N1*local_n0, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    if( id == 0 ){
        /*
        for( i = 0; i < dim; ++i ){
            printf( "data[%d] = (%f, %f)\n", i, creal(res[i]), cimag(res[i]) );
            //printf( "data[%d] = %f\n", i, res[i] );
        }
        printf( "\n\n" );
        */
        for( i = 0; i < dim; ++i ){
            printf( "vec[%d] = (%f, %f)\n", i, creal(vec[i]), cimag(vec[i]) );
            //printf( "vec[%d] = %f\n", i, vec[i] );
        }
        printf( "\n\n" );
    }
     
    fftw_destroy_plan( planT );
    fftw_destroy_plan( planU );

    fftw_free( data );
     
    MPI_Finalize();
}
