#include <stdio.h>
#include <stdint.h>
#include <complex.h>
#include <fftw3-mpi.h>

/* BEGIN FWHT METHODS */
// source: http://www.jjj.de/fxt/fxtbook.pdf
static void sumdiff( fftw_complex *a, fftw_complex *b){
    fftw_complex t = *a - *b;
    *a += *b;
    *b = t;
}

void short_walsh_wak_dif_2( fftw_complex *vec){
    sumdiff( &vec[0], &vec[1] );
}

void short_walsh_wak_dif_4( fftw_complex *vec){
    fftw_complex t0, t1, t2, t3;
    
    t0 = vec[0]; t1 = vec[1];
    t2 = vec[2]; t3 = vec[3];
    sumdiff( &t0, &t2 ); sumdiff( &t1, &t3 );
    sumdiff( &t0, &t1 ); sumdiff( &t2, &t3 );
    vec[0] = t0; vec[1] = t1;
    vec[2] = t2; vec[3] = t3;
}

void short_walsh_wak_dif_4s( fftw_complex *vec, uint64_t s ){
    fftw_complex t0, t1, t2, t3;
    
    {
        uint64_t x = 0;
        t0 = vec[x]; x += s;
        t1 = vec[x]; x += s;
        t2 = vec[x]; x += s;
        t3 = vec[x];
    }
    sumdiff( &t0, &t2 ); sumdiff( &t1, &t3 );
    sumdiff( &t0, &t1 ); sumdiff( &t2, &t3 );
    {
        uint64_t x = 0;
        vec[x] = t0; x += s;
        vec[x] = t1; x += s;
        vec[x] = t2; x += s;
        vec[x] = t3; 
    }
}

void short_walsh_wak_dif_8( fftw_complex *vec ){
    fftw_complex t0, t1, t2, t3, t4, t5, t6, t7;
    
    t0 = vec[0]; t1 = vec[1];
    t2 = vec[2]; t3 = vec[3];
    t4 = vec[4]; t5 = vec[5];
    t6 = vec[6]; t7 = vec[7];
    sumdiff( &t0, &t4 ); sumdiff( &t1, &t5 );
    sumdiff( &t2, &t6 ); sumdiff( &t3, &t7 );
    sumdiff( &t0, &t2 ); sumdiff( &t1, &t3 );
    sumdiff( &t4, &t6 ); sumdiff( &t5, &t7 );
    sumdiff( &t0, &t1 ); sumdiff( &t2, &t3 );
    sumdiff( &t4, &t5 ); sumdiff( &t6, &t7 );
    vec[0] = t0; vec[1] = t1;
    vec[2] = t2; vec[3] = t3;
    vec[4] = t4; vec[5] = t5;
    vec[6] = t6; vec[7] = t7;
}

void short_walsh_wak_dif_16( fftw_complex *vec ){
    fftw_complex t0, t1, t2, t3, t4, t5, t6, t7;
    fftw_complex t8, t9, t10, t11, t12, t13, t14, t15;
    
    t0 = vec[0]; t1 = vec[1]; t2 = vec[2]; t3 = vec[3];
    t4 = vec[4]; t5 = vec[5]; t6 = vec[6]; t7 = vec[7];
    t8 = vec[8]; t9 = vec[9]; t10 = vec[10]; t11 = vec[11];
    t12 = vec[12]; t13 = vec[13]; t14 = vec[14]; t15 = vec[15];
    sumdiff( &t0, &t8 ); sumdiff( &t1, &t9 );
    sumdiff( &t2, &t10 ); sumdiff( &t3, &t11 ); sumdiff( &t4, &t12 );
    sumdiff( &t5, &t13 ); sumdiff( &t6, &t14 ); sumdiff( &t7, &t15 );
    sumdiff( &t0, &t4 ); sumdiff( &t1, &t5 ); sumdiff( &t2, &t6 );
    sumdiff( &t3, &t7 ); sumdiff( &t8, &t12 ); sumdiff( &t9, &t13 );
    sumdiff( &t10, &t14 ); sumdiff( &t11, &t15 ); sumdiff( &t0, &t2 );
    sumdiff( &t1, &t3 ); sumdiff( &t4, &t6 ); sumdiff( &t5, &t7 );
    sumdiff( &t8, &t10 ); sumdiff( &t9, &t11 ); sumdiff( &t12, &t14 );
    sumdiff( &t13, &t15 ); sumdiff( &t0, &t1 ); sumdiff( &t2, &t3 );
    sumdiff( &t4, &t5 ); sumdiff( &t6, &t7 ); sumdiff( &t8, &t9 );
    sumdiff( &t10, &t11 ); sumdiff( &t12, &t13 ); sumdiff( &t14, &t15 );
    vec[0] = t0; vec[1] = t1; vec[2] = t2; vec[3] = t3;
    vec[4] = t4; vec[5] = t5; vec[6] = t6; vec[7] = t7; 
    vec[8] = t8; vec[9] = t9; vec[10] = t10; vec[11] = t11;
    vec[12] = t12; vec[13] = t13; vec[14] = t14; vec[15] = t15;
}

void walsh_wak_dif_4( fftw_complex *vec, uint64_t nQ ){
    const uint64_t n = (1UL << nQ);
    uint64_t i0, mQ, r, j;

    if( n <= 2 ){
        if( n == 2 ) short_walsh_wak_dif_2( vec );
        return;
    }

    for( mQ = nQ; mQ > 3; mQ -= 2){
        uint64_t m = (1UL << mQ);
        uint64_t m4 = (m >> 2);
        for( r = 0; r < n; r += m){
            for( j = 0; j < m4; j++)
                short_walsh_wak_dif_4s( vec + j + r, m4 );
        }
    }

    if( nQ & 1 ){ //n is not a power of 4, need radix-8 step
        for( i0 = 0; i0 < n; i0 += 8 )
            short_walsh_wak_dif_8( vec + i0 );
    } else {
        for( i0 = 0; i0 < n; i0 += 4 )
            short_walsh_wak_dif_4( vec + i0 );
    }
}

void FWHT2( fftw_complex *vec, uint64_t nQ ){
    uint64_t mQ, t1, t2;

    if( nQ <= 13){ //want (2**nQ)*sizeof(double complex) <= L1-cache 
        walsh_wak_dif_4( vec, nQ );
        return;
    }

    for( mQ = 1; mQ <= nQ; ++mQ){
        const uint64_t m = (1UL << mQ);
        const uint64_t mh = (m >> 1);
        for( t1=0, t2=mh; t1 < mh; ++t1, ++t2 ) sumdiff( &vec[t1], &vec[t2] );
    }

    //Recursion
    short_walsh_wak_dif_2( vec + 2 ); //mQ == 1
    short_walsh_wak_dif_4( vec + 4 ); //mQ == 2
    short_walsh_wak_dif_8( vec + 8 ); //mQ == 3
    short_walsh_wak_dif_16( vec + 16 ); //mQ == 4
    for( mQ = 5; mQ < nQ; ++mQ ) FWHT2( vec + (1UL << mQ), mQ );
}

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
/* END FWHT METHODS */


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
