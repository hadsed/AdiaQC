#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include "parson.h"

typedef struct {
    uint64_t nQ;
    uint64_t L;
    double dt;
    double T;
} scalars_t;

typedef struct {
    uint64_t ndx;
    double mag;
} data_t;

//TODO: trivial to thread these methods 
void scaleVec( fftw_complex *vec, double s, uint64_t N ){
    uint64_t i;

    for( i = 0; i < N; i++ ){
        vec[i] *= s;
    }
}

void expMatTimesVec( fftw_complex *vec, const double *mat, fftw_complex cc, uint64_t N, fftw_complex scale ){
    uint64_t i;

    for( i = 0; i < N; i++ ){
        vec[i] *= scale * cexp( cc*mat[i] );
    }
}

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

//TODO: determine whether creating/destroying each call slows things
void FWHTP( fftw_complex *vec, ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t local_n0, ptrdiff_t local_n1 ){
    uint64_t j;
    fftw_plan planT, planU;

    planT = fftw_mpi_plan_many_transpose( N0, N1, 2, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
        (double *)vec, (double *)vec, MPI_COMM_WORLD, FFTW_ESTIMATE );
    planU = fftw_mpi_plan_many_transpose( N1, N0, 2, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
        (double *)vec, (double *)vec, MPI_COMM_WORLD, FFTW_ESTIMATE );

    for( j = 0; j < local_n0; ++j){
        //FWHT2( &vec[j*N1], n1, N1 );
        FWHT2( &vec[j*N1], n1 );
    }
    fftw_execute( planT );

    for( j = 0; j < local_n1; ++j){
        //FWHT2( &vec[j*N0], n0, N0 );
        FWHT2( &vec[j*N0], n0 );
    }
    fftw_execute( planU );

    fftw_destroy_plan( planT );
    fftw_destroy_plan( planU );
}
/* END FWHT METHODS */

//from: http://www.dreamincode.net/forums/topic/61496-find-n-max-elements-in-unsorted-array/
//author: baavgai
//date: 2014-09-08, 16:30
void addLarger( double value, uint64_t indx, data_t *list, uint64_t size ){
    uint64_t i = 0;
    while( i < size-1 && value > list[i+1].mag ){
        list[i].mag = list[i+1].mag;
        list[i].ndx = list[i+1].ndx;
        i++;
    }
    list[i].ndx = indx;
    list[i].mag = value;
}

void findLargestP( data_t *listN, const double complex *list, uint64_t size, uint64_t N, uint64_t base ){
    uint64_t i;
    double temp;
    
    for( i = 0; i < N; i++ ){
        //addLarger( cabs( list[i] ), base+i, listN, N );
        addLarger( cabs( list[i]*list[i] ), base+i, listN, N );
    }
    for( i = N; i < size; i++ ){
        temp = cabs( list[i]*list[i] );
        //if( temp > *mag_list ){
        if( temp > listN[0].mag ){
            addLarger( temp, base+i, listN, N );
        }
    }
}

void findLargestL( data_t *listN, const data_t *list, uint64_t size, uint64_t N ){
    uint64_t i;
    double temp;
    
    for( i = 0; i < N; i++ ){
        addLarger( list[i].mag, list[i].ndx, listN, N );
    }
    for( i = N; i < size; i++ ){
        //if( temp > *mag_list ){
        if( list[i].mag > listN[0].mag ){
            addLarger( list[i].mag, list[i].ndx, listN, N );
        }
    }
}


/*------------------------------------------------------------------
- Config file format - simplify, don't need xml, but like the structure
{
    "scalars" : {
        "nq" : 3,
        "lrgs" : 4,
        "t" : 10.0,
        "dt" : 0.1
    },
    "coefficients" : {
        "alpha" : [0.112, 0.234, 0.253],
        "beta" : [0.453, 0.533, -0.732, 0.125, -0.653, 0.752],
        "delta" : [1.0, 1.0, 1.0]
    }
}
------------------------------------------------------------------*/
int main( int argc, char **argv ){
    double *hz, *hhxh;     /* hamiltonian components */
    double *al, *be, *de; 
    fftw_complex *psi; //*work;   /* State vector */
    fftw_complex factor;
    scalars_t sc = { 3, 4, 10.0, 0.1 };
    uint64_t dim, i, j, k, bcount, N;
    uint64_t ldim;
    uint64_t testi, testj;
    int dzi, dzj; //TODO: consider using smaller vars for flags and these
    int id, np;
    MPI_Status status;
    MPI_Datatype SCALARS_T, DATA_T, oldtypes[2];
    MPI_Aint offsets[2], extent; //extent2;
    int blockcounts[2];
    ptrdiff_t N0, N1, n0, n1;
    ptrdiff_t alloc_local, local_N0, local_0_start;
    ptrdiff_t local_N1, local_1_start;
    data_t *large, *largest;

    MPI_Init( &argc, &argv );
    fftw_mpi_init();
     
    MPI_Comm_rank( MPI_COMM_WORLD, &id );
    MPI_Comm_size( MPI_COMM_WORLD, &np );
    
    /* Set up MPI data types for defined structs */
    {
        offsets[0] = 0;
        oldtypes[0] = MPI_UNSIGNED_LONG;
        blockcounts[0] = 2;
        MPI_Type_extent( MPI_UNSIGNED_LONG, &extent );
        offsets[1] = 2*extent;
        oldtypes[1] = MPI_DOUBLE;
        blockcounts[1] = 2;
        /*
        MPI_Type_extent( MPI_DOUBLE, &extent2 );
        offsets[2] = 2 * (extent1 + extent2);
        oldtypes[2] = MPI_INT;
        blockcounts[2] = 1;
        */
        MPI_Type_struct( 2, blockcounts, offsets, oldtypes, &SCALARS_T );
        MPI_Type_commit( &SCALARS_T );

        offsets[0] = 0;
        blockcounts[0] = 1;
        offsets[1] = extent;
        blockcounts[1] = 1;
        MPI_Type_struct( 2, blockcounts, offsets, oldtypes, &DATA_T );
        MPI_Type_commit( &DATA_T );
    }

    if( id == 0 ){
        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                             Parse configuration file
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        //TODO: going to need logic to handle incomplete config files
        if( argc < 2 ){
            fprintf( stderr, "Need a json configuration file. Terminating...\n" );
            return 1;
        }

        /* Parse file and populate applicable data structures */
        JSON_Value *root_value = NULL;
        JSON_Object *root_object;
        JSON_Array *array;

        root_value = json_parse_file_with_comments( argv[1] );
        root_object = json_value_get_object( root_value );

        sc.nQ = (uint64_t) json_object_dotget_number( root_object, "scalars.nq" );
        //sc.prnt = json_object_dotget_boolean( root_object, "scalars.print" );
        sc.L = (uint64_t) json_object_dotget_number( root_object, "scalars.lrgs" );
        sc.T = json_object_dotget_number( root_object, "scalars.t" );
        sc.dt = json_object_dotget_number( root_object, "scalars.dt" );

        al   = (double *)malloc( sc.nQ*sizeof(double) );
        de   = (double *)malloc( sc.nQ*sizeof(double) );
        be   = (double *)malloc( (sc.nQ*(sc.nQ-1)/2)*sizeof(double) );

        array = json_object_dotget_array( root_object, "coefficients.alpha" );
        if( array != NULL ){
            for( i = 0; i < json_array_get_count(array); i++ ){
                al[i] = -json_array_get_number( array, i );
            }
        }

        array = json_object_dotget_array( root_object, "coefficients.beta" );
        if( array != NULL ){
            for( i = 0; i < json_array_get_count(array); i++ ){
                be[i] = -json_array_get_number( array, i );
            }
        }

        array = json_object_dotget_array( root_object, "coefficients.delta" );
        if( array != NULL ){
            for( i = 0; i < json_array_get_count(array); i++ ){
                de[i] = -json_array_get_number( array, i );
            }
        }

        json_value_free( root_value );
    }
    MPI_Bcast( &sc, 1, SCALARS_T, 0, MPI_COMM_WORLD );
    if( id != 0 ){
        al   = (double *)malloc( sc.nQ*sizeof(double) );
        de   = (double *)malloc( sc.nQ*sizeof(double) );
        be   = (double *)malloc( (sc.nQ*(sc.nQ-1)/2)*sizeof(double) );
    }
    MPI_Bcast( al, sc.nQ, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( de, sc.nQ, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( be, (sc.nQ*(sc.nQ-1)/2), MPI_DOUBLE, 0, MPI_COMM_WORLD );

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Compute the Hamiltonian and state vector for the simulation
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
        Create state vector and initialize to 1/sqrt(2^n)*(|00...0> + ... + |11...1>)
    */
    dim = 1 << sc.nQ;
    n0 = sc.nQ/2;
    n1 = sc.nQ - n0;
    N0 = 1 << n0; //this should be the number of processors?
    N1 = 1 << n1; //also defines the order of the WHT product
    factor = 1.0/sqrt( dim );
    
    alloc_local = fftw_mpi_local_size_2d_transposed( N0, N1, MPI_COMM_WORLD,
                                                     &local_N0, &local_0_start,
                                                     &local_N1, &local_1_start );
    local_0_start *= N1; 
    psi = fftw_alloc_complex( alloc_local );
    ldim = local_N0*N1;
    hz   = (double *)calloc( ldim, sizeof(double) );
    hhxh = (double *)calloc( ldim, sizeof(double) );

    /*
        Assemble Hamiltonian and state vector
    */
    for( k = 0; k < ldim; k++ ){
        bcount = 0;
        for( i = 0; i < sc.nQ; i++ ){
            testi = 1 << (sc.nQ - i - 1);
            dzi = ( ((k+local_0_start)/testi) & 1 ) ? 1 : -1;
            //test = k/testi & 1; //be careful, uint64_t vs int
            //dzi = (-test & -1) | (-(!test) & 1);

            hz[k] += al[i] * dzi;
            hhxh[k] += de[i] * dzi;

            //for( j = i; j < sc.nQ; j++ ){
            for( j = i+1; j < sc.nQ; j++ ){
                testj = 1 << (sc.nQ - j - 1);
                dzj = ( ((k+local_0_start)/testj) & 1 ) ? 1 : -1;
                //test = k/testj & 1; //be careful, uint64_t vs int
                //dzj = (-test & -1) | (-(!test) & 1);

                hz[k] += be[bcount] * dzi * dzj;
                bcount++;
            }
        }
            
        psi[k] = factor;
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                            Run the Simulation
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    fftw_complex cz, cx;
    double t;
    N = (uint64_t)(sc.T / sc.dt);
    for( i = 0; i < N; i++ ){
        t = i*sc.dt;
        //t0 = (i-1)*dt;

        //Time-dependent coefficients
        cz = (-sc.dt * I)*t/(2.0*sc.T);
        cx = (-sc.dt * I)*(1 - t/sc.T);

        //Evolve system
        expMatTimesVec( psi, hz, cz, ldim, 1.0 ); //apply Z part
        FWHTP( psi, N0, N1, n0, n1, local_N0, local_N1 );
        expMatTimesVec( psi, hhxh, cx, ldim, factor ); //apply X part
        FWHTP( psi, N0, N1, n0, n1, local_N0, local_N1 );
        expMatTimesVec( psi, hz, cz, ldim, factor ); //apply Z part
        
        //scaleVec( psi, 1.0/dim, dim );
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Check solution and clean up
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    large = (data_t *)calloc( sc.L, sizeof(data_t) );
    if( id == 0 )
        largest = (data_t *)calloc( np*sc.L, sizeof(data_t) );

    findLargestP( large, psi, ldim, sc.L, local_0_start );

    for( i = 0; i < sc.L; ++i ){
        printf( "%d: psi[%llu] = %f\n", id, large[i].ndx, large[i].mag );
    }

    /*
    //TODO: get this working
    MPI_Gather( large, sc.L, DATA_T, largest, sc.L, DATA_T, 0, MPI_COMM_WORLD );
    findLargestL( large, largest, np*sc.L, sc.L ); //zero out large?
    if( id == 0 ){
        for( i = 0; i < sc.L; ++i ){
            //printf( "|psi[%d]| = %f\n", large[i].ndx, large[i].mag );
            printf( "|psi[%llu]| = %f\n", large[i].ndx, large[i].mag );
		//cabs( psi[largest[L-1-i]]*psi[largest[L-1-i]] ) );
        }
    }
    */

    /*
        Free work space.
    */
    fftw_free( psi );
    free( hz ); free( hhxh );
    free( al ); free( be ); free( de );
    free( large );
    if( id == 0 )
        free( largest );

    MPI_Type_free( &SCALARS_T );
    MPI_Type_free( &DATA_T );
    MPI_Finalize();

    return 0;
}
