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

void scaleVec( fftw_complex *vec, double s, uint64_t N ){
    uint64_t i;

    for( i = 0UL; i < N; i++ ){
        vec[i] *= s;
    }
}

void expMatTimesVec( fftw_complex *vec, const double *mat, fftw_complex cc, uint64_t N, fftw_complex scale ){
    uint64_t i;

    for( i = 0UL; i < N; i++ ){
        vec[i] *= scale * cexp( cc*mat[i] );
    }
}

/* BEGIN FWHT METHODS */
void FWHT2( fftw_complex *vec, uint64_t nQ ){
    uint64_t i, iter;
    uint64_t temp;
    uint64_t dim = 1UL << nQ;
    fftw_complex veci;
    
    iter = 1UL;
    while( iter <= nQ ){
        temp = 1UL << (nQ-iter);
        
        i = 0UL;
        while( i < dim ){
            veci = vec[i];
            vec[i] += vec[i + temp];
            vec[i + temp] = veci - vec[i + temp];
            
            //test = ((i+1)/temp) & 1; //be careful, uint64_t vs int
            //TODO: issue using minus sign here?
            i += 1UL + ( -((i+1UL)/temp & 1UL) & temp );
        }
        
        ++iter;
    }
}

//void FWHTP( fftw_complex *vec, ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t local_n0, ptrdiff_t local_n1 ){
void FWHTP( fftw_complex *vec, ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t n0, ptrdiff_t n1,
            ptrdiff_t local_n0, ptrdiff_t local_n1, fftw_plan planT, fftw_plan planU ){
    uint64_t j;
    //fftw_plan planT, planU;

    /*
    planT = fftw_mpi_plan_many_transpose( N0, N1, 2, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
        (double *)vec, (double *)vec, MPI_COMM_WORLD, FFTW_ESTIMATE );
    planU = fftw_mpi_plan_many_transpose( N1, N0, 2, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
        (double *)vec, (double *)vec, MPI_COMM_WORLD, FFTW_ESTIMATE );
    */

    for( j = 0UL; j < local_n0; ++j){
        FWHT2( &vec[j*N1], n1 );
    }
    fftw_execute( planT );

    for( j = 0UL; j < local_n1; ++j){
        FWHT2( &vec[j*N0], n0 );
    }
    fftw_execute( planU );

    //fftw_destroy_plan( planT );
    //fftw_destroy_plan( planU );
}
/* END FWHT METHODS */

//from: http://www.dreamincode.net/forums/topic/61496-find-n-max-elements-in-unsorted-array/
//author: baavgai
//date: 2014-09-08, 16:30
void addLarger( double value, uint64_t indx, data_t *list, uint64_t size ){
    uint64_t i = 0UL;
    while( i < size-1UL && value > list[i+1UL].mag ){
        list[i].mag = list[i+1UL].mag;
        list[i].ndx = list[i+1UL].ndx;
        i++;
    }
    list[i].ndx = indx;
    list[i].mag = value;
}

void findLargestP( data_t *listN, const double complex *list, uint64_t size, uint64_t N, uint64_t base ){
    uint64_t i;
    double temp;
    
    for( i = 0UL; i < N; i++ ){
        //addLarger( cabs( list[i] ), base+i, listN, N );
        addLarger( cabs( list[i]*list[i] ), base+i, listN, N );
    }
    for( i = N; i < size; i++ ){
        temp = cabs( list[i]*list[i] );
        //if( temp > *mag_list ){
        if( temp > listN[0UL].mag ){
            addLarger( temp, base+i, listN, N );
        }
    }
}

void findLargestL( data_t *listN, const data_t *list, uint64_t size, uint64_t N ){
    uint64_t i;
    double temp;
    
    for( i = 0UL; i < N; i++ ){
        addLarger( list[i].mag, list[i].ndx, listN, N );
    }
    for( i = N; i < size; i++ ){
        //if( temp > *mag_list ){
        if( list[i].mag > listN[0UL].mag ){
            addLarger( list[i].mag, list[i].ndx, listN, N );
        }
    }
}


/*------------------------------------------------------------------
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
    scalars_t sc = { 3UL, 4UL, 10.0, 0.1 };
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
    fftw_plan planT, planU; //transpose
    double t;
    N = (uint64_t)(sc.T / sc.dt);

    planT = fftw_mpi_plan_many_transpose( N0, N1, 2, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
        (double *)psi, (double *)psi, MPI_COMM_WORLD, FFTW_ESTIMATE ); //transpose
    planU = fftw_mpi_plan_many_transpose( N1, N0, 2, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
        (double *)psi, (double *)psi, MPI_COMM_WORLD, FFTW_ESTIMATE ); //transpose

    for( i = 0; i < N; i++ ){
        t = i*sc.dt;
        //t0 = (i-1)*dt;

        //Time-dependent coefficients
        cz = (-sc.dt * I)*t/(2.0*sc.T);
        cx = (-sc.dt * I)*(1 - t/sc.T);

        //Evolve system
        expMatTimesVec( psi, hz, cz, ldim, 1.0 ); //apply Z part
        //FWHTP( psi, N0, N1, n0, n1, local_N0, local_N1 );
        FWHTP( psi, N0, N1, n0, n1, local_N0, local_N1, planT, planU ); //transpose
        expMatTimesVec( psi, hhxh, cx, ldim, factor ); //apply X part
        //FWHTP( psi, N0, N1, n0, n1, local_N0, local_N1 );
        FWHTP( psi, N0, N1, n0, n1, local_N0, local_N1, planT, planU ); //transpose
        expMatTimesVec( psi, hz, cz, ldim, factor ); //apply Z part
        
        //scaleVec( psi, 1.0/dim, dim );
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Check solution and clean up
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    large = (data_t *)calloc( sc.L, sizeof(data_t) );
    if( id == 0 )
        largest = (data_t *)malloc( np*sc.L*sizeof(data_t) );

    findLargestP( large, psi, ldim, sc.L, local_0_start );

    MPI_Gather( large, sc.L, DATA_T, largest, sc.L, DATA_T, 0, MPI_COMM_WORLD );
    if( id == 0 ){
        for( i = 0; i < sc.L; ++i ){
            large[i].ndx = 0;
            large[i].mag = 0;
        }
        findLargestL( large, largest, np*sc.L, sc.L ); //zero out large?
        for( i = 0; i < sc.L; ++i ){
            printf( "|psi[%llu]| = %f\n", large[i].ndx, large[i].mag );
        }
    }

    /*
        Free work space.
    */
    fftw_destroy_plan( planT ); //transpose
    fftw_destroy_plan( planU ); //transpose
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
