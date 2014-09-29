/*
    Concepts: Simulation of quantum annealing
    Processors: 1
*/

#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
//#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
//#include <mpi.h>
#include "parson.h"
//TODO: include orbtimer?


inline void scaleVec( fftw_complex *vec, const double s, const uint64_t N ){
    uint64_t i;

    for( i = 0; i < N; i++ ){
        vec[i] *= s;
    }
}


inline void expMatTimesVec( fftw_complex *vec, const double *mat, const fftw_complex cc, const uint64_t N ){
    uint64_t i;

    for( i = 0; i < N; i++ ){
        vec[i] *= cexp( cc*mat[i] );
    }
}


//from: http://www.dreamincode.net/forums/topic/61496-find-n-max-elements-in-unsorted-array/
//author: baavgai
//date: 2014-09-08, 1630
/*
void addLarger( double value, double *list, const uint64_t size ){
    uint64_t i = 0;
    while( i < size-1 && value > list[i+1] ){
        list[i++] = list[i+1];
    }
    list[i] = value;
}

//find max by magnitude
void findLargest( double *listN, fftw_complex *list, const uint64_t size, const uint64_t N ){
    uint64_t i;
    double temp;
    
    //listNth = (int*) calloc (nLargest,sizeof(int));
    for( i = 0; i < N; i++ ){
        addLarger( cabs(list[i]), listN, N );
    }
    for( i = N; i < size; i++ ){
        temp = cabs( list[i] );
        if( temp > *listN ){
            addLarger( temp, listN, N );
        }
    }
    //dumpArray(listNth, nLargest);
    //free(listNth);
}
*/

void addLarger( double value, uint64_t indx, uint64_t *list, double *mag_list, const uint64_t size ){
    uint64_t i = 0;
    while( i < size-1 && value > mag_list[i+1] ){
        mag_list[i] = mag_list[i+1];
        list[i] = list[i+1];
        i++;
    }
    list[i] = indx;
    mag_list[i] = value;
}

void findLargest( uint64_t *listN, const fftw_complex *list, const uint64_t size, const uint64_t N ){
    uint64_t i;
    double temp;
    
    double *mag_list = (double *)calloc( N, sizeof(double) );
    for( i = 0; i < N; i++ ){
        //mag_list[i] = cabs( list[i] );
        addLarger( cabs( list[i] ), i, listN, mag_list, N );
    }
    for( i = N; i < size; i++ ){
        temp = cabs( list[i] );
        if( temp > *mag_list ){
            addLarger( temp, i, listN, mag_list, N );
        }
    }
    //dumpArray(listNth, nLargest);
    free( mag_list );
}


/*------------------------------------------------------------------
- Config file format - simplify, don't need xml, but like the structure
{
    "scalars" : {
        "nq" : 3,
        "lrgs" : 4,
        "print" : true
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
    fftw_complex *psi;   /* State vector */
    double T = 10.0, dt = 0.1;
    fftw_complex factor;
    uint64_t i, j, k;
    uint64_t nQ=3, N, L=4;
    uint64_t testi, testj, dim;
    int *fft_dims;
    int dzi, dzj;
    int flag, flag2, prnt=0;
    fftw_plan plan;
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                         Parse configuration file
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //DEBUG
    
    nQ = 1;
    prnt = 1;
    al = (double *)malloc( nQ*sizeof(double) );
    de = (double *)malloc( nQ*sizeof(double) );
    al[0] = 1.0;
    de[0] = 1.0;
    
    //TODO: going to need logic to handle incomplete config files
    /*
    JSON_Value *root_value = NULL;
    JSON_Object *root_object;
    JSON_Array *array;

    if( argc < 2 ){
        fprintf( stderr, "Need a json configuration file. Terminating...\n" );
        return 1;
    }
    root_value = json_parse_file_with_comments( argv[1] );
    root_object = json_value_get_object( root_value );

    nQ = (uint64_t) json_object_dotget_number( root_object, "scalars.nq" );
    prnt = json_object_dotget_boolean( root_object, "scalars.print" );
    L = (uint64_t) json_object_dotget_number( root_object, "scalars.lrgs" );
    T = json_object_dotget_number( root_object, "scalars.t" );
    dt = json_object_dotget_number( root_object, "scalars.dt" );

    al = (double *)malloc( nQ*sizeof(double) );
    de = (double *)malloc( nQ*sizeof(double) );
    if( nQ > 1) be = (double *)malloc( (nQ*(nQ-1)/2)*sizeof(double) );

    array = json_object_dotget_array( root_object, "coefficients.alpha" );
    if( array != NULL ){
        for( i = 0; i < json_array_get_count(array); i++ ){
            al[i] = -json_array_get_number( array, i );
        }
    }

    if( nQ > 1 ){
        array = json_object_dotget_array( root_object, "coefficients.beta" );
        if( array != NULL ){
            for( i = 0; i < json_array_get_count(array); i++ ){
                be[i] = -json_array_get_number( array, i );
            }
        }
    }

    array = json_object_dotget_array( root_object, "coefficients.delta" );
    if( array != NULL ){
        for( i = 0; i < json_array_get_count(array); i++ ){
            de[i] = -json_array_get_number( array, i );
        }
    }

    json_value_free( root_value );
    */

    /*
    //DEBUG
    printf("al = ");
    for( i = 0; i < nQ; i++ ){
        printf("%f ", al[i] );
    }
    printf("\n\n");

    printf("de = ");
    for( i = 0; i < nQ; i++ ){
        printf("%f ", de[i] );
    }
    printf("\n\n");

    printf("be = ");
    for( i = 0; i < (nQ*(nQ-1))/2; i++ ){
        printf("%f ", be[i] );
    }
    printf("\n\n");
    */

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Compute the Hamiltonian and state vector for the simulation
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
        Create state vector and initialize to 1/sqrt(2^n)*(|00...0> + ... + |11...1>)
        TODO: keep track of local size and local base
    */
    dim = 1 << nQ;
    //factor = 1.0/sqrt( dim );
    factor = 1.0;
    
    fft_dims = (int *)malloc( nQ*sizeof(int) );
    psi  = (fftw_complex *)malloc( (dim)*sizeof(fftw_complex) );
    hz   = (double *)malloc( (dim)*sizeof(double) );
    hhxh = (double *)malloc( (dim)*sizeof(double) );

    for( i = 0; i < nQ; i++ ){
        fft_dims[i] = 2;
    }

    plan = fftw_plan_dft( nQ, fft_dims, psi, psi, FFTW_FORWARD, FFTW_MEASURE );

    /*
        Assemble Hamiltonian and state vector
    */
    flag2 = 1;
    uint64_t bcount = 0;
    for( i = 0; i < nQ; i++ ){
        testi = 1 << (nQ - i - 1);

        flag = 1;
        for( j = i; j < nQ; j++ ){
            testj = 1 << (nQ - j - 1);
            dzi = dzj = -1;

            //TODO: when parallel, only go through local indices and such
            for( k = 0; k < dim; k++ ){
                dzi = (k % testi == 0) ? (-1)*dzi : dzi;
                if( flag ){
                    hz[k] += al[i] * dzi;
                    hhxh[k] += de[i] * dzi;
                } else {
                    dzj = (k % testj == 0) ? (-1)*dzj : dzj;
                    //hz[k] += be[i][j] * dzi * dzj;
                    hz[k] += be[bcount] * dzi * dzj;
                    //bcount++;
                }
                if( flag2 ){
                    psi[k] = factor;
                }
            }
            if( !flag ){
                bcount++;
            }
            flag = 0;
            flag2 = 0;
        }
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                            Run the Simulation
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    fftw_complex cz, cx;
    double t;
    N = (uint64_t)(T / dt);
    //plan = fftw_plan_dft( nQ, fft_dims, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE );
    //plan = fftw_plan_dft( nQ, fft_dims, psi, psi, FFTW_FORWARD, FFTW_MEASURE );
    for( i = 0; i < N; i++ ){
        t = i*dt;
        //t0 = (i-1)*dt;

        //Time-dependent coefficients
        cz = (-dt * I)*t/(2.0*T);
        cx = (-dt * I)*(1 - t/T);

        //Evolve system
        expMatTimesVec( psi, hz, cz, dim ); //apply Z part
        fftw_execute( plan );
        //_fwht( psi, 0, dim );
        expMatTimesVec( psi, hhxh, cx, dim ); //apply X part
        fftw_execute( plan );
        //_fwht( psi, 0, dim );
        expMatTimesVec( psi, hz, cz, dim ); //apply Z part
        scaleVec( psi, 1.0/dim, dim );
    }
    fftw_destroy_plan( plan );
    //TODO: need to figure out where to do this, beginning or here
    scaleVec( psi, 1.0/sqrt(dim), dim );

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Check solution and clean up
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //TODO: locally, collect all local largests on one
    //      node, find k largest from that subset
    //TODO: report locations in state vector
    if( prnt && nQ < 6 ){
        for( i = 0; i < dim; i++ ){
            printf( "psi[%d] = (%f, %f)\t%f\n", i, creal(psi[i]), cimag(psi[i]), cabs(psi[i]*psi[i]) );
            //printf( "psi[%d] = (%f, %f)\t%f\n", i, creal(psi[i]), cimag(psi[i]), cabs(psi[i]) );
        }
    } else {
        uint64_t *largest = (uint64_t *)calloc( L, sizeof(uint64_t) );
        findLargest( largest, psi, dim, L );
        for( i = 0; i < L; i++ ){
            printf( "psi[%d] = (%f, %f)\t%f\n", largest[i], creal( psi[ largest[i] ] ), cimag( psi[ largest[i] ] ), cabs( psi[ largest[i] ] ) );
        }
        free( largest );
    }

    /*
        Free work space.
    */
    fftw_free( psi );
    free( fft_dims );
    free( hz );
    free( hhxh );

    return 0;
}
