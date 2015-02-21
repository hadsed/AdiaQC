#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include "support/parson.h"
#include <omp.h>
#include <time.h> //for timing purposes


typedef struct {
    unsigned long size,resident,share,text,lib,data,dt;
} statm_t;

void read_off_memory_status(statm_t* result){
    unsigned long dummy;
    const char* statm_path = "/proc/self/statm";
    FILE *f = fopen(statm_path,"r");  

    if(!f){
        fprintf(stderr, statm_path);
        abort();
    }  

    if(7 != fscanf(f,"%ld %ld %ld %ld %ld %ld %ld",
                    &result->size,&result->resident,
                    &result->share,&result->text,&result->lib,
                    &result->data,&result->dt)){
        fprintf(stderr, statm_path);
        abort();
    }

    fclose(f);
}

void expMatTimesVecS( double complex *vec, const double *mat, double complex cc, uint64_t N, double complex scale ){
    uint64_t i;

    #pragma omp parallel for
    for( i = 0UL; i < N; i++ ){
        vec[i] *= scale * cexp( cc*mat[i] );
    }
}

void expMatTimesVec( double complex *vec, const double *mat, double complex cc, uint64_t N ){
    uint64_t i;

    #pragma omp parallel for
    for( i = 0UL; i < N; i++ ){
        vec[i] *= cexp( cc*mat[i] );
    }
}

void FWHT( double complex *vec, uint64_t nQ, const uint64_t dim ){
    uint64_t stride, base, j;

    //Cycle through stages with different butterfly strides
    for( stride = dim / 2; stride >= 1; stride >>= 1 ){   
            //Butterfly index within subvector of (2 * stride) size
            #pragma omp parallel for private(base)
            for( j = 0; j < dim/2; j++ ){   
                base = j - (j & (stride-1));

                uint64_t i0 = base + j +      0;  
                uint64_t i1 = base + j + stride;

                double complex T1 = vec[i0];
                double complex T2 = vec[i1];
                vec[i0] = T1 + T2; 
                vec[i1] = T1 - T2; 
            }
    }   
}

//from: http://www.dreamincode.net/forums/topic/61496-find-n-max-elements-in-unsorted-array/
//author: baavgai
//date: 2014-09-08, 16:30
void addLarger( double value, uint64_t indx, uint64_t *list, double *mag_list, uint64_t size ){
    uint64_t i = 0;
    while( i < size-1 && value > mag_list[i+1] ){
        mag_list[i] = mag_list[i+1];
        list[i] = list[i+1];
        i++;
    }
    list[i] = indx;
    mag_list[i] = value;
}

void findLargest( uint64_t *listN, const double complex *list, uint64_t size, uint64_t N ){
    uint64_t i;
    double temp;
    
    double *mag_list = (double *)calloc( N, sizeof(double) );
    for( i = 0; i < N; i++ ){
        addLarger( cabs( list[i] ), i, listN, mag_list, N );
    }
    for( i = N; i < size; i++ ){
        temp = cabs( list[i] );
        if( temp > *mag_list ){
            addLarger( temp, i, listN, mag_list, N );
        }
    }
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
    double complex *psi;   /* State vector */
    double complex factor;
    double T = 10.0, dt = 0.1;
    uint64_t i, j, k, bcount;
    uint64_t nQ=3UL, N, L=4UL, dim;
    int test, prnt=0;
    uint64_t testi, testj;
    int dzi, dzj;
    //clock_t tend, tbegin = clock();
    struct timeval tend, tbegin;
    double delta;
    
    //tbegin = clock();
    gettimeofday( &tbegin, NULL );

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                         Parse configuration file
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //TODO: going to need logic to handle incomplete config files
    if( argc < 2 ){
        fprintf( stderr, "Need a json configuration file. Terminating...\n" );
        return 1;
    }

    /* Parse file and populate applicable data structures */
    {
        JSON_Value *root_value = NULL;
        JSON_Object *root_object;
        JSON_Array *array;

        root_value = json_parse_file_with_comments( argv[1] );
        root_object = json_value_get_object( root_value );

        nQ = (uint64_t) json_object_dotget_number( root_object, "scalars.nq" );
        prnt = json_object_dotget_boolean( root_object, "scalars.print" );
        L = (uint64_t) json_object_dotget_number( root_object, "scalars.lrgs" );
        T = json_object_dotget_number( root_object, "scalars.t" );
        dt = json_object_dotget_number( root_object, "scalars.dt" );

        al   = (double *)malloc( nQ*sizeof(double) );
        de   = (double *)malloc( nQ*sizeof(double) );
        be   = (double *)malloc( (nQ*(nQ-1)/2)*sizeof(double) );

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

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Compute the Hamiltonian and state vector for the simulation
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
        Create state vector and initialize to 1/sqrt(2^n)*(|00...0> + ... + |11...1>)
    */
    dim = 1 << nQ;
    factor = 1.0/sqrt( dim );
    
    psi  = (double complex *)malloc( (dim)*sizeof(double complex) );
    hz   = (double *)calloc( (dim),sizeof(double) );
    hhxh = (double *)calloc( (dim),sizeof(double) );

    /*
        Assemble Hamiltonian and state vector
    */
    #pragma omp parallel for private(i, j, bcount, testi, testj, dzi, dzj)
    for( k = 0; k < dim; k++ ){
        bcount = 0;
        for( i = 0; i < nQ; i++ ){
            testi = 1 << (nQ - i - 1);
            dzi = ( (k/testi) & 1 ) ? 1 : -1;

            hz[k] += al[i] * dzi;
            hhxh[k] += de[i] * dzi;

            for( j = i+1; j < nQ; j++ ){
                testj = 1 << (nQ - j - 1);
                dzj = ( (k/testj) & 1 ) ? 1 : -1;

                hz[k] += be[bcount] * dzi * dzj;
                bcount++;
            }
        }
            
        psi[k] = factor;
    }
    free( al ); free( be ); free( de );

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                            Run the Simulation
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    double complex cz, cx;
    double t;
    N = (uint64_t)(T / dt);
    for( i = 0; i < N; i++ ){
        t = i*dt;
        //t0 = (i-1)*dt;

        //Time-dependent coefficients
        cz = (-dt * I)*t/(2.0*T);
        //cz = (-dt * I)*t/(T);
        cx = (-dt * I)*(1 - t/T);

        //Evolve system
        expMatTimesVec( psi, hz, cz, dim ); //apply Z part
        FWHT( psi, nQ, dim );
        expMatTimesVec( psi, hhxh, cx, dim ); //apply X part
        FWHT( psi, nQ, dim );
        expMatTimesVecS( psi, hz, cz, dim, 1.0/dim ); //apply Z part
        
        /*
        FWHT( psi, nQ, dim );
        expMatTimesVecS( psi, hhxh, cx, dim, factor ); //apply Z part
        FWHT( psi, nQ, dim );
        expMatTimesVecS( psi, hz, cz, dim, factor ); //apply Z part
        */
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Check solution and clean up
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
    for( i = 0; i < dim; i++ ){
        printf( "|psi[%d]| = %f\n", 
		i, cabs( psi[i]*psi[i] ) );
    }
    */
    uint64_t *largest = (uint64_t *)calloc( L, sizeof(uint64_t) );
    findLargest( largest, psi, dim, L );
    for( i = 0; i < L; ++i ){
        printf( "|psi[%d]| = %.8f\n",
		largest[L-1-i],
		cabs( psi[largest[L-1-i]]*psi[largest[L-1-i]] ) );
    }
    statm_t res;
    read_off_memory_status( &res );
    free( largest );

    /*
        Free work space.
    */
    free( psi );
    free( hz );
    free( hhxh );

    //tend = clock();
    gettimeofday( &tend, NULL );
    delta = ((tend.tv_sec - tbegin.tv_sec)*1000000u + tend.tv_usec - tbegin.tv_usec)/1.e6;
    //printf( "Total time: %f s\n", (double)(tend - tbegin)/CLOCKS_PER_SEC );
    printf( "Total time: %f s\n", delta );
    printf( "Memory used: %ld kB\n", res.resident );

    return 0;
}
