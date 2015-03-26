/*
    Concepts: Simulation of quantum annealing
    Processors: 1
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include "support/parson.h"
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


// from: http://stackoverflow.com/questions/111928/is-there-a-printf-converter-to-print-in-binary-format
// author: EvilTeach, Daniel Lyons
/*
const char *byte_to_binary( int x ){
    static char b[9];
    b[0] = '\0';

    int z;
    for( z = 128; z > 0; z >>= 1 ){
        strcat( b, ( (x & z) == z ) ? "1" : "0" );
    }

    return b;
}
*/

void scaleVec( double complex *vec, double s, uint64_t N ){
    uint64_t i;

    for( i = 0; i < N; i++ ){
        vec[i] *= s;
    }
}

void expMatTimesVec( double complex *vec, const double *mat, double complex cc, uint64_t N ){
    uint64_t i;

    for( i = 0; i < N; i++ ){
        vec[i] *= cexp( cc*mat[i] );
    }
}

/*
void FWHT( double complex *vec, const uint64_t nQ, const uint64_t dim ){
    uint64_t i, iter;
    uint64_t temp;
    int test;
    double complex veci;
    
    iter = 1UL;
    while( iter <= nQ ){
        temp = 1UL << (nQ-iter);
        
        i = 0UL;
        while( i < dim ){
            veci = vec[i];
            vec[i] += vec[i + temp];
            vec[i + temp] = veci - vec[i + temp];
            
            //test = ((i+1)/temp) & 1; //be careful, uint64_t vs int
            //i += 1 + ((-test) & temp);
            i += 1UL + ( -((i+1UL)/temp & 1UL) & temp );
            //if( (i+1)/temp % 2 == 1 )
            //    i += temp + 1;
            //else
            //    ++i;
        }
        
        ++iter;
    }
}
*/

void FWHT( double complex *vec, uint64_t nQ, const uint64_t dim ){
    uint64_t stride, base, j;

    //Cycle through stages with different butterfly strides
    for( stride = dim / 2; stride >= 1; stride >>= 1 ){   
            //Butterfly index within subvector of (2 * stride) size
            for( j = 0; j < dim/2; j++ ){   
                base = j - (j & (stride-1));

                int i0 = base + j +      0;  
                int i1 = base + j + stride;

                double complex T1 = vec[i0];
                double complex T2 = vec[i1];
                vec[i0] = T1 + T2; 
                vec[i1] = T1 - T2; 
            }
    }   
}

//from: http://www.dreamincode.net/forums/topic/61496-find-n-max-elements-in-unsorted-array/
//author: baavgai
//date: 2014-09-08, 1630
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
    double complex *psi;   /* State vector */
    double complex factor;
    double T = 10.0, dt = 0.1;
    uint64_t i, j, k, bcount;
    uint64_t nQ=3, N, L=4, dim;
    int test, prnt=0;
    uint64_t testi, testj;
    int dzi, dzj;
    clock_t tend, tbegin = clock();
    
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
        TODO: keep track of local size and local base
    */
    dim = 1 << nQ;
    factor = 1.0/sqrt( dim );
    
    psi  = (double complex *)malloc( (dim)*sizeof(double complex) );
    hz   = (double *)calloc( (dim),sizeof(double) );
    hhxh = (double *)calloc( (dim),sizeof(double) );
    /*
    for( i = 0; i < dim; ++i ){
        printf( "hz[%llu][%llu] = %f\n", i, i, hz[i] );
    }
    for( i = 0; i < dim; ++i ){
        printf( "hhxh[%llu][%llu] = %f\n", i, i, hhxh[i] );
    }
    */

    /*
        Assemble Hamiltonian and state vector
    */
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
        cx = (-dt * I)*(1 - t/T);

        //Evolve system
        expMatTimesVec( psi, hz, cz, dim ); //apply Z part
        FWHT( psi, nQ, dim );
        expMatTimesVec( psi, hhxh, cx, dim ); //apply X part
        FWHT( psi, nQ, dim );
        expMatTimesVec( psi, hz, cz, dim ); //apply Z part
        
        /* 
            TODO: can probably get some minor speedup by incorporating this 
                  into expMatTimesVec if needed 
        */
        scaleVec( psi, 1.0/dim, dim );
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Check solution and clean up
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
    if( prnt && nQ < 6 ){
        for( i = 0; i < dim; i++ ){
            printf( "|psi[%d]| = %f\n", 
		    i, cabs( psi[i]*psi[i] ) );
        }
    } else {
    */
        uint64_t *largest = (uint64_t *)calloc( L, sizeof(uint64_t) );
        findLargest( largest, psi, dim, L );
        for( i = 0; i < L; ++i ){
            //printf( "psi[%d] = (%f, %f)\t%f\n",
            printf( "|psi[%d]| = %0.16f\n",
		    largest[L-1-i],
		    cabs( psi[largest[L-1-i]]*psi[largest[L-1-i]] ) );
        }
        statm_t res;
        read_off_memory_status( &res );
        free( largest );
    //}

    /*
        Free work space.
    */
    free( psi );
    free( hz );
    free( hhxh );

    tend = clock();
    printf( "Total time: %f s\n", (double)(tend - tbegin)/CLOCKS_PER_SEC );
    printf( "Memory used: %ld kB\n", res.resident );

    return 0;
}
