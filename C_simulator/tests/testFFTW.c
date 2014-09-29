#include <fftw3.h>
#include <stdio.h>
#include <math.h>

int main(){
    int result[8] = { 4, 2, 0, -2, 0, 2, 0, 2 };
    int i;
    double vec[8] = { 1, 0, 1, 0, 0, 1, 1, 0 };

    fftw_complex in[8];
    fftw_complex out[8];
    fftw_plan p1, p2;

    int dim[3] = { 2, 2, 2 };
    
    printf( "Expected soln:\n" );
    for( i = 0; i < 8; i++ ){
        in[i][0] = vec[i];
        in[i][1] = 0;

        //printf( " %d", result[i] );
        printf( " %d", (int)vec[i] );
    }
    printf( "\n" );

    //p = fftw_plan_dft( 3, dim, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
    p1 = fftw_plan_dft( 3, dim, in, in, FFTW_FORWARD, FFTW_ESTIMATE );
    //p2 = fftw_plan_dft( 3, dim, in, in, FFTW_BACKWARD, FFTW_ESTIMATE );
    fftw_execute(p1);
    fftw_execute(p1);
    fftw_destroy_plan(p1);
    //fftw_destroy_plan(p2);
    
    printf( "Calculated soln:\n" );
    for( i = 0; i < 8; i++ ){
        //printf( " %d", (int)out[i][0] );
        printf( " %d", (int)(in[i][0]/8) );
    }
    printf( "\n" );

    return 0;
}
