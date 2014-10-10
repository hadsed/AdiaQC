#include "parson.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(){
    int i, nQ, lrgs;
    double t, dt;
    double *al, *be, *de;
    JSON_Value *root_value = NULL;
    JSON_Object *root_object;
    //JSON_Object *scalar_object;
    //JSON_Object *coeff_object;
    JSON_Array *array;

    root_value = json_parse_file_with_comments( "config.json" );
    root_object = json_value_get_object( root_value );
    //scalar_object = json_value_get_object( root_object,  );
    //coeff_object = json_value_get_object( root_value );

    nQ = (int)json_object_dotget_number( root_object, "scalars.NQ" );
    lrgs = (int)json_object_dotget_number( root_object, "scalars.LRGS" );
    t = json_object_dotget_number( root_object, "scalars.T" );
    dt = json_object_dotget_number( root_object, "scalars.DT" );

    al = (double *)malloc( nQ*sizeof(double) );
    be = (double *)malloc( ((nQ*(nQ-1))/2)*sizeof(double) );
    de = (double *)malloc( nQ*sizeof(double) );

    array = json_object_dotget_array( root_object, "coefficients.ALPHA" );
    if( array != NULL ){
        for( i = 0; i < json_array_get_count(array); i++ ){
            al[i] = json_array_get_number( array, i );
        }
    }

    array = json_object_dotget_array( root_object, "coefficients.BETA" );
    if( array != NULL ){
        for( i = 0; i < json_array_get_count(array); i++ ){
            be[i] = json_array_get_number( array, i );
        }
    }

    array = json_object_dotget_array( root_object, "coefficients.DELTA" );
    if( array != NULL ){
        for( i = 0; i < json_array_get_count(array); i++ ){
            de[i] = json_array_get_number( array, i );
        }
    }

    json_value_free(root_value);

    printf("nQ = %d\n", nQ);
    printf("lrgs = %d\n", lrgs);
    printf("t = %f\n", t);
    printf("dt = %f\n\n", dt);

    for( i = 0; i < nQ; i++ ){
        printf("al[%d] = %f ", i, al[i]);
    }
    printf("\n\n");

    for( i = 0; i < (nQ*nQ-nQ)/2; i++ ){
        printf("be[%d] = %f ", i, be[i]);
    }
    printf("\n\n");

    for( i = 0; i < nQ; i++ ){
        printf("de[%d] = %f ", i, de[i]);
    }
    printf("\n\n");

    free(al); free(be); free(de);

    return 0;
}
