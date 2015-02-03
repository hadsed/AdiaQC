#include <stdio.h>
#include <stdlib.h>

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
    &result->size,&result->resident,&result->share,&result->text,&result->lib,&result->data,&result->dt))
  {
    fprintf(stderr, statm_path);
    abort();
  }
  fclose(f);
}

int main(){
    statm_t res;

    read_off_memory_status( &res );
    printf( "Test: %ld\n", res.resident );

    return 0;
}
