#include <omp.h>
int main(void) { return omp_pause_resource_all(omp_pause_hard); }
