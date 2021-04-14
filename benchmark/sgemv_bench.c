/*
 *https://forums.developer.nvidia.com/t/cublas-vs-cblas-sgemv-benchmarking-matrix-vector-operations-on-gpu-and-cpu/14878
 */
#include <bench.h>

#include <cblas.h>

int main(int argc, char** argv)
{

double time1, timeg;

int nbIter = 10000;

int m;

int n = 128;

for (int j = 0; j < 16; ++j) {

m = 16 << j;

// n = m;

printf("-------------\nEvaluating %i iterations for a matrix %ix%i\n", nbIter, m, n);

float *mat, *x, *y;

float *data = (float*) malloc(sizeof(float) * m * n);

for (int i = 0; i < m*n; ++i)

  data[i] = ((float)i) / ((float)(m * n));


mat = (float*) malloc(m * n * sizeof(float));

x = (float*) malloc(n*sizeof(float));

y = (float*) malloc(m*sizeof(float));

memcpy(mat, data, m * n * sizeof(float));

memcpy(x, data, n * sizeof(float));

memcpy(y, data, m * sizeof(float));

timeg = 0.;
  
for (int i = 0; i < nbIter; ++i)
{
begin();

  cblas_sgemv(CblasColMajor, CblasTrans, n, m, 1, mat, n, x, 1, 1, y, 1);

end();
timeg += getsec();

}
printf("CPU Time: %10.8f (secs)\n", timeg/(double)nbIter );

free(mat);

free(x);

free(y);

free(data);

}

m = 128;

for (int j = 0; j < 16; ++j) {

n = 16 << j;

// n = m;

printf("-------------\nEvaluating %i iterations for a matrix %ix%i\n", nbIter, m, n);

float *mat, *x, *y;

float *data = (float*) malloc(sizeof(float) * m * n);

for (int i = 0; i < m*n; ++i)

  data[i] = ((float)i) / ((float)(m * n));


mat = (float*) malloc(m * n * sizeof(float));

x = (float*) malloc(n*sizeof(float));

y = (float*) malloc(m*sizeof(float));

memcpy(mat, data, m * n * sizeof(float));

memcpy(x, data, n * sizeof(float));

memcpy(y, data, m * sizeof(float));

timeg = 0.;
  
for (int i = 0; i < nbIter; ++i)
{
begin();

  cblas_sgemv(CblasColMajor, CblasTrans, n, m, 1, mat, n, x, 1, 1, y, 1);

end();
timeg += getsec();

}
printf("CPU Time: %10.8f (secs)\n", timeg/(double)nbIter );

free(mat);

free(x);

free(y);

free(data);

}


for (int j = 0; j < 12; ++j) {

m = 16 << j;

n = m;

printf("-------------\nEvaluating %i iterations for a matrix %ix%i\n", nbIter, m, n);

float *mat, *x, *y;

float *data = (float*) malloc(sizeof(float) * m * n);

for (int i = 0; i < m*n; ++i)

  data[i] = ((float)i) / ((float)(m * n));


mat = (float*) malloc(m * n * sizeof(float));

x = (float*) malloc(n*sizeof(float));

y = (float*) malloc(m*sizeof(float));

memcpy(mat, data, m * n * sizeof(float));

memcpy(x, data, n * sizeof(float));

memcpy(y, data, m * sizeof(float));

timeg = 0.;
  
for (int i = 0; i < nbIter; ++i)
{
begin();

  cblas_sgemv(CblasColMajor, CblasTrans, n, m, 1, mat, n, x, 1, 1, y, 1);

end();
timeg += getsec();

}
printf("CPU Time: %10.8f (secs)\n", timeg/(double)nbIter );

free(mat);

free(x);

free(y);

free(data);

}
}

