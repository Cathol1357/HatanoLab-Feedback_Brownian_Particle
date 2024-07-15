//ポアソンノイズを出力する関数
#include "poisson_noise.h"
#include <stdlib.h>
#include <math.h>

int poisson_noise(double lambda) {
    double L = exp(-lambda);
    double p = 1.0;
    int k = 0;

    do {
        k++;
        p *= (double)rand() / RAND_MAX;
    } while (p > L);

    return k - 1;
}
