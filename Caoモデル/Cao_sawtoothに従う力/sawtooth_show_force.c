/*----------------------------------------------------------------------------------------------------
Information and flux in a feedback controlled Brownian ratchet（Cao, 2009)に記載の
saw-toothポテンシャルを描くプログラム

周期境界条件V(x) = V(x + L)を課している。


------------------------------------------------------------------------------------------------------*/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>


// 粒子の位置を周期的に補正する関数
double periodic_position(double x, double L) {
    while (x < 0.0) x += L;
    while (x > L ) x -= L;
    return x;
}



double sawtooth_force(double x, double V_0, double a, double L){

    //周期的境界条件
    x = periodic_position(x, L);

    //条件を二個に分ける
    if(0.0 <= (x / L) && (x / L) <= a){
        printf("Detected!\n");
        return -(V_0/(a * L));
    }else if(a < (x/L) && (x/L) <= 1.0){
        return  (V_0/((1 - a) * L));
    }

    return -1;
}



int main() {
    int x_steps = 1000;
    int lambda = 0;

    double L = 1.0;      //箱の長さ
    double V_0 = 5.0;       //ポテンシャルの高さ
    double a = (1.0/3.0);       //ポテンシャルの山の位置
    double k_BT = 1.0;      //ボルツマン定数と温度

    double x_min = 0.0;  // xの最小値
    double x_max = 3.0;   // xの最大値
    double dx = (x_max - x_min) / x_steps;

    double force[x_steps];



    FILE *file = fopen("sawtooth_force.d", "w");

    if (file == NULL) {
        printf("ファイルを開くことができませんでした。\n");
        return 1;
    }

    for (int i = 0; i <= x_steps; i++) {        
        double x = x_min + i * dx;

        force[i] = sawtooth_force(x, V_0, a, L);
        fprintf(file, "%f %f\n", x, force[i]);
    }

    fclose(file);

    return 0;
}
