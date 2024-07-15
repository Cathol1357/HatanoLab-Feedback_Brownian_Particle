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

//saw_toothポテンシャルに従う値を返す
double saw_tooth_potential(double x, double V_0, double a, double L, double noize_level){


    //周期的境界条件
    x = periodic_position(x, L);

    //二つの場合に条件を分ける
    if(0.0 <= (x / L) && (x / L) <= a){
        return noize_level * (x * V_0)/(a * L);
    }else if(a < (x/L) && (x/L) <= 1.0){
        return (1 - noize_level) * (V_0 - (V_0 / (1 - a)) * ((x / L) - a));
    }

    return -1;  //変な事が起こった時に返す値
}



int main() {
    int x_steps = 1000;
    int lambda = 0;

    double L = 1.0;      //箱の長さ
    double V_0 = 5.0;       //ポテンシャルの高さ
    double a = (1.0/3.0);       //ポテンシャルの山の位置
    double k_BT = 1.0;      //ボルツマン定数と温度

    double noize_level = 0.0;       //ノイズの大きさを制御（論文中ではp、0 <= p < 1）

    double x_min = 0.0;  // xの最小値
    double x_max = 3.0;   // xの最大値
    double dx = (x_max - x_min) / x_steps;

    double potential[x_steps];
    double force[x_steps];

    //FILE *file = fopen("sawtooth_force.d", "w");


    FILE *file = fopen("sawtooth_potential.d", "w");

    if (file == NULL) {
        printf("ファイルを開くことができませんでした。\n");
        return 1;
    }

    for (int i = 0; i <= x_steps; i++) {        
        double x = x_min + i * dx;

        potential[i] = saw_tooth_potential(x, V_0, a, L, noize_level);
        fprintf(file, "%f %f\n", x, potential[i]);
    }

    fclose(file);

    return 0;
}
