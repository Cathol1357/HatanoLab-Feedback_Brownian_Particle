#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include<time.h>


//確率epsilonで粒子の位置を正しく測定し、
//確率1 - epsilonで粒子の位置を間違えて測定してしまう関数
double binary_measurement(double Xt, double epsilon){
    double error = 0.05; //測定エラーの大きさ
    double random_value = (double)rand() / RAND_MAX; //0.0から1.0までの乱数を作成
    //イプシロン以下であれば正しい値を測定できる
    if (random_value < epsilon) {
        return Xt;  //正しい値
    } else {
        if(rand() % 2 == 0){
            return Xt - error;
        }else{
            return Xt + error;
        }
    }
}