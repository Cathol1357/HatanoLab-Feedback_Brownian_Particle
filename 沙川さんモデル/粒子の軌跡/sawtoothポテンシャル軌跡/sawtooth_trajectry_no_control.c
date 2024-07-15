#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include<time.h>

// 正規乱数を生成する関数
double generateGaussianNoise(double mean, double stdDev) {
    static int haveSpare = 0;
    static double spare;
    
    if (haveSpare) {
        haveSpare = 0;
        return mean + stdDev * spare;
    }
    
    haveSpare = 1;
    double u, v, s;
    do {
        u = (rand() / ((double) RAND_MAX)) * 2.0 - 1.0;
        v = (rand() / ((double) RAND_MAX)) * 2.0 - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    
    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return mean + stdDev * (u * s);
}

// 粒子の位置を周期的に補正する関数
double periodic_position(double x, double L) {
    while (x < -L / 2.0) x += L;
    while (x > L / 2.0) x -= L;
    return x;
}

double saw_tooth_potential(double x, int lambda) {
    double L = 1.0;                 // 粒子のいる箱の長さ
    double l = (3.0 * L) / 10.0;    // ポテンシャルの山がある位置を決定する定数
    double K = 1.0;                 // ポテンシャルの高さ

    //周期的境界条件
    x = periodic_position(x, L);

    // 一パターン目
    if (lambda == 0) {
        if (-(L / 2.0) <= x && x < -(L / 2.0) + l) {
            return K * (x + (L / 2.0)) / l;
        } else if (-(L / 2.0) + l <= x && x < L / 2.0) {
            return -K * (x - (L / 2.0)) / (L - l);
        } else {
            return -1;  // 範囲外の値
        }
    }

    // 2パターン目
    if (lambda == 1) {
        if (-(L / 2.0) <= x && x < -(L / 2.0) + l) {
            return -K * (x + (L / 2.0) - l) / (L - l);
        } else if (-(L / 2.0) + l <= x && x < -(L / 2.0) + 2.0 * l) {
            return K * (x + (L / 2.0) - l) /l;
        } else if (-(L / 2.0) + 2 * l <= x && x < L / 2.0) {
            return -K * (x - (L / 2.0) - l) / (L - l);
        } else {
            return -1;  // 範囲外の値
        }
    }

    return -1;  // lambdaが0でも1でもない場合のデフォルト値
}

double calc_sawtooth_potential_force(double x, double lambda){
    double L = 1.0;                 // 粒子のいる箱の長さ
    double l = (3.0 * L) / 10.0;    // ポテンシャルの山がある位置を決定する定数
    double K = 1.0;                 // ポテンシャルの高さ

    //周期的境界条件
    x = periodic_position(x, L);


    // 一パターン目
    if (lambda == 0) {
        if (x == -(L/2.0) || -(L / 2.0) <= x && x < -(L / 2.0) + l) {       //なぜか-L/2の時だけ値が条件式から外れてしまうので追加してみた
            return K/l;
        } else if (-(L / 2.0) + l <= x && x < L / 2.0) {
            return -K/l;
        } else {
            return -1;  // 範囲外の値
        }
    }

    // 2パターン目
    if (lambda == 1) {
        if (-(L / 2.0) <= x && x < -(L / 2.0) + l) {
            return -K / (L - l);
        } else if (-(L / 2.0) + l <= x && x < -(L / 2.0) + 2.0 * l) {
            return K / l;
        } else if (-(L / 2.0) + 2 * l <= x && x < L / 2.0) {
            return -K / (L - l);
        } else {
            return -1;  // 範囲外の値
        }
    }

    return -1;  // lambdaが0でも1でもない場合のデフォルト値
}

//ランダムに（ほぼ1/2の確率で）lambdaの値を変更するための関数
int Flashing_Rachet(){
    double random = generateGaussianNoise(0, 1);
    //今平均が0のガウス乱数を生成したので
    if(random < 0){
        return 0;
    }else{
        return 1;
    }
}

//測定結果yをみてコントロールパラメータを変化させる関数
int Feedback_control(double y, double L, double l){
    //なぜか<=にするとうまくif文の中に入ってくれないので放置
    if(-L/2.0 < y && y < -L/2.0 + l){
        printf("Changed!\n");
        return 1;
    }else{
        return 0;
    }
}


int main() {
    // パラメータ設定

    /*-------シミュレーション全体に関わるパラメータ-------*/

    double t_ini = 0.0;                  //シミュレーション開始時刻
    double t_fin = 1.0;                 //シミュレーション終了時刻
    int steps = 10000;                   // シミュレーションステップ数
    double dt = (t_fin - t_ini)/steps;   //時間刻み幅

    /*-------終了-------*/

    /*-------OU Processに関わるパラメータ-------*/

    double sigma = 0.5;         //ゆらぎの大きさ
    double X0 = 0.0;            //初期位置

    /*-------終了-------*/

    /*-------粒子の系設定に関わるパラメータ-------*/

    //今はランダムでラムダが切り替わるようにしているよ
    int lambda = 0;             //コントロールパラメータ（0か1）
    double L = 1.0;             //箱の長さ
    double l = (3.0 * L) / 10.0; //ポテンシャルの山の存在位置
    double force[steps];        //力の大きさ

    double x_min = - L / 2.0;  // xの最小値
    double x_max = L / 2.0;   // xの最大値
    double dx = L / steps;      //x刻みはば

    /*-------終了-------*/

    /*------フィードバック制御に関わるパラメータ------*/

    double y;       //測定結果
    double tau_0 = 0.05; //測定を行うタイミング（時間）
    double m = 0.0;     //ある時刻間で測定を行う回数
    int feedback_count = 0;         //フィードバックが起こった回数を記録

    /*-------フィードバックパラメータ終了-------*/

    //シミュレーション回数
    int num_simulations = 10;
    int bound_count = 0;            //境界を超えた回数をカウント（軌跡プロット用）


    // 乱数シードを初期化
    srand(time(NULL));


    /*-------複数回シミュレーションを開始-------*/


    for(int sim = 0; sim < num_simulations; sim++){
        feedback_count = 0;

        /*-------各シミュレーションごとにフォルダを作成-------*/

       const char *folder = "output";      // 出力フォルダ名
        char command[256];                 // フォルダ作成（既に存在する場合はエラーを無視）
        sprintf(command, "mkdir -p %s", folder);
        system(command);

        // 異なる乱数シードを設定
        srand(time(NULL) + rand());

        // ファイル名を作成
        char filename[256];
        sprintf(filename, "%s/output%04d.d", folder, sim);

        // ファイルオープン
        FILE *fp = fopen(filename, "w");
        if (fp == NULL) {
            fprintf(stderr, "ファイル %s を開くことができませんでした。\n", filename);
            return 1;
        }

        /*-------フォルダ作成終了-------*/



        /*-------シミュレーションごとに初期化-------*/

        double Xt = X0;         //初期位置 
        m = 0;                  //測定を行った回数
        lambda = 0;
        bound_count = 0;

        /*-------初期化処理終了-------*/





        /*-------シミュレーション実行-------*/

        for (int i = 0; i <= steps - 1; i++) {  //L/2を超えないようにしてみた


            double t = i * dt;

            //位置Xに対して周期的境界条件をかける

            if(Xt < x_min ){
                Xt += L;
                //境界カウントを-1する
                bound_count--;
            }else if(x_max  < Xt){
                Xt -= L;
                //境界カウントを+1する
                bound_count++;
            }
            
            //lambdaは0で一定にする

            force[i] = calc_sawtooth_potential_force(Xt, lambda);       //力を計算

            //データを出力
            //fprintf(fp, "%.6f %.6f \n", t, Xt);
            //x全範囲バージョンのデータ
            fprintf(fp, "%.6f %.6f \n", t, Xt + (bound_count) * L);

            double dWt = sqrt(dt) * generateGaussianNoise(0, 1);        //ガウス乱数を計算
            Xt = Xt - force[i] * dt + sigma * dWt;      //位置を更新




        }

        // ファイルクローズ
        fclose(fp);

        /*-------シミュレーション終わり-------*/

        printf("シミュレーション結果が '%s' に出力されました。\n", filename);
    }

    /*-------複数回シミュレーション終了-------*/


    return 0;
}