/*-----------------------------------------------------------

バイナリーエラーを出す測定系に対して相互情報量と到達時間を計算するプログラム
全然動かない

--------------------------------------------------------*/

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

//確率epsilonで粒子の位置を正しく測定し、
//確率1 - epsilonで粒子の位置を間違えて測定してしまう関数
double binary_measurement(double Xt, double epsilon){
    double error = 0.05; //測定エラーの大きさ
    double random_value = (double)rand() / RAND_MAX; //0.0から1.0までの乱数を作成
    //イプシロン以下であれば正しい値を測定できる
    if (random_value < epsilon) {
        return Xt;  //正しい値
    } else {

        return -Xt;
    }
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
    double K = 3.0;                 // ポテンシャルの高さ

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
    double K = 3.0;                 // ポテンシャルの高さ

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

//lambdaの値を0なら1に、1なら0に変化させる関数
int Flashing_Rachet(int lambda){
    if(lambda % 2 == 0){
        return 1;
    }else{
        return 0;
    }
}

//測定結果yをみてコントロールパラメータを変化させる関数
int Feedback_control(double y, double L, double l){
    //なぜか<=にするとうまくif文の中に入ってくれないので放置
    if(-L/2.0 < y && y < -L/2.0 + l){
        //printf("Changed!\n");
        return 1;
    }else{
        return 0;
    }
}


int main(){
    // パラメータ設定

    /*-------シミュレーション全体に関わるパラメータ-------*/

    double t_ini = 0.0;                  //シミュレーション開始時刻
    double t_fin = 10.0;                 //シミュレーション終了時刻
    int steps = 1000;                   // シミュレーションステップ数
    double dt = (t_fin - t_ini)/steps;   //時間刻み幅

    /*-------終了-------*/

    /*-------OU Processに関わるパラメータ-------*/

    double sigma = 1.0;         //ゆらぎの大きさ
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

    double t;               //時間パラメータ



    /*-------終了-------*/

    /*------フィードバック制御に関わるパラメータ------*/

    double y;       //測定結果
    double tau_0 = 0.05; //測定を行うタイミング（時間）
    double m = 0.0;     //ある時刻間で測定を行う回数
    int feedback_count = 0;         //フィードバックが起こった回数を記録

    /*-------フィードバックパラメータ終了-------*/

    //シミュレーション回数
    int num_simulations =  100;      //1000にするとメモリ足りない
    int bound_count = 0;            //境界を超えた回数をカウント（軌跡プロット用）

    /*-------アンサンブル平均計算用-------*/

    //軌跡を保存しておく配列（分散と平均計算用）
    double X_trajectry[steps + 1][num_simulations];
    double V_trajectory[steps + 1][num_simulations];
    double variances[steps + 1];        //x^2のアンサンブル平均の計算結果を保存する
    double mean[steps + 1];         //粒子の位置のアンサンブル平均
    double velocity[steps + 1];     //粒子の速度のアンサンブル平均
    double ave_velocity;            //速度出力用


    //相互情報量
    double mutual_information = 0.0;
    //粒子を正しく測定できる確率（測定可能確率）（0 <= epsilon <= 1)
    double epsilon = 0.1;

    //初期到達時刻
    double ave_time = 0.0;
    double first_pass_time[num_simulations];
    for(int sim = 0; sim < num_simulations; sim++){
        first_pass_time[sim] = -1.0;
    }

    //目標地点
    double pass_position = 6.0;


    /*-------終了-------*/




    /*------バイナリーの場合の相互情報量を計算-------*/
    mutual_information =  - epsilon * log(epsilon) - (1 - epsilon) * log(1 - epsilon);

    /*-------終了-------*/


    // 乱数シードを初期化
    srand(time(NULL));


    /*-------複数回シミュレーションを開始-------*/


    for(int sim = 0; sim < num_simulations; sim++){
        feedback_count = 0;



        // 異なる乱数シードを設定
        srand(time(NULL) + rand());



        /*-------フォルダ作成終了-------*/



        /*-------シミュレーションごとに初期化-------*/

        double Xt = X0;                 //初期位置 
        X_trajectry[0][sim] = X0;      //初期位置を設定
        m = 0;                          //測定を行った回数
        t = 0.0;
        lambda = 0;
        bound_count = 0;


        /*-------初期化処理終了-------*/






        /*-------シミュレーション実行-------*/

        for (int i = 0; i <= steps; i++){ 
            //時間を更新
            t = i * dt;

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
            

            /*-------測定とフィードバック制御-------*/

            if(m * tau_0 < t && t < (m + 1) * tau_0){
                //測定を行う
                //測定において、確率epslonで本来の位置とは正負逆として判断してしまうものとする
                y = binary_measurement(Xt, epsilon);

                //フィードバックコントロールを行う
                lambda = Feedback_control(y, L, l);

                //フィードバックが起こった回数を記録
                if(lambda == 1) feedback_count++;

                m++;        //測定回数を増やす
            }

            /*-------測定とフィードバック制御終わり-------*/

            force[i] = calc_sawtooth_potential_force(Xt, lambda);       //力を計算

            //データを出力
            //fprintf(fp, "%.6f %.6f \n", t, Xt);
            //x全範囲バージョンのデータ
            //fprintf(fp, "%.6f %.6f \n", t, Xt + (bound_count) * L);

            double dWt = sqrt(dt) * generateGaussianNoise(0, 1);        //ガウス乱数を計算
            V_trajectory[i][sim] = Xt;                                   //一個前の位置を保存しておく

            Xt = Xt - force[i] * dt + sigma * dWt;      //位置を更新
            X_trajectry[i][sim] = Xt + (bound_count) * L;                   //位置を保存
            //printf("今の位置 : %.6f\n", X_trajectry[i][sim]);
            V_trajectory[i][sim] = (V_trajectory[i][sim] - Xt)/dt;          //速度を計算して保存する
            printf("今の速度 : %.6f\n", V_trajectory[i][sim]);

            //初期到達時刻を判定
            if(pass_position < X_trajectry[i][sim] && first_pass_time[sim] < 0){
                first_pass_time[sim] = t;
                printf("%d 回目の到達時間 : %.6f\n", sim + 1, first_pass_time[sim]);
            }
            


        }


        /*-------分散を計算する-------*/

        for (int i = 0; i <= steps; i++) {
            double t = i * dt;      //到達時間計算用
            //平均を初期化
            mean[i] = 0.0;
            //各刻み幅ごとに位置の値をたし上げる（一旦mean[i]に全部保存）
            for (int sim = 0; sim < num_simulations; sim++) {
                mean[i] += X_trajectry[i][sim];

                velocity[i] += V_trajectory[i][sim];
            }
        
            //シミュレーション数で割って位置のアンサンブル平均を出す
            mean[i] = (mean[i])/ num_simulations;


            velocity[i] = (velocity[i]/num_simulations);




            //x^2の値を初期化
            double variance = 0.0;
            for (int sim = 0; sim < num_simulations; sim++) {
                //x^2を、足し上げる
                variance += X_trajectry[i][sim] * X_trajectry[i][sim];
            }

        //シミュレーション数で割って分散の値を出す
        variance = variance / num_simulations;
        variances[i] = variance;
        }



        /*-------計算終了-------*/



        /*-------シミュレーション終わり-------*/

    }

    /*-------複数回シミュレーション終了-------*/

    /*-------平均到達時間計算-------*/

     int count = 0.0;
    for(int sim = 0;sim < num_simulations; sim++){

        if(first_pass_time[sim] > 0){
            ave_time += first_pass_time[sim];
            count++;
        }
    }

    
    ave_velocity = 0.0;

    for(int i = 0;i < steps;i++){
        ave_velocity += velocity[i];
    }

    ave_velocity /= steps;
    

    ave_time /= count;

    printf("目標地点に到達したシミュレーションの数 : %d\n", count);

    /*------平均時間計算終了-------*/

    /*-------ファイル出力準備-------*/

    /*-------フォルダ作成-------*/

    const char *folder = "output";      // 出力フォルダ名
    char command[256];                 // フォルダ作成（既に存在する場合はエラーを無視）
    sprintf(command, "mkdir -p %s", folder);
    system(command);

    /*
    // ファイル名を作成
    char filename[256];
    sprintf(filename, "%s/output_velocity.d", folder);

    // ファイルオープン
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "ファイル %s を開くことができませんでした。\n", filename);
        return 1;
    }


    fprintf(fp, "Mutual Information (Under feedback) I_c : %.6f", mutual_information);

    fclose(fp);
    */
    

    printf("\n シミュレーション回数 : %d\n", num_simulations);
    printf("相互情報量 : %.6f\n", mutual_information);
    printf("速度の平均 : %.6f\n", ave_velocity);
    printf("正しい値を測定できる確率 : %.6f\n", epsilon);
    printf("目標地点 %.1f までの平均到達時間 : %.6f\n",pass_position, ave_time);


    return 0;
}