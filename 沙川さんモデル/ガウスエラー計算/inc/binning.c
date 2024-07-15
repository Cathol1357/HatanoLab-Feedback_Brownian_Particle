#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEFAULT_NUM_BINS 10

void calculate_bins(double *data, int *bins, int num_data, double min_val, double max_val, int num_bins) {
    double bin_width = (max_val - min_val) / num_bins;
    for (int i = 0; i < num_data; i++) {
        bins[i] = (int)((data[i] - min_val) / bin_width);
        if (bins[i] >= num_bins) bins[i] = num_bins - 1;
        if (bins[i] < 0) bins[i] = 0;
    }
}


double mutual_information(int *x_bins, int *y_bins, int num_data, int num_bins) {
    int **joint_histogram = (int **)malloc(num_bins * sizeof(int *));
    int *x_histogram = (int *)calloc(num_bins, sizeof(int));
    int *y_histogram = (int *)calloc(num_bins, sizeof(int));

    for (int i = 0; i < num_bins; i++) {
        joint_histogram[i] = (int *)calloc(num_bins, sizeof(int));
    }

    // ヒストグラムの計算
    for (int i = 0; i < num_data; i++) {
        joint_histogram[x_bins[i]][y_bins[i]]++;
        x_histogram[x_bins[i]]++;
        y_histogram[y_bins[i]]++;
    }

    double mi = 0.0;
    for (int i = 0; i < num_bins; i++) {
        for (int j = 0; j < num_bins; j++) {
            if (joint_histogram[i][j] > 0) {
                double p_xy = (double)joint_histogram[i][j] / num_data;
                double p_x = (double)x_histogram[i] / num_data;
                double p_y = (double)y_histogram[j] / num_data;
                mi += p_xy * log(p_xy / (p_x * p_y));
            }
        }
    }

    // メモリの解放
    for (int i = 0; i < num_bins; i++) {
        free(joint_histogram[i]);
    }
    free(joint_histogram);
    free(x_histogram);
    free(y_histogram);

    return mi / log(2);  // 2を底とするログを取るため
}


double calculate_mutual_information(const char *filename, int num_bins) {
    // ファイルからデータを読み込み
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return -1.0;
    }

    // ファイルからデータの行数をカウント
    int num_data = 0;
    char ch;
    while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') num_data++;
    }
    rewind(file);  // ファイルポインタを先頭に戻す

    // データを格納するためのメモリを動的に確保
    double *X = (double *)malloc(num_data * sizeof(double));
    double *Y = (double *)malloc(num_data * sizeof(double));
    if (X == NULL || Y == NULL) {
        perror("Error allocating memory");
        fclose(file);
        return -1.0;
    }

    // ファイルからデータを読み込み
    for (int i = 0; i < num_data; i++) {
        if (fscanf(file, "%lf %lf", &X[i], &Y[i]) != 2) {
            perror("Error reading data");
            free(X);
            free(Y);
            fclose(file);
            return -1.0;
        }
    }
    fclose(file);

    // XとYのビンを格納する配列
    int *x_bins = (int *)malloc(num_data * sizeof(int));
    int *y_bins = (int *)malloc(num_data * sizeof(int));
    if (x_bins == NULL || y_bins == NULL) {
        perror("Error allocating memory");
        free(X);
        free(Y);
        return -1.0;
    }

    // XとYの最小値と最大値を決定
    double x_min = X[0], x_max = X[0];
    double y_min = Y[0], y_max = Y[0];
    for (int i = 1; i < num_data; i++) {
        if (X[i] < x_min) x_min = X[i];
        if (X[i] > x_max) x_max = X[i];
        if (Y[i] < y_min) y_min = Y[i];
        if (Y[i] > y_max) y_max = Y[i];
    }

    // ビニングを計算
    calculate_bins(X, x_bins, num_data, x_min, x_max, num_bins);
    calculate_bins(Y, y_bins, num_data, y_min, y_max, num_bins);

    // 相互情報量を計算
    double mi = mutual_information(x_bins, y_bins, num_data, num_bins);

    // メモリの解放
    free(X);
    free(Y);
    free(x_bins);
    free(y_bins);

    return mi;
}





//一応メイン関数（動作チェック）
/*
int main() {
    const char *filename = "test.d";
    int num_bins = 10;

    double mi = calculate_mutual_information(filename, num_bins);
    printf("Mutual Information: %lf\n", mi);

    return 0;
}
*/

