#ifndef BINNING_H


//各区間に入っているデータの個数を計算する
void calculate_bins(double *data, int *bins, int num_data, double min_val, double max_val, int num_bins);

//相互情報量を計算する
double mutual_information(int *x_bins, int *y_bins, int num_data, int num_bins);

//指定したファイルに対して、指定した区間はば（num_bins）における相互情報量の値を計算する
double calculate_mutual_information(const char *filename, int num_bins);

#endif
