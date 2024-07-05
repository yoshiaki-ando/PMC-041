/*
 * pmc_observation.h
 *
 *  Created on: 2024/06/07
 *      Author: ando
 */

#ifndef PMC_OBSERVATION_H_
#define PMC_OBSERVATION_H_

#include <string>

constexpr int Num_observed_data_altitude { 101 };  /* 高度方向のデータ数 */
constexpr int Num_observed_data_latitude { 44*2 }; /* 緯度のデータ数 */

constexpr double Altitude_observed_min { 0.0 };     /* 観測データの最低高度 */
constexpr double Altitude_observed_max { 100.0e3 }; /* 観測データの最高高度 */
constexpr double DAltitude_observed { 1.0e3 };      /* 観測データの高度刻み */

constexpr int Idx_background { 90 }; /* 背景光強度を求める最低高度のインデックス */

extern std::string obs_data_dir;

/*
 * 観測データを取得する
 */
void get_observation_data(
    double ***intensity,
    double **latlon,
    char **argv
    );

/*
 * 指定した緯度・経度に近いインデックスを取得
 */
int select_latlon(double **latlon, const double Latitude, const double Longitude);

/* フィッティング係数の導出 */
void ObtainFittingCoefficient(
    double **Rayleigh,          /* シミュレーションしたデータ */
    const int i_alpha,          /* 観測データのうち、使用する緯度・経度のインデックス */
    double ***Observed_data,    /* 観測データ */
    double *SolarRayIntensity,  /* 係数 = 太陽光の強さ */
    double *BackgroundIntensity /* オフセット = 背景光の明るさ */
    );

#endif /* PMC_OBSERVATION_H_ */
