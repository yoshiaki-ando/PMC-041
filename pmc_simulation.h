/*
 * pmc_simulation.h
 *
 *  Created on: 2024/05/27
 *      Author: ando
 */

#ifndef PMC_SIMULATION_H_
#define PMC_SIMULATION_H_

#include <string>

#include <Msis21.h>
#include <Vector3d.h>

#include "Fixed_parameter.h"
#include "Date.h"
#include "pmc.h"

/* 衛星視線上で大気圏上界面間を、分けて積分する区間 */
constexpr double Integral_Interval { 200.0e3 };
constexpr double PMC_Integral_Interval { 1.0e3 };


/*********************
 * 設定するパラメタ *
 *********************/

/* フィッティングパラメタを計算する高度 */
constexpr double Altitude_min { 20.e3 };
constexpr double Altitude_max { 43.e3 };
constexpr int N_alt { 24 }; /* 43 - 20 + 1 = 24 */
//constexpr double Altitude_min { 0.e3 };
//constexpr double Altitude_max { 100.e3 };
//constexpr int N_alt { 101 }; /* 43 - 20 + 1 = 24 */
constexpr double dAlt { (Altitude_max - Altitude_min) / (N_alt-1) };

constexpr int Maximum_convergence { 800 }; /* PMC最適化での反復回数上限 */

class Parameters{
private:
public:
  double wavelength;
  AndoLab::Msis21 *ptr_msis;
};

/* 収束しないときの反復回数 */
extern int number_of_iteration;

/* ジョブID */
extern std::string process_id;

/********************** 関数プロトタイプ **********************/

/* 引数から緯度、経度、日付を取得 */
void get_arg(const int argc, char **argv, double &latitude, double &longitude,
    Date &date,
    AndoLab::Msis21 &msis, std::string &data_dir,
    std::string &process_id);

/* 緯度・経度から、計算に使う座標(ひまわり方向を +x方向)へ変換 */
void convert_coordinate(
    const double latitude, /* [deg] */
    const double longitude, /* [deg] */
    AndoLab::Vector3d <double> &r0 /* [m] converted coordinate at altitude 0km */
    );

void calculate_intensity(
    Date Day_of_Year,
    const double Lambda,
    const double Lower_Altitude,
    const double Step_Altitude,
    const int Number_of_Altitude,
    const double alpha,
    AndoLab::Msis21 &msis,
    PMC &pmc,
    double *intensity             /* 返り値。各高度の光強度 */
    );

/* ひまわりから見て、北極から角度αで、高度Altitudeの球面とのtangential point
 * αが正なら東経140°以下側、負なら東経140°以上 */
AndoLab::Vector3d <double> tangential_point(double Altitude_in_m, double Angle_from_North_Pole_in_rad);

/* ひまわりから r の点を見た時、大気圏上空(高度100km)との交点2点を計算 */
AndoLab::Vector3d <double> *Across_point_atmosphere(AndoLab::Vector3d <double> r);

/* 点r から d方向にいったとき、高度 H で交差する点 */
AndoLab::Vector3d <double> cross_point_at_altitude(
    AndoLab::Vector3d <double> r,
    AndoLab::Vector3d <double> d, /* 単位ベクトルとすること */
    const double H /* [m] 高度 */
    );

/* 太陽を示すベクトル */
AndoLab::Vector3d <double> r_s(Date);

/*
 * r1 → r2 の経路中、pmc の分布と近くなる点が存在するか
 */
bool This_position_is_near_PMC(PMC &pmc, AndoLab::Vector3d <double> &r1, AndoLab::Vector3d <double> &r2 );

double intensity_integral(
    AndoLab::Vector3d <double> &r_from,
    AndoLab::Vector3d <double> &r_to,
    Date &date,
    const double Wavelength,
    const double th, /* 散乱角 */
    AndoLab::Msis21 &msis,
    double &delta,
    PMC &pmc,
    Parameters &param,
    const int Whether_inside_or_outside_PMC_region,
    const double integral_interval
    );

double beta(AndoLab::Vector3d <double> r, void* p);
double sigma_t(const double lambda, AndoLab::Msis21 *msis);
double sigma(const double th, const double lambda, AndoLab::Msis21 *msis);

inline double lat2theta(double latitude){
  return (90 - latitude)*Deg2rad;
}
inline double phi(double theta){
  return acos( Radius_of_Earth / Rgeo / sin(theta) );
}

void ObtainFittingPMC(
    Date date,
    PMC *pmc,
    AndoLab::Msis21 &msis,
    const double alpha,
    const int i_alpha,
    double ***Observed_data,
    double **rayleigh_intensity,
    const double *SolarRayIntensity,
    const double *BackgroundIntensity,
    double &shift_distance);

std::complex <double> interpolated_m(double wavelength_in_m);

#endif /* PMC_SIMULATION_H_ */
