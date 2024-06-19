/*
 * main.cpp
 *
 *  Created on: 2024/05/27
 *      Author: ando
 */
#include <iostream>
#include <fstream>
#include <filesystem>
#include <complex>

#include <Msis21.h>
#include <memory_allocate.h>
#include <Mie_scattering.h>

#include "Fixed_parameter.h"
#include "Date.h"
#include "Geocoordinate.h"
#include "pmc.h"
#include "pmc_simulation.h"
#include "pmc_observation.h"

std::ofstream ofs_log; /* ログ出力 */
std::ofstream ofs_param; /* 計算したパラメタの出力 */
bool logging { false };

int main(int argc, char **argv){

  wrap_msisinit_(); /* MSIS初期化 */

  std::filesystem::create_directory("data");
  ofs_param.open("data/param.txt", std::ios::out);
  ofs_log.open("data/log.txt", std::ios::out);

  /************************************************************
   * 初期化
   ************************************************************/

  Date date; /* (自作の日付クラス)解析をする日付、UT */
  double latitude, longitude; /* 解析をする(観測点の)緯度・経度 */

  AndoLab::Msis21 msis;

  /* 引数から、観測点の緯度・経度、(msisクラスへの)日付・時刻を取得 */
  get_arg(argc, argv, latitude, longitude, date, msis);

  /* 緯度・経度から、計算に使う座標(ひまわり方向を +x方向)への変換
   *
   * 緯度・経度・高度 0m の位置を r0 とする。地表面のtangential point
   */
  AndoLab::Vector3d <double> r0;
  convert_coordinate(latitude, longitude, r0);
  Geocoordinate Geo_r0( r0 ); /* ひまわり座標・地理座標を扱う変換 */
  std::cout << latitude << ", " << longitude << " ==> "
      << Geo_r0.latitude() << ", " << Geo_r0.longitude() << std::endl;

  double alpha = Geo_r0.alpha(); /* ひまわりから見た角度 */

  /************************************************************
   * 初期化（ここまで）
   ************************************************************/


  /* レイリー散乱のみの計算(こちらはフィッティングする緯度経度のみ計算) */
  double **Rayleigh = AndoLab::allocate_memory2d(Num_Lambda, N_alt, 0.0);

  PMC pmc[3]; /* ダミーPMCの設定 */
  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::ofstream ofs("data/rayleigh_" + std::to_string(j_lambda) + ".dat");
//    std::cout << "Lambda: " << j_lambda << std::endl;
    calculate_intensity(date, Lambda[j_lambda],
        Altitude_min, dAlt, N_alt,
        alpha, msis, pmc[j_lambda], Rayleigh[j_lambda]);

    for(int i = 0; i < N_alt; i++){
      ofs << (Altitude_min + i*dAlt) * 1e-3 << " " << Rayleigh[j_lambda][i] << "\n";
    }
    ofs.close();
  }

  /* 観測データの読み込み */
  /* 観測データの取得 */
  double ***Observed_data
  = AndoLab::allocate_memory3d(Num_Lambda, Num_observed_data_latitude, Num_observed_data_altitude, 0.0);
  double **Obsrvd_LatLon = AndoLab::allocate_memory2d(2, Num_observed_data_latitude, 0.0);
  get_observation_data(Observed_data, Obsrvd_LatLon);

  const int i_alpha = select_latlon( Obsrvd_LatLon, latitude, longitude );

  for(int jLambda = 0; jLambda < Num_Lambda; jLambda++){
    std::ofstream ofs_obs("data/obs_" + std::to_string(jLambda) + ".dat");
    for(int iAlt = 0; iAlt <= 100; iAlt++){
      ofs_obs << iAlt << " " << Observed_data[jLambda][i_alpha][iAlt] << std::endl;
    }
    ofs_obs.close();
  }

  /* フィッティング係数の導出 */
  double *SolarRayIntensity = new double [Num_Lambda];   /* 係数 = 太陽光の強さ */
  double *BackgroundIntensity = new double [Num_Lambda]; /* オフセット = 背景光の明るさ */

  ObtainFittingCoefficient(Rayleigh, i_alpha, Observed_data, SolarRayIntensity, BackgroundIntensity);

//  /* (一時的)フィッティングの様子を見るための出力 */
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    std::ofstream ofs("data/rayleigh_" + std::to_string(j_lambda) + ".dat");
//    calculate_intensity(date, Lambda[j_lambda],
//        0.0e3, dAlt, 100,
//        alpha, msis, pmc[j_lambda], intensity);
//
//    for(int i = 0; i < 100; i++){
//      ofs << i << " " << SolarRayIntensity[j_lambda]*intensity[i] + BackgroundIntensity[j_lambda] << "\n";
//    }
//    ofs.close();
//  }
//  /* (一時的なルーチン、ここまで) */

//  /* (一時的)適当な位置にPMCをおいて計算する */
//  double *intensity = new double [100];
//  logging = true;
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
//    AndoLab::Vector3d <double> rc = tangential_point( 83.0e3, alpha );
//    pmc[j_lambda].set( Lambda[j_lambda], m_ice, rc, 0.2e3, 300.e3, 80.0e-9, 1.4, 25.4e6 );
//
//    std::ofstream ofs("data/pmc_" + std::to_string(j_lambda) + ".dat");
//    calculate_intensity(date, Lambda[j_lambda],
//        0.0e3, dAlt, 100,
//        alpha, msis, pmc[j_lambda], intensity);
//
//    for(int i = 0; i < 100; i++){
//      ofs << i << " " << SolarRayIntensity[j_lambda]*intensity[i] + BackgroundIntensity[j_lambda] << "\n";
//    }
//    ofs.close();
//    logging = false;
//  }
//  delete [] intensity;
//  /* (一時的なルーチン、ここまで) */

  /* PMCパラメタの最適化をする */
  ObtainFittingPMC(date, pmc, msis, alpha, i_alpha, Observed_data, SolarRayIntensity, BackgroundIntensity);

  double *intensity = new double [100];
  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){

    std::ofstream ofs("data/optimized_pmc_" + std::to_string(j_lambda) + ".dat");
    calculate_intensity(date, Lambda[j_lambda],
        0.0e3, dAlt, 100,
        alpha, msis, pmc[j_lambda], intensity);

    for(int i = 0; i < 100; i++){
      ofs << i << " " << SolarRayIntensity[j_lambda]*intensity[i] + BackgroundIntensity[j_lambda] << "\n";
    }
    ofs.close();
  }
  delete [] intensity;

  AndoLab::deallocate_memory2d(Rayleigh);
  AndoLab::deallocate_memory3d(Observed_data);
  AndoLab::deallocate_memory2d(Obsrvd_LatLon);
  delete [] SolarRayIntensity;
  delete [] BackgroundIntensity;
  ofs_param.close();
  ofs_log.close();

  return 0;
}


