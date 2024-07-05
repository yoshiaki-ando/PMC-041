/*
 * obtainfittingpmc.cpp
 *
 *  Created on: 2024/06/18
 *      Author: ando
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <nlopt.hpp>

#include <Msis21.h>
#include <memory_allocate.h>

#include "pmc_simulation.h"

AndoLab::Vector3d <double> rotate_alpha(AndoLab::Vector3d <double> r, const double alpha);
AndoLab::Vector3d <double> rotate_theta(AndoLab::Vector3d <double> r, const double theta);

AndoLab::Vector3d <double> set_center_pmc(
    const double pmc_altitude, const double horizontal_shift, const double alpha){

  const double R = Radius_of_Earth + pmc_altitude;
  double cos_b = R / Rgeo;
  double sin_b = std::sqrt( 1.0 - cos_b*cos_b );
  AndoLab::Vector3d <double> r(R*cos_b, 0.0, R*sin_b );
  return rotate_alpha( rotate_theta(r, horizontal_shift/R ), alpha );
}

extern std::ofstream ofs_log;

class PMC_FittingParameter{
private:
public:
  int i_alpha;
  double alpha;
  PMC *pmc;
  const Date *ptr_date;
  AndoLab::Msis21 *ptr_msis;
  const double *SolarRayIntensity;
  const double *BackgroundIntensity;
  double ***Observed_data;
  double **rayleigh_intensity;
};

constexpr double N0_base { 20.0e6 };

/* 同定に用いる評価関数。この誤差を最小化する */
double ErrPMC(const std::vector <double> &Optimized_param, std::vector <double> &grad, void *Params){

  double **intensity = AndoLab::allocate_memory2d(Num_Lambda, Num_PMC_Layer, 0.0);

  PMC_FittingParameter *p = (PMC_FittingParameter*)Params;

//  AndoLab::Vector3d <double> rc = tangential_point( Optimized_param[0], p->alpha ); /* PMCの中心 */
  AndoLab::Vector3d <double> rc = set_center_pmc(Optimized_param[0], Optimized_param[5], p->alpha); /* PMCの中心 */


  double sq_err { 0.0 }; /* 二乗誤差 */
  const int IdxLowerAlt { int( std::round( Lower_PMC / dAlt ) ) };

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
    p->pmc[j_lambda].set( Lambda[j_lambda], m_ice, rc,
        Optimized_param[1],
        Optimized_param[2],
        Optimized_param[3],
        Optimized_param[4],
        N0_base ); /* 密度は固定、比例係数を求める */
    /* PMCの高度のみを求める */
    calculate_intensity( *( p->ptr_date ), Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
        p->alpha, *(p->ptr_msis), p->pmc[j_lambda], intensity[j_lambda]);
  }

  /* 数密度係数を求める
   *
   */
  double nume { 0.0 };
  double deno { 0.0 };
  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    double nume_lambda { 0.0 };
    double deno_lambda { 0.0 };
    for(int idx_alt = 0; idx_alt < Num_PMC_Layer; idx_alt++){
      double I_PMC = intensity[j_lambda][idx_alt] - p->rayleigh_intensity[j_lambda][idx_alt];
      nume_lambda += I_PMC * ( p->Observed_data[j_lambda][p->i_alpha][IdxLowerAlt + idx_alt] -
          p->SolarRayIntensity[j_lambda] * p->rayleigh_intensity[j_lambda][idx_alt] -
          p->BackgroundIntensity[j_lambda] );
      deno_lambda += I_PMC * I_PMC;
    }
    nume += p->SolarRayIntensity[j_lambda] * nume_lambda;
    deno += p->SolarRayIntensity[j_lambda] * p->SolarRayIntensity[j_lambda] * deno_lambda;
  }
  double density_coefficient { nume / deno };
//  std::cout << density_coefficient << std::endl;

//    std::ofstream ofs("data/diff_" + std::to_string(j_lambda) + ".dat");
  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    /* 対数値の二乗誤差の計算 */
    for(int idx_alt = 0; idx_alt < Num_PMC_Layer; idx_alt++){
      double I_PMC = intensity[j_lambda][idx_alt] - p->rayleigh_intensity[j_lambda][idx_alt];
      double calc_result = p->SolarRayIntensity[j_lambda] *
          ( p->rayleigh_intensity[j_lambda][idx_alt] + density_coefficient * I_PMC ) +
          p->BackgroundIntensity[j_lambda];
      double log_diffs = std::log( p->Observed_data[j_lambda][p->i_alpha][IdxLowerAlt + idx_alt] / calc_result );
      sq_err += log_diffs * log_diffs;
//      ofs << ( Lower_PMC + idx_alt*dAlt ) * m2km << " " << log_diffs
//          << " " << std::log( p->Observed_data[j_lambda][p->i_alpha][IdxLowerAlt + idx_alt] )
//          << " " << std::log( p->SolarRayIntensity[j_lambda] * intensity[idx_alt] + p->BackgroundIntensity[j_lambda] )
//          << "\n";
    }
//    ofs.close();
  }

  AndoLab::deallocate_memory2d( intensity );

  p->pmc[0].N0( density_coefficient * N0_base ); /* 数密度係数で、N0を設定 */

  /* 最適化している様子を出力 */
  std::cout << "(" << process_id << ") "
      << number_of_iteration  << " "
      << std::sqrt( sq_err ) << " " /* Error */
      << ( p->pmc[0].rc().abs() - Radius_of_Earth ) * m2km << "km " /* Center Altitude */
      << p->pmc[0].sig_z() * m2km << "km " /* SD in altitude */
      << p->pmc[0].sig_r() * m2km << "km " /* SD in horizontal */
      << p->pmc[0].r0() * 1e9 << "nm " /* Particle radius */
      << p->pmc[0].sigma() << " "    /* SD parameter of log-normal radius distribution */
      << p->pmc[0].N0() * 1e-6 << "cm^-3 " /* Number density at the maximum (at the center) */
      << Optimized_param[5] * m2km << "km " /* horizontal shift */
      << std::endl;
  ofs_log << "(" << process_id << ") "
      << number_of_iteration  << " "
      << std::sqrt( sq_err ) << " " /* Error */
      << ( p->pmc[0].rc().abs() - Radius_of_Earth ) * m2km << " " /* Center Altitude */
      << p->pmc[0].sig_z() * m2km << " " /* SD in altitude */
      << p->pmc[0].sig_r() * m2km << " " /* SD in horizontal */
      << p->pmc[0].r0() * 1e9 << " " /* Particle radius */
      << p->pmc[0].sigma() << " "    /* SD parameter of log-normal radius distribution */
      << p->pmc[0].N0() * 1e-6 << " " /* Number density at the maximum (at the center) */
      << Optimized_param[5] * m2km << " " /* horizontal shift */
      << std::endl;

  number_of_iteration++;
  if ( number_of_iteration > Maximum_convergence ){
    throw nlopt::forced_stop();
  }
  return std::sqrt( sq_err );
}


void ObtainFittingPMC(
    /*
     ********* PMCの最適化 *********
     * NLOpt を用いて PMCパラメタを同定する
     */
    Date date,
    PMC *pmc,
    AndoLab::Msis21 &msis,
    const double alpha,
    const int i_alpha,
    double ***Observed_data,
    double **rayleigh_intensity,
    const double *SolarRayIntensity,
    const double *BackgroundIntensity,
    double &shift_distance){


  /*
   * アルゴリズムと、パラメタ数の設定
   * Optimized_param[0] : [m] 高度 Lower_PMC - Upper_PMC
   * Optimized_param[1] : [m] σ_z 0.1km - 5.0km
   * Optimized_param[2] : [m] σ_h (horizontal) 10km - 1500km
   * Optimized_param[3] : [m] r0 モード半径 10nm - 180nm
   * Optimized_param[4] : [-] σ 粒径分布の標準偏差 0.8-2.0
   * ×(高速化のため削除) Optimized_param[5] : [m^-3] N0 粒子数 1e5 - 200e6
   * Optimized_param[5] : [m] r_s  +y方向に対して右ねじにPMC中心位置をこの距離だけ回転(r_s = (R0 + z_pmc)*θ )
   * 0 - 750km
   *
   */
  constexpr int NUM_OPTIMIZED_PARAMETER { 6 };

  std::vector <double> Optimized_param(NUM_OPTIMIZED_PARAMETER);
  Optimized_param[0] = 82.5057e3;
  Optimized_param[1] = 1.10469e3;
  Optimized_param[2] = 750.0e3;
  Optimized_param[3] = 66.5536e-9;
  Optimized_param[4] = 1.39691;
//  Optimized_param[5] = 22.3356e6;
  Optimized_param[5] = 150e3;

  /* 最大値・最小値の設定 */
  std::vector <double> LowerBound { Lower_PMC, 0.1e3, 10.0e3,    10.0e-9, 0.8, 0.0 };
  std::vector <double> UpperBound { Upper_PMC, 5.0e3, 2000.0e3, 180.0e-9, 2.0, 1500.0e3 };

  std::vector <double> grad(NUM_OPTIMIZED_PARAMETER);

  /* わたすパラメタの設定 */
  PMC_FittingParameter Param;
  Param.alpha = alpha;
  Param.pmc = pmc;
  Param.ptr_date = &date;
  Param.ptr_msis = &msis;
  Param.SolarRayIntensity = SolarRayIntensity;
  Param.BackgroundIntensity = BackgroundIntensity;
  Param.Observed_data = Observed_data;
  Param.i_alpha = i_alpha;
  Param.rayleigh_intensity = rayleigh_intensity;

//  ErrPMC( Optimized_param, grad, (void*)(&Param) ); /* お試し計算 */

  /* NLOpt による最適化 */
  nlopt::opt opt( nlopt::LN_NELDERMEAD, NUM_OPTIMIZED_PARAMETER ); /* アルゴリズムと同定パラメタ数 */
  opt.set_min_objective( ErrPMC, (void*)(&Param) ); /*　評価関数と補助パラメタ */
  opt.set_xtol_rel( 1.0e-2 ); /* Torelance */
  opt.set_lower_bounds(LowerBound);
  opt.set_upper_bounds(UpperBound);

  double minf;

  try {
    number_of_iteration = 0;
    nlopt::result result = opt.optimize( Optimized_param, minf );
  } catch (std::exception &e){
    std::cout << "Optimization of PMC parameters (" << process_id << ") failed. : " << e.what() << std::endl;
    ofs_log << "Optimization of PMC parameters (" << process_id << ") failed. : " << e.what() << std::endl;
  }

  /* pmcクラスに残っているパラメタ */
  ofs_log << "# pmc class\n"
      << number_of_iteration  << " "
      << ( pmc[0].rc().abs() - Radius_of_Earth ) * m2km << " " /* Center Altitude [km] */
      << pmc[0].sig_z() * m2km << " " /* SD in altitude [km] */
      << pmc[0].sig_r() * m2km << " " /* SD in horizontal [km] */
      << pmc[0].r0() * 1e9 << " " /* Particle radius [nm] */
      << pmc[0].sigma() << " "    /* SD parameter of log-normal radius distribution [-] */
      << pmc[0].N0() * 1e-6 << " " /* Number density at the maximum (at the center) [x10^6 m^-3] */
      << std::endl;

  /* 最適化されたパラメタ */
  ofs_log << "# opt param\n"
      << Optimized_param[0] * m2km << " "
      << Optimized_param[1] * m2km << " "
      << Optimized_param[2] * m2km << " "
      << Optimized_param[3] * 1e9 << " "
      << Optimized_param[4] << " "
      << Optimized_param[5] * m2km << " "
      << std::endl;

  /* 最適化された値をPMCに設定する */
//  AndoLab::Vector3d <double> rc = tangential_point( Optimized_param[0], alpha );
  AndoLab::Vector3d <double> rc = set_center_pmc(Optimized_param[0], Optimized_param[5], alpha);

  const double N0 = pmc[0].N0();

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
    pmc[j_lambda].set( Lambda[j_lambda], m_ice, rc,
        Optimized_param[1], /* σ_z */
        Optimized_param[2], /* σ_h */
        Optimized_param[3], /* r0 */
        Optimized_param[4], /* σ(対数正規分布) */
        N0); /* N0 */
  }


}


