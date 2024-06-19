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

int Num_iteration { 0 };

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
};

double ErrPMC(const std::vector <double> &Optimized_param, std::vector <double> &grad, void *Params){

  double *intensity = new double [Num_PMC_Layer];

  PMC_FittingParameter *p = (PMC_FittingParameter*)Params;

  AndoLab::Vector3d <double> rc = tangential_point( Optimized_param[0], p->alpha );

  double sq_err { 0.0 }; /* 二乗誤差 */
  const int IdxLowerAlt { int( std::round( Lower_PMC / dAlt ) ) };

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
    p->pmc[j_lambda].set( Lambda[j_lambda], m_ice, rc,
        Optimized_param[1],
        Optimized_param[2],
        Optimized_param[3],
        Optimized_param[4],
        Optimized_param[5]);
    calculate_intensity( *( p->ptr_date ), Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
        p->alpha, *(p->ptr_msis), p->pmc[j_lambda], intensity);

//    std::ofstream ofs("data/diff_" + std::to_string(j_lambda) + ".dat");
    for(int idx_alt = 0; idx_alt < Num_PMC_Layer; idx_alt++){
      double log_diffs = std::log( p->Observed_data[j_lambda][p->i_alpha][IdxLowerAlt + idx_alt] )
          - std::log( p->SolarRayIntensity[j_lambda] * intensity[idx_alt] + p->BackgroundIntensity[j_lambda] );
      sq_err += log_diffs * log_diffs;
//      ofs << ( Lower_PMC + idx_alt*dAlt ) * m2km << " " << log_diffs
//          << " " << std::log( p->Observed_data[j_lambda][p->i_alpha][IdxLowerAlt + idx_alt] )
//          << " " << std::log( p->SolarRayIntensity[j_lambda] * intensity[idx_alt] + p->BackgroundIntensity[j_lambda] )
//          << "\n";
    }
//    ofs.close();
  }

  delete [] intensity;
  std::cout << Num_iteration  << " : err = " << std::sqrt( sq_err )
      << " : Alt = " << ( p->pmc[0].rc().abs() - Radius_of_Earth ) * m2km << " km, "
      << "sig_z = " << p->pmc[0].sig_z() * m2km << " km, "
      << "sig_r = " << p->pmc[0].sig_r() * m2km << " km, "
      << "r0 = " << p->pmc[0].r0() * 1e9 << " nm, "
      << "sig = " << p->pmc[0].sigma() << " , "
      << "N0 = " << p->pmc[0].N0() * 1e-6 << " x10^6 m^-3"
      << std::endl;
  ofs_log << Num_iteration  << " : err = " << std::sqrt( sq_err )
      << " : Alt = " << ( p->pmc[0].rc().abs() - Radius_of_Earth ) * m2km << " km, "
      << "sig_z = " << p->pmc[0].sig_z() * m2km << " km, "
      << "sig_r = " << p->pmc[0].sig_r() * m2km << " km, "
      << "r0 = " << p->pmc[0].r0() * 1e9 << " nm, "
      << "sig = " << p->pmc[0].sigma() << " , "
      << "N0 = " << p->pmc[0].N0() * 1e-6 << " x10^6 m^-3"
      << std::endl;

  Num_iteration++;
  return std::sqrt( sq_err );
}

void ObtainFittingPMC(
    Date date,
    PMC *pmc,
    AndoLab::Msis21 &msis,
    const double alpha,
    const int i_alpha,
    double ***Observed_data,
    const double *SolarRayIntensity,
    const double *BackgroundIntensity){

  /********** PMCの最適化 **********/

  /*
   * アルゴリズムと、パラメタ数の設定
   * Optimized_param[0] : [m] 高度 Lower_PMC - Upper_PMC
   * Optimized_param[1] : [m] σ_z 0.1km - 5.0km
   * Optimized_param[2] : [m] σ_h (horizontal) 10km - 1500km
   * Optimized_param[3] : [m] r0 モード半径 10nm - 180nm
   * Optimized_param[4] : [-] σ 粒径分布の標準偏差 0.8-2.0
   * Optimized_param[5] : [m^-3] N0 粒子数 1e5 - 200e6
   *
   */
  constexpr int NUM_OPTIMIZED_PARAMETER { 6 };

  std::vector <double> Optimized_param(NUM_OPTIMIZED_PARAMETER);
  Optimized_param[0] = 82.47e3;
  Optimized_param[1] = 1.13766e3;
  Optimized_param[2] = 750.0e3;
  Optimized_param[3] = 67.3754e-9;
  Optimized_param[4] = 1.39474;
  Optimized_param[5] = 21.33e6;

  /* 最大値・最小値の設定 */
  std::vector <double> lb { Lower_PMC, 0.1e3, 10.0e3, 10.0e-9, 0.8, 1.0e5 };
  std::vector <double> ub { Upper_PMC, 5.0e3, 750.0e3, 180.0e-9, 2.0, 200.0e6 };

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

//  ErrPMC( Optimized_param, grad, (void*)(&Param) ); /* お試し計算 */

  /* 最適化 */
  nlopt::opt opt( nlopt::LN_NELDERMEAD, NUM_OPTIMIZED_PARAMETER );
  opt.set_min_objective( ErrPMC, (void*)(&Param) );
  opt.set_xtol_rel( 1.0e-2 );
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  double minf;
  nlopt::result result = opt.optimize( Optimized_param, minf );

  ofs_log << "(pmc class)"
      << " : Alt = " << ( pmc[0].rc().abs() - Radius_of_Earth ) * m2km << " km, "
      << "sig_z = " << pmc[0].sig_z() * m2km << " km, "
      << "sig_r = " << pmc[0].sig_r() * m2km << " km, "
      << "r0 = " << pmc[0].r0() * 1e9 << " nm, "
      << "sig = " << pmc[0].sigma() << " , "
      << "N0 = " << pmc[0].N0() * 1e-6 << " x10^6 m^-3"
      << std::endl;

  ofs_log << "(opt param)"
      << " : Alt = " << Optimized_param[0] * m2km << " km, "
      << "sig_z = " << Optimized_param[1] * m2km << " km, "
      << "sig_r = " << Optimized_param[2] * m2km << " km, "
      << "r0 = " << Optimized_param[3] * 1e9 << " nm, "
      << "sig = " << Optimized_param[4] << " , "
      << "N0 = " << Optimized_param[5] * 1e-6 << " x10^6 m^-3"
      << std::endl;

  /* 最適化された値をPMCに設定する */
  AndoLab::Vector3d <double> rc = tangential_point( Optimized_param[0], alpha );

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
    pmc[j_lambda].set( Lambda[j_lambda], m_ice, rc,
        Optimized_param[1], /* σ_z */
        Optimized_param[2], /* σ_h */
        Optimized_param[3], /* r0 */
        Optimized_param[4], /* σ(対数正規分布) */
        Optimized_param[5]); /* N0 */
  }

}


