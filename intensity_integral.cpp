#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include <Vector3d.h>
#include <Msis21.h>
#include <vector_gauss_quadrature.h>


#include "Fixed_parameter.h"
#include "Date.h"
#include "Across_pts.h"
#include "Geocoordinate.h"
#include "pmc_simulation.h"


double mie_beta(AndoLab::Vector3d <double> r, void *parameters);


/* rの点に太陽光が入射しているか */
bool is_illuminated(AndoLab::Vector3d <double> r, Date date){

  double rQ = AndoLab::abs( r_s(date) * (r * r_s(date)) );

  return ( r%r_s(date) >= 0.0 ) ||
      ( (r%r_s(date) < 0.0) && (rQ >= Radius_of_Earth) );
}

bool is_in_shadow_zone(AndoLab::Vector3d <double> r, Date date){
  return !is_illuminated(r,date);
}

/* r の点から太陽方向へ移動して大気圏界面に達する点 */
AndoLab::Vector3d <double> upper_atmosphere_to_solar(
    AndoLab::Vector3d <double> r, Date date){

  return cross_point_at_altitude(r, r_s(date).n(),
      Altitude_of_Atmosphere );
}

bool is_inside_PMC_region(AndoLab::Vector3d <double> r){
  return ( (r.abs() > R_lower_pmc) && (r.abs() < R_upper_pmc) );
}

constexpr double Reflection_Coefficient { 0.3 };
constexpr int M_alpha { 30 };
constexpr int N_beta { 30 };
constexpr double Dbeta { 2.0*M_PI / N_beta };
constexpr double R0 { Radius_of_Earth };


double intensity_integral(
    AndoLab::Vector3d <double> &r_from,
    AndoLab::Vector3d <double> &r_to,
    Date &date,
    const double lambda,
    const double th, /* 散乱角 */
    AndoLab::Msis21 &msis,
    double &delta,
    PMC &pmc,
    Parameters &param,
    const int region,
    const double interval){

  std::string function_name( "[intensity_integral] " );
  bool pmc_calc = pmc.is_this_real() && (region == INSIDE_PML_LAYER) && This_position_is_near_PMC( pmc, r_from, r_to );

  double I { 0.0 };

  /* 寄与を計算するベクトル */
  AndoLab::Vector3d <double> length = r_to - r_from;

  /* 分割点数。interval毎に光学的深さを計算する。PMCがあれば細かな積分区間とする */
  const int Num = int( length.abs() / interval + 1 );
  AndoLab::Vector3d <double> dr = length / double(Num); /* 分割距離 */

  AndoLab::Vector3d <double> pre_rp = r_from;

  /* 反射のための各種設定 */
  AndoLab::Vector3d <double> yv( 0.0, 1.0, 0.0 );
  AndoLab::Vector3d <double> zv( 0.0, 0.0, 1.0 );
  AndoLab::Vector3d <double> r_gs = (r_from - r_to).n();

  for(int i_dr = 0; i_dr < Num; i_dr++){
//    std::cout << "intensity_integral : " << i_dr << " / " << Num << std::endl;
    AndoLab::Vector3d <double> r_p = r_from + (i_dr+0.5)*dr;

    /* 視線方向の光学的深さの計算
     * 視線から遠い方向へと計算してゆけば
     * 先の計算が使える */

    delta += vector_gauss_quadrature(beta, (void*)(&param), pre_rp, r_p); /* 大気 */
    if ( pmc_calc ){
      delta += vector_gauss_quadrature(mie_beta, (void*)(&pmc), pre_rp, r_p); /* 大気 */
    }

    if ( is_illuminated(r_p, date) ){

      /* r の点がシャドウ領域でなければ
       * 1. 太陽方向の大気圏界面の点 r0 を見つける
       * 2. exp( - delta(r0,r) - delta(r,r2geo) ) * σ(θ) */
      AndoLab::Vector3d <double> r_i =
          upper_atmosphere_to_solar(r_p, date);

      /* 大気圏内に入射から散乱点までの光学的深さ (レイリー) */
      double optical_depth = vector_gauss_quadrature(beta, (void*)(&param), r_i, r_p);

      /* 大気圏内に入射から散乱点までの光学的深さ (PMC) */
      /* 2024/06/20 → 正確に計算したが、殆ど影響なし。誤差は 0.01%以下 */
      if ( pmc_calc ){
        AndoLab::Vector3d <double> r_pml_l[2], r_pml_u[2]; /* 探索するPMC層の上下限 */
        AndoLab::Vector3d <double> integral_from[2], integral_to[2]; /* 探索結果に基づいた積分範囲 */
        int Num_integral_region { 0 }; /* 積分すべき領域数 */
        Across_pts PtsAcrossLowerPMC(r_p, r_i, Lower_PMC);
        Across_pts PtsAcrossUpperPMC(r_p, r_i, Upper_PMC);

        if ( region == LOWER_THAN_PML_LAYER ){
          /* PMC領域は1つのみ(上下限が一つずつ) */
          if ( (PtsAcrossLowerPMC.num == 1) && (PtsAcrossUpperPMC.num == 1) ){
            Num_integral_region = 1;
            integral_from[0] = PtsAcrossLowerPMC.r[0];
            integral_to[0] = PtsAcrossUpperPMC.r[0];
          } else {
            std::cerr
            << "Error (" << process_id
            << ", intensity_integral) : integral region is not correctly found! (Lower than PMC layer)"
            << "[Lower = " << PtsAcrossLowerPMC.num
            << ", Upper = " << PtsAcrossUpperPMC.num << "]"
            << std::endl;
            exit(0);
          }

        } else if ( region == INSIDE_PML_LAYER ){
          /* PMC領域は1つ、または2つ
           * 1つのときは 上限が一つ見つかる → PMC層中の点からPMC層上限で大気圏を抜ける */
          if ( (PtsAcrossLowerPMC.num == 0) && (PtsAcrossUpperPMC.num == 1) ){
            Num_integral_region = 1;
            integral_from[0] = r_p;
            integral_to[0] = PtsAcrossUpperPMC.r[0];

          } else if ( (PtsAcrossLowerPMC.num == 2) && (PtsAcrossUpperPMC.num == 1) ){
            /* 2つのときは 下限が二つ、上限が一つ見つかる
             * → PMC層中の点からPMC層下限を通って、PMC層下の領域へ行き、
             * 別のPMC層下限から上限へ抜けて大気圏を抜ける */
            Num_integral_region = 2;
            integral_from[0] = r_p;
            integral_to[0]   = PtsAcrossLowerPMC.r[0];

            integral_from[1] = PtsAcrossLowerPMC.r[1];
            integral_to[1]   = PtsAcrossUpperPMC.r[0];
          } else {
            std::cerr
            << "Error (intensity_integral) : integral region is not correctly found! (inside PML layer) "
            << "[Lower = " << PtsAcrossLowerPMC.num
            << ", " << PtsAcrossUpperPMC.num << "]"
            << std::endl;
            exit(0);
          }

        } else {/* HIGHER_THAN_PML_LAYER */
          if ( PtsAcrossUpperPMC.num == 0) {
            /* 通常はPMC層を通らず上に抜ける */

          } else if ( (PtsAcrossLowerPMC.num == 0) && (PtsAcrossUpperPMC.num == 2) ){
            /* 領域1つ：一旦PMC層に入り、そのまま抜けてゆく。上限が二つ見つかる */
            Num_integral_region = 1;
            integral_from[0] = PtsAcrossUpperPMC.r[0];
            integral_to[0] = PtsAcrossUpperPMC.r[1];

          } else if ( (PtsAcrossLowerPMC.num == 2) && (PtsAcrossUpperPMC.num == 2) ){
            /* 領域2つ：PMC層に入り、PMC層下の領域へ行って、再度PMC層を通って抜けてゆく。
             * 上下限が二つずつ見つかる */
            Num_integral_region = 2;
            integral_from[0] = PtsAcrossUpperPMC.r[0];
            integral_to[0] = PtsAcrossLowerPMC.r[0];
            integral_from[1] = PtsAcrossLowerPMC.r[1];
            integral_to[1] = PtsAcrossUpperPMC.r[1];

          } else {
            std::cerr
            << "Error (intensity_integral) : integral region is not correctly found! (Higher than PML layer) "
            << "[Lower = " << PtsAcrossLowerPMC.num
            << ", " << PtsAcrossUpperPMC.num << "]"
            << std::endl;
            exit(0);
          }

        }
        for(int n = 0; n < Num_integral_region; n++){
//          std::cout << "From " << ( integral_from[n].abs() - Radius_of_Earth )*m2km
//              << " To " << ( integral_to[n].abs() - Radius_of_Earth )*m2km << std::endl;
          optical_depth += vector_gauss_quadrature(mie_beta, (void*)(&pmc), integral_from[n], integral_to[n]);
        }
      }

      Geocoordinate gr_p( r_p );
      msis.calc_at( gr_p.altitude(), gr_p.latitude(), gr_p.longitude() );

      /* 散乱係数 */
      double total_beta = sigma(th, lambda, &msis) * msis.N();
      if ( pmc_calc ){
        total_beta += pmc.bn_th() * pmc.dist( r_p );
      }

      /* 反射の検討 */
      double reflection = 0.0; /* 単一分子当たりの反射光 */
      //      const double alpha_max = std::acos( R0 / r_p.r() );
      //      const double Dalpha = alpha_max / M_alpha;
      //      const double theta_s = r_p.theta();
      //      const double phi_s = r_p.phi();
      //
      //      for(int m = 0; m < M_alpha; m++){
      //        const double alpha = (m + 0.5) * Dalpha;
      //        const double Smn = R0*R0 * std::sin(alpha) * Dalpha * Dbeta;
      //
      //        for(int n = 0; n < N_beta; n++){
      //          AndoLab::Vector3d <double> rp_R_mn(R0, alpha, (n+0.5)*Dbeta, AndoLab::coordinate::Spherical);
      //          AndoLab::Vector3d <double> r_R_mn = ( rp_R_mn.rotate(theta_s, yv) ).rotate(phi_s, zv); /* 反射点 */
      //
      //          /* 反射点への太陽光入射判定 */
      //          if ( r_R_mn%r_s(date) <= 0.0 ){
      //            continue;
      //          }
      //
      //          AndoLab::Vector3d <double> r_Rs_mn = r_p - r_R_mn; /* 反射点→散乱点ベクトル */
      //          AndoLab::Vector3d <double> r_Ri_mn =
      //            upper_atmosphere_to_solar(r_R_mn, date); /* 反射点に入射する光の大気圏入射点 */
      //
      //          /* 反射点到達光 */
      //          double I_r_mn = std::exp( - Optical_depth(lambda, date, r_Ri_mn, r_R_mn, num_pmc, pmc) );
      //          double Ip_S_mn = Smn / M_PI / r_Rs_mn.abs() / r_Rs_mn.abs()
      //            * I_r_mn * r_R_mn.n()%r_p.n() * r_Rs_mn.n()%r_R_mn.n(); /* 無損失散乱光 */
      //          /* 損失散乱光 */
      //          double I_S_mn = Ip_S_mn * std::exp( - Optical_depth(lambda, date, r_p, r_R_mn, num_pmc, pmc) );
      //          double theta_s_mn = AndoLab::angle_between( r_Rs_mn, r_gs  ); /* 反射光散乱角 */
      //          double tmp_ = I_S_mn * sigma(theta_s_mn, lambda, msis);
      //          //reflection += I_S_mn * sigma(theta_s_mn, lambda, msis);
      //          reflection += tmp_;
      //
      //        }
      //      }

      /* 反射の検討(ここまで) */

      I += std::exp(-delta) *
          (std::exp( -optical_depth ) * total_beta + reflection * msis.N() )
          * dr.abs();
    }
    pre_rp = r_p;
  }

  delta += vector_gauss_quadrature(beta, (void*)(&param), pre_rp, r_to); /* 大気 */
  if ( pmc_calc ){
    delta += vector_gauss_quadrature(mie_beta, (void*)(&pmc), pre_rp, r_to); /* PMC */
  }

  return I;
}
