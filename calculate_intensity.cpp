/*
 * レイリー散乱の計算
 */
#include <iostream>
#include <fstream>
#include <string>

#include "Fixed_parameter.h"
#include "Across_pts.h"
#include "pmc_simulation.h"

constexpr int j_pmc { 0 }; /* 将来的な複数PMC対応のための予備 */

void calculate_intensity(
    Date Day_of_Year,
    const double lambda,
    const double Lower_Altitude,
    const double Step_Altitude,
    const int Number_of_Altitude,
    const double alpha,
    AndoLab::Msis21 &msis,
    PMC &pmc,
    double *intensity             /* 返り値。各高度の光強度 */
){

  std::string function_name( "[calculate_intensity] " );

  /* ガウス積分に渡すクラス化したパラメタ
   * 本当はもっと上で定義できる */
  Parameters param;
  param.ptr_msis = &msis;  /* MSISE */
  param.wavelength = lambda;

  AndoLab::Vector3d <double> Solar_d = r_s(Day_of_Year); /* 太陽の方向 */

  for(int iAlt = 0; iAlt < Number_of_Altitude; iAlt++){
//    std::cout << iAlt << " / " << Number_of_Altitude << " : " << std::flush;
    intensity[iAlt] = 0.0;
    double alt = Lower_Altitude + iAlt*Step_Altitude;

    /* 高度alt, 北極からの角度αのtangential point */
    AndoLab::Vector3d <double> r_tangential = tangential_point(alt, alpha);
    /* 衛星からの視線と、大気圏上界との交点 */
    AndoLab::Vector3d <double> *Pts_atmos = Across_point_atmosphere(r_tangential);


    AndoLab::Vector3d <double> Pts_upper_pmc[2], Pts_lower_pmc[2];
    /* PMC領域上面との交点 */
    if ( alt < Upper_PMC ){
      Across_pts PtsAcrossUpperPMC(Pts_atmos[0], Pts_atmos[1], Upper_PMC);
      PtsAcrossUpperPMC.copy(Pts_upper_pmc);
    }
    /* PMC領域下面との交点 */
    if ( alt < Lower_PMC ){
      Across_pts PtsAcrossLowerPMC(Pts_atmos[0], Pts_atmos[1], Lower_PMC);
      PtsAcrossLowerPMC.copy(Pts_lower_pmc);
    }

    /* 散乱角
     * Pts_atmos[0] - Pts_atmos[1] は、衛星の方向
     * -Solar_d は太陽光の入射方向
     */
    double th =
        angle_between( Pts_atmos[0] - Pts_atmos[1], -1.0*Solar_d);

    /* Mie散乱のための前処理(2) */
    if ( pmc.is_this_real() ){
        pmc.calc_normalized_beta_th(th);
    }

    if ( alt < Lower_PMC ){

//      std::cout << "Lower than PMC region." << std::endl;

      double delta2r = 0.0;
      /* 大気圏上部からPMC層上部までの計算
       * Pts_atmos[0]     -> Pts_upper_pmc[0] …散乱点はレイリー */
      intensity[iAlt] +=
          intensity_integral(Pts_atmos[0], Pts_upper_pmc[0],
              Day_of_Year, lambda,
              th, msis, delta2r, pmc, param, HIGHER_THAN_PML_LAYER, Integral_Interval);

      /* PMC層
       * Pts_upper_pmc[0] -> Pts_lower_pmc[0] …散乱点はPMC層内 */
//      std::cout << "region2" << std::endl;
      if ( pmc.is_this_real() &&  This_position_is_near_PMC( pmc, Pts_upper_pmc[0], Pts_lower_pmc[0] ) ){
        /* PMCがあるときは細かな刻みで積分する */
        intensity[iAlt] +=
            intensity_integral(Pts_upper_pmc[0], Pts_lower_pmc[0],
                Day_of_Year, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, PMC_Integral_Interval);
      } else {
        intensity[iAlt] +=
            intensity_integral(Pts_upper_pmc[0], Pts_lower_pmc[0],
                Day_of_Year, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, Integral_Interval);
      }

      /* PMC層下部から、PMC層下部まで
       * Pts_lower_pmc[0] -> Pts_lower_pmc[1] …散乱点はレイリー */
//      std::cout << "region3" << std::endl;
      intensity[iAlt] +=
          intensity_integral(Pts_lower_pmc[0], Pts_lower_pmc[1],
              Day_of_Year, lambda,
              th, msis, delta2r, pmc, param, LOWER_THAN_PML_LAYER, Integral_Interval);

      /* PMC層
       * Pts_lower_pmc[1] -> Pts_upper_pmc[1] …散乱点はPMC層内 */
//      std::cout << "region4" << std::endl;
      if ( pmc.is_this_real() &&  This_position_is_near_PMC( pmc, Pts_lower_pmc[1], Pts_upper_pmc[1] ) ){
        /* PMCがあるときは細かな刻みで積分する */
        intensity[iAlt] +=
            intensity_integral(Pts_lower_pmc[1], Pts_upper_pmc[1],
                Day_of_Year, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, PMC_Integral_Interval);
      } else {
        intensity[iAlt] +=
            intensity_integral(Pts_lower_pmc[1], Pts_upper_pmc[1],
                Day_of_Year, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, Integral_Interval);      }

      /* PMC層上部から大気圏上部まで
       * Pts_upper_pmc[1] -> Pts_atmos[1]     …散乱点はレイリー */
//      std::cout << "region5" << std::endl;
      intensity[iAlt] +=
          intensity_integral(Pts_upper_pmc[1], Pts_atmos[1],
              Day_of_Year, lambda,
              th, msis, delta2r, pmc, param, HIGHER_THAN_PML_LAYER, Integral_Interval);

    } else if ( alt < Upper_PMC ){ /* PMC内 */

      double delta2r = 0.0;
      /* 大気圏上部からPMC層上部までの計算
       * Pts_atmos[0]     -> Pts_upper_pmc[0] …散乱点はレイリー */
      intensity[iAlt] +=
          intensity_integral(Pts_atmos[0], Pts_upper_pmc[0],
              Day_of_Year, lambda,
              th, msis, delta2r, pmc, param, HIGHER_THAN_PML_LAYER, Integral_Interval);

      /* PMC層
       * Pts_upper_pmc[0] -> Pts_upper_pmc[1] …散乱点はPMC層内 */
      if ( pmc.is_this_real() &&  This_position_is_near_PMC( pmc, Pts_upper_pmc[0], Pts_upper_pmc[1] ) ){
        /* PMCがあるときは細かな刻みで積分する */
        intensity[iAlt] +=
            intensity_integral(Pts_upper_pmc[0], Pts_upper_pmc[1],
                Day_of_Year, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, PMC_Integral_Interval);
      } else {
        intensity[iAlt] +=
            intensity_integral(Pts_upper_pmc[0], Pts_upper_pmc[1],
                Day_of_Year, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, Integral_Interval);
      }

      /* PMC層上部から大気圏上部まで
       * Pts_upper_pmc[1] -> Pts_atmos[1]   …散乱点はレイリー   */
      intensity[iAlt] +=
          intensity_integral(Pts_upper_pmc[1], Pts_atmos[1],
              Day_of_Year, lambda,
              th, msis, delta2r, pmc, param, HIGHER_THAN_PML_LAYER, Integral_Interval);

    } else {


      double delta2r = 0.0;
      /* Pts_atmos[0] -> Pts_atmos[1]   …散乱点はレイリー   */
      intensity[iAlt] =
          intensity_integral(Pts_atmos[0], Pts_atmos[1],
              Day_of_Year, lambda,
              th, msis, delta2r, pmc, param, HIGHER_THAN_PML_LAYER, Integral_Interval);

    }

//    std::cout << std::endl;
  }

}
