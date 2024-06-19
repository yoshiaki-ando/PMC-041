/*
 * get_arg.cpp
 *
 *  Created on: 2024/02/28
 *      Author: ando
 *
 * 引数から、解析する緯度・経度、日付・時刻を得る
 *
 * 引数の仕様
 * (実行ファイル) 緯度 経度 YYMMDD hhmm
 * argv[0]        [1]  [2]  [3]    [4]
 */


#include <iostream>
#include <string>

#include <Msis21.h>

#include "Date.h"
#include "pmc_simulation.h"

void get_arg(const int argc, char **argv, double &latitude, double &longitude,
    Date &date, AndoLab::Msis21 &msis){

  if ( argc < 4 ){
    std::cerr << "Insufficient arguments." << std::endl;
    exit(1);
  }

  /* 緯度、経度 設定 */
  latitude = std::stod( std::string( argv[1] ) );
  longitude = std::stod( std::string( argv[2] ) );


  std::string YYMMDD = argv[3];
  std::string hhmm = argv[4];

  int year = std::stoi( YYMMDD.substr(0,2) ); /* 使わない */
  int month = std::stoi( YYMMDD.substr(2,2) ); /* 月 */
  int day = std::stoi( YYMMDD.substr(4,2) );
  int hour = std::stoi( hhmm.substr(0,2) );
  int minute = std::stoi( hhmm.substr(2,2) );

  msis.date( 2000+year, month, day, hour*60 + minute );

  date.doy( msis.doy() );
  date.minute_day( hour*60 + minute );

  return;
}

//  std::string str_hour;
//  if ( argc > 1 ){
//    str_hour = argv[1];
//    Minute_of_Day = int( std::stod(str_hour) * 60 + 0.5 );
//  }

