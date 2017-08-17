/*
 * RKF45法で二階の常微分方程式を解く．
 * テーマ：二体問題．
 * 線形代数ライブラリのEigenを使用．
 * Boostも使用．
 *
 * by みーくん．
 */

#include "ERK.hpp"
#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <string>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <boost/multiprecision/cpp_dec_float.hpp>


namespace mp = boost::multiprecision;

using multiFloat = mp::cpp_dec_float_100;

// パラメータ
const multiFloat e("0.9999999999999999");

//時刻に関するパラメータ
multiFloat dt("1.0e-1");
const multiFloat t_limit("100.0");

const multiFloat RTol("10e-10");
const multiFloat ATol("10e-10");
const multiFloat t_min("10e-50");
const multiFloat t_max("0.1");

//インターバル
constexpr int INTV = 1;

//movie
constexpr int sim = 0;

int main(){
    Eigen::Matrix<multiFloat, 4, 1> x(1 - e, 0, 0, mp::sqrt((1 + e) / (1 - e)));
    Eigen::Matrix<multiFloat, 4, 1> x8, x9;
    //std::cout << std::fixed << std::setprecision(std::numeric_limits<multiFloat>::digits10 + 1);
    std::cout << std::fixed << std::setprecision(35);
    std::cerr << std::fixed << std::setprecision(35);

    multiFloat t{};

    mino2357::RKF78<multiFloat> rkf78(ATol, RTol);

    for(int i=0; t<t_limit; i++){
        if(i%INTV == 0){
            //std::cerr << t << " " << mp::log10(dt) << std::endl;
            std::cout << t << " " << mp::log10(dt) << " " << x(0,0) << " " << x(1,0) << " " << x(2,0) << " " << x(3,0) << std::endl;
        }
        rkf78.Integrate(t, dt, x);
    }

}
