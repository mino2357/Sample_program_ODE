/*
 * RKF法で二階の常微分方程式を解く．
 * テーマ：多体問題．
 * 線形代数ライブラリのEigenを使用．
 * Boostも使用．
 *
 * by みーくん．
 */


#include "param.hpp"
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

int main(){
    Eigen::Matrix<multiFloat, N, 1> x;
    // initial values
    x(0, 0) = 1;    x(1, 0) = 1;
    x(2, 0) = 4;    x(3, 0) = 1;
    x(4, 0) = 4;    x(5, 0) = 5;
    x(6, 0) = 0;    x(7, 0) = 0;
    x(8, 0) = 0;    x(9, 0) = 0;
    x(10, 0) = 0;   x(11, 0) = 0;
    
    Eigen::Matrix<multiFloat, N, 1> x8, x9;
    
    //std::cout << std::fixed << std::setprecision(std::numeric_limits<multiFloat>::digits10 + 1);
    std::cout << std::fixed << std::setprecision(6);
    std::cerr << std::fixed << std::setprecision(6);

    multiFloat t{};

    mino2357::RKF78<multiFloat> rkf78(ATol, RTol);

    for(int i=0; t<t_limit; i++){
        if(i%INTV == 0){
            std::cerr << t << " " << mp::log10(dt) << std::endl;
            for(int j=0; j<N; ++j){
                std::cout << x(j,0) << " ";
            }
            std::cout << std::endl;
        }
        rkf78.Integrate(t, dt, x);
    }
}
