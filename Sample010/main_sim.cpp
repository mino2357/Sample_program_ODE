/*
 * RKF45法で二階の常微分方程式を解く．
 * テーマ：2重振り子．
 * 線形代数ライブラリのEigenを使用．
 * Boostも使用．
 *
 * by みーくん．
 */

#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <string>

#include <cstring>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <stdlib.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace mp = boost::multiprecision;

using multiFloat = mp::cpp_dec_float_100;

// パラメータ
const multiFloat e("0.99999999");
//初期角
const multiFloat x_init("3.0");
const multiFloat y_init("3.0");
//初角速度
const multiFloat x_dot_init("0.0");
const multiFloat y_dot_init("0.0");

//時刻に関するパラメータ
multiFloat dt("1.0e-20");
const multiFloat t_limit("20.0");

const multiFloat e_tol("10e-10");
const multiFloat t_min("10e-50");
const multiFloat t_max("0.1");

//インターバル
constexpr int INTV = 1;

//movie
constexpr int sim = 0;

//R^4からR^4への関数．
template <typename T = multiFloat>
Eigen::Matrix<T, 4, 1> func(const Eigen::Matrix<multiFloat, 4, 1>& u){
    T x       = u(0, 0);
    T x_dot   = u(1, 0);
    T y       = u(2, 0);
    T y_dot   = u(3, 0);

    T R = x * x + y * y;

    return Eigen::Matrix<T, 4, 1> {
        x_dot,
        - x / mp::sqrt(R * R * R),
        y_dot,
        - y / mp::sqrt(R * R * R)
    };
}

int main(){
    Eigen::Matrix<multiFloat, 4, 1> x(1 - e, 0, 0, mp::sqrt((1 + e) / (1 - e)));
    Eigen::Matrix<multiFloat, 4, 1> x4, x5;
    std::cout << std::fixed << std::setprecision(std::numeric_limits<multiFloat>::digits10 + 1);

    Eigen::Matrix<multiFloat ,4, 1> k1, k2, k3, k4, k5, k6;

    multiFloat t{};

    for(std::size_t i{}; t<t_limit; ++i){
    
        if(i%INTV == 0){
            std::cout << t << " " << mp::log10(dt) << " " << x(0, 0) << " " << x(1, 0) << " " << x(2, 0) << " " << x(3, 0) << std::endl;
        }

        //RK4法で常微分方程式を解く．
        k1 = func<>(x);
        k2 = func<>(x + dt / 4 * k1);
        k3 = func<>(x + dt / 32 * (3 * k1 + 9 * k2));
        k4 = func<>(x + dt / 2197 * (1932 * k1 - 7200 * k2 + 7296 * k3));
        k5 = func<>(x + dt * (static_cast<multiFloat>(439)/216 * k1 - 8 * k2 + static_cast<multiFloat>(3680)/513 * k3 - static_cast<multiFloat>(845)/4104 * k4));
        k6 = func<>(x + dt * (- static_cast<multiFloat>(8)/27 * k1 + 2 * k2 - static_cast<multiFloat>(3544)/2565 * k3 + static_cast<multiFloat>(1859)/4104 * k4 - static_cast<multiFloat>(11)/40 * k5));
        
        x4 = x + dt * (static_cast<multiFloat>(25)/216 * k1 + static_cast<multiFloat>(1408)/2565 * k3 + static_cast<multiFloat>(2197)/4104 * k4 - static_cast<multiFloat>(1)/5 * k5);
        x5 = x + dt * (static_cast<multiFloat>(16)/135 * k1 + static_cast<multiFloat>(6656)/12825 * k3 + static_cast<multiFloat>(28561)/56430 * k4 - static_cast<multiFloat>(9)/50 * k5 + static_cast<multiFloat>(2)/55 * k6);
       
        auto temp  = (static_cast<multiFloat>(25)/216 * k1 + static_cast<multiFloat>(1408)/2565 * k3 + static_cast<multiFloat>(2197)/4104 * k4 - static_cast<multiFloat>(1)/5 * k5) - (static_cast<multiFloat>(16)/135 * k1 + static_cast<multiFloat>(6656)/12825 * k3 + static_cast<multiFloat>(28561)/56430 * k4 - static_cast<multiFloat>(9)/50 * k5 + static_cast<multiFloat>(2)/55 * k6);
        auto R = mp::sqrt(temp(0, 0) * temp(0, 0) + temp(1, 0) * temp(1, 0) + temp(2, 0) * temp(2, 0) + temp(3, 0) * temp(3, 0));

        auto delta = static_cast<multiFloat>("0.86") * static_cast<multiFloat>(mp::pow(static_cast<multiFloat>(e_tol/(2 * R)), static_cast<multiFloat>("0.25")));

        //この辺は適当にチューニング．
        
        if(delta < static_cast<multiFloat>("0.1")){
            dt = static_cast<multiFloat>("0.1") * dt;
            if(dt < t_min){
                dt = t_min;
            }
        }else if(delta > 4){
            dt = 4 * dt;
            if(dt > static_cast<multiFloat>("0.1")){
                dt = t_max;
            }
        }else{
            dt = delta * dt;
            if(dt > static_cast<multiFloat>("0.1")){
                dt = t_max;
            }
            if(dt < t_min){
                dt = t_min;
            }
        }

        //std::cout << dt << std::endl;

        x  = x4;

        t += dt;

        //std::cout << x4 << std::endl;
    }
    
    //pclose(gp);
}
