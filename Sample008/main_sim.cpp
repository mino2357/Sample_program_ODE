/*
 * RKF45法で二階の常微分方程式を解く．
 * テーマ：2重振り子．
 * 線形代数ライブラリのEigenを使用．
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
const multiFloat g( "1.0");//9.8
const multiFloat m1("1.0");
const multiFloat m2("1.0");
const multiFloat L1("1.0");
const multiFloat L2("1.0");
//初期角
const multiFloat theta1_init("3.0");
const multiFloat theta2_init("3.0");
//初角速度
const multiFloat omega1_init("0.0");
const multiFloat omega2_init("0.0");

//時刻に関するパラメータ
multiFloat dt("0.0001");
const multiFloat t_limit("10.0");

const multiFloat e_tol("10e-10");
const multiFloat t_min("10e-50");

//インターバル
constexpr int INTV = 10;

//movie
constexpr int sim = 0;

//R^4からR^4への関数．
template <typename T = multiFloat>
Eigen::Matrix<T, 4, 1> func(const Eigen::Matrix<multiFloat, 4, 1>& x){
    T theta1 = x(0, 0);
    T eta1   = x(1, 0);
    T theta2 = x(2, 0);
    T eta2   = x(3, 0);

    return Eigen::Matrix<T, 4, 1> {
        eta1,
        (- m1 * g * mp::sin(theta1) - m2 * (g * mp::sin(theta1) + L2 * eta2 * eta2 * mp::sin(theta1 - theta2) + (L1 * eta1 * eta1 * mp::sin(theta1 - theta2) - g * mp::sin(theta2)) * mp::cos(theta1 - theta2))) / (L1 * (m1 + m2 * (mp::sin(theta1 - theta2) * (mp::sin(theta1 - theta2))))),
        eta2,
        ((m1 + m2) * (L1 * eta1 * eta1 * mp::sin(theta1 - theta2) - g * mp::sin(theta2) + g * mp::sin(theta1) * mp::cos(theta1 - theta2)) + m2 * L2 * eta2 * eta2 * mp::cos(theta1 - theta2) * mp::sin(theta1 - theta2)) / (L2 * (m1 + m2 * mp::sin(theta1 - theta2) * mp::sin(theta1 - theta2)))
    };
}

template <typename T = multiFloat>
T potentialEnergy(Eigen::Matrix<T, 4, 1>& x){
    T theta1 = x(0, 0);
    T eta1   = x(1, 0);
    T theta2 = x(2, 0);
    T eta2   = x(3, 0);

    return 0.5 * m1 * L1 * L1 * eta1 * eta1 + 0.5 * m2 * (L1 * L1 * eta1 * eta1 + L2 * L2 * eta2 * eta2 + 2. * L1 * L2 * eta1 * eta2 * mp::cos(theta1 - theta2));
}

template <typename T = multiFloat>
T kineticEnergy(Eigen::Matrix<T, 4, 1>& x){
    T theta1 = x(0, 0);
    T theta2 = x(2, 0);

    return - m1 * g * L1 * mp::cos(theta1) - m2 * g * (L1 * mp::cos(theta1) + L2 * mp::cos(theta2));
}

int main(){
    Eigen::Matrix<multiFloat, 4, 1> x(theta1_init, omega1_init, theta2_init, omega2_init);
    Eigen::Matrix<multiFloat, 4, 1> x4, x5;
    std::cout << std::fixed << std::setprecision(std::numeric_limits<multiFloat>::digits10 + 1);

    Eigen::Matrix<multiFloat ,4, 1> k1, k2, k3, k4, k5, k6;

    multiFloat t{};

    /*********************************************************/
    /*            ここから計算には関係ない                   */
    /*********************************************************/
/*
    std::string dir_name = "L1-" + L1.str();
    dir_name += "-L2-" + L2.str();
    dir_name += "-m1-" + m1.str();
    dir_name += "-m2-" + m2.str();
    dir_name += "-g-" + g.str();
    dir_name += "-theta1_init-" + theta1_init.str();
    dir_name += "-theta2_init-" + theta2_init.str();
    
    char c_dir_name[dir_name.size() + 1];
    std::strcpy(c_dir_name, dir_name.c_str());
    
    std::string str =  "mkdir -p " + dir_name;

    char com_dir[str.size() + 1];
    std::strcpy(com_dir, str.c_str());
    
    if(sim){
        //system
        if(system(com_dir)){
            std::cout << dir_name + "was not created." << std::endl;
        }
    }
*/
    // gnuplot

    /*
    FILE *gp = popen("gnuplot -persist","w");
    fprintf(gp, "set xr [-2.1:2.1]\n");
    fprintf(gp, "set yr [%f:%f]\n", -2.1, 2.1);
    fprintf(gp, "set size square\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "unset key\n");
    if(sim) fprintf(gp, "set term png\n");
    */


    /*********************************************************/
    /*            ここまで計算には関係ない                   */
    /*********************************************************/

    auto initPot = kineticEnergy<>(x) + potentialEnergy<>(x);

    //std::cout << initPot << std::endl;

    for(std::size_t i{}; t<t_limit; ++i){
        //std::cout << t << " " << x(0, 0) << " " << x(1, 0) << " " << x(2, 0) << " " << x(3, 0) << std::endl;
        //std::cout << t << " " << x(0, 0) << " " << x(1, 0) << " " << x(2, 0) << " " << x(3, 0) << " " << kineticEnergy<>(x) + potentialEnergy<>(x) << std::endl;
        //std::cout << t << " " << std::log10(dt) << " " << std::log10(std::abs(kineticEnergy<>(x) + potentialEnergy<>(x) - initPot)) << std::endl;
        //std::cout << t << " " << mp::log10(dt) << " " << mp::log10(mp::abs(kineticEnergy<>(x) + potentialEnergy<>(x) - initPot)) << std::endl;

    
        //描画
        if(i%INTV == 0){
            //std::cout << t << " " << kineticEnergy<>(x) << " " << potentialEnergy<>(x) << " " << kineticEnergy<>(x) + potentialEnergy<>(x) << std::endl;
            std::cout << t << " " << mp::log10(dt) << " " << mp::log10(mp::abs(kineticEnergy<>(x) + potentialEnergy<>(x) - initPot)) << std::endl;
            /*
            if(sim) fprintf(gp, "set output '%s/%06d.png'\n", c_dir_name, static_cast<int>(i/INTV));
            fprintf(gp, "plot '-' w lp lw 3 pt 7 ps 3\n");
            fprintf(gp, "0.0 0.0\n");
            fprintf(gp, "%f %f\n", L1 * std::sin(x(0, 0)), - L1 * std::cos(x(0,0)));
            fprintf(gp, "%f %f\n", L1 * std::sin(x(0, 0)) + L2 * std::sin(x(2,0)), - L1 * std::cos(x(0, 0)) - L2 * std::cos(x(2, 0)));
            fprintf(gp, "e\n");
            if(sim) fprintf(gp, "set output\n");
            fflush(gp);
            */
        }
    

        //RKF45法で常微分方程式を解く．
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

        //auto delta = static_cast<multiFloat>(0.86) * mp::pow(e_tol * dt / (2 * R), static_cast<multiFloat>(1)/4);
        auto delta = static_cast<multiFloat>("0.86") * static_cast<multiFloat>(mp::pow(static_cast<multiFloat>(e_tol/(2 * R)), static_cast<multiFloat>("0.25")));

        //この辺は適当にチューニング．最適な刻み幅制御は自分もよくわかっていない．
        
        if(delta < static_cast<multiFloat>("0.1")){
            dt = static_cast<multiFloat>("0.1") * dt;
            if(dt < t_min){
                dt = t_min;
            }
        }else if(delta > 4){
            dt = 4 * dt;
            if(dt > static_cast<multiFloat>("0.1")){
                dt = static_cast<multiFloat>("0.1");
            }
        }else{
            dt = delta * dt;
            if(dt > static_cast<multiFloat>("0.1")){
                dt = static_cast<multiFloat>("0.1");
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
