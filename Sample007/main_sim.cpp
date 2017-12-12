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
// パラメータ
constexpr double g           = 1.0;//9.8
constexpr double m1          = 1.0;
constexpr double m2          = 1.0;
constexpr double L1          = 1.0;
constexpr double L2          = 1.0;
//初期角
constexpr double theta1_init = 3.0;
constexpr double theta2_init = 3.0;
//初角速度
constexpr double omega1_init = 0.;
constexpr double omega2_init = 0.;

//時刻に関するパラメータ
double dt                    = 10e-3;
constexpr double t_limit     = 10000.0;

constexpr double e_tol = 10e-10;
constexpr double t_min = 10e-6;

//インターバル
constexpr int INTV = 100;

//movie
constexpr int sim = 0;

//R^4からR^4への関数．
template <typename T = double>
Eigen::Matrix<T, 4, 1> func(const Eigen::Matrix<double, 4, 1>& x){
    T theta1 = x(0, 0);
    T eta1   = x(1, 0);
    T theta2 = x(2, 0);
    T eta2   = x(3, 0);

    return Eigen::Matrix<T, 4, 1> {
        eta1,
        (- m1 * g * std::sin(theta1) - m2 * (g * std::sin(theta1) + L2 * eta2 * eta2 * std::sin(theta1 - theta2) + (L1 * eta1 * eta1 * std::sin(theta1 - theta2) - g * std::sin(theta2)) * std::cos(theta1 - theta2))) / (L1 * (m1 + m2 * (std::sin(theta1 - theta2) * (std::sin(theta1 - theta2))))),
        eta2,
        ((m1 + m2) * (L1 * eta1 * eta1 * std::sin(theta1 - theta2) - g * std::sin(theta2) + g * std::sin(theta1) * std::cos(theta1 - theta2)) + m2 * L2 * eta2 * eta2 * std::cos(theta1 - theta2) * std::sin(theta1 - theta2)) / (L2 * (m1 + m2 * std::sin(theta1 - theta2) * std::sin(theta1 - theta2)))
    };
}

template <typename T = double>
T potentialEnergy(Eigen::Matrix<T, 4, 1>& x){
    T theta1 = x(0, 0);
    T eta1   = x(1, 0);
    T theta2 = x(2, 0);
    T eta2   = x(3, 0);

    return 0.5 * m1 * L1 * L1 * eta1 * eta1 + 0.5 * m2 * (L1 * L1 * eta1 * eta1 + L2 * L2 * eta2 * eta2 + 2. * L1 * L2 * eta1 * eta2 * std::cos(theta1 - theta2));
}

template <typename T = double>
T kineticEnergy(Eigen::Matrix<T, 4, 1>& x){
    T theta1 = x(0, 0);
    T theta2 = x(2, 0);

    return - m1 * g * L1 * std::cos(theta1) - m2 * g * (L1 * std::cos(theta1) + L2 * std::cos(theta2));
}

int main(){
    Eigen::Matrix<double, 4, 1> x(theta1_init, omega1_init, theta2_init, omega2_init);
    Eigen::Matrix<double, 4, 1> x4, x5;
    std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    Eigen::Matrix<double ,4, 1> k1, k2, k3, k4, k5, k6;

    double t{};

    /*********************************************************/
    /*            ここから計算には関係ない                   */
    /*********************************************************/

    std::string dir_name = "L1-" + std::to_string(L1);
    dir_name += "-L2-" + std::to_string(L2);
    dir_name += "-m1-" + std::to_string(m1);
    dir_name += "-m2-" + std::to_string(m2);
    dir_name += "-g-" + std::to_string(g);
    dir_name += "-theta1_init-" + std::to_string(theta1_init);
    dir_name += "-theta2_init-" + std::to_string(theta2_init);

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

    // gnuplot
    FILE *gp = popen("gnuplot -persist","w");
    fprintf(gp, "set xr [-2.1:2.1]\n");
    fprintf(gp, "set yr [%f:%f]\n", -2.1, 2.1);
    fprintf(gp, "set size square\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "unset key\n");
    if(sim) fprintf(gp, "set term png\n");
    
    /*********************************************************/
    /*            ここまで計算には関係ない                   */
    /*********************************************************/

    auto initPot = kineticEnergy<>(x) + potentialEnergy<>(x);

    for(std::size_t i{}; t<t_limit; ++i){
        //std::cout << t << " " << x(0, 0) << " " << x(1, 0) << " " << x(2, 0) << " " << x(3, 0) << " " << kineticEnergy<>(x) + potentialEnergy<>(x) << std::endl;
        //std::cout << t << " " << std::log10(dt) << " " << std::log10(std::abs(kineticEnergy<>(x) + potentialEnergy<>(x) - initPot)) << std::endl;

        //描画
        if(i%INTV == 0){
            //std::cout << t << " " << kineticEnergy<>(x) << " " << potentialEnergy<>(x) << " " << kineticEnergy<>(x) + potentialEnergy<>(x) << std::endl;
            std::cout << t << " " << std::log10(dt) << " " << std::log10(std::abs(kineticEnergy<>(x) + potentialEnergy<>(x) - initPot)) << std::endl;
            if(sim) fprintf(gp, "set output '%s/%06d.png'\n", c_dir_name, static_cast<int>(i/INTV));
            fprintf(gp, "plot '-' w lp lw 3 pt 7 ps 3\n");
            fprintf(gp, "0.0 0.0\n");
            fprintf(gp, "%f %f\n", L1 * std::sin(x(0, 0)), - L1 * std::cos(x(0,0)));
            fprintf(gp, "%f %f\n", L1 * std::sin(x(0, 0)) + L2 * std::sin(x(2,0)), - L1 * std::cos(x(0, 0)) - L2 * std::cos(x(2, 0)));
            fprintf(gp, "e\n");
            if(sim) fprintf(gp, "set output\n");
            fflush(gp);
        }

        //RKF45法で常微分方程式を解く．
        k1 = func<>(x);
        k2 = func<>(x + dt / 4. * k1);
        k3 = func<>(x + dt / 32. * (3. * k1 + 9. * k2));
        k4 = func<>(x + dt / 2197. * (1932. * k1 - 7200. * k2 + 7296. * k3));
        k5 = func<>(x + dt * (439./216. * k1 - 8 * k2 + 3680./513. * k3 - 845./4104 * k4));
        k6 = func<>(x + dt * (- 8./27. * k1 + 2. * k2 - 3544./2565. * k3 + 1859./4104. * k4 - 11./40. * k5));
        
        x4 = x + dt * (25./216. * k1 + 1408./2565. * k3 + 2197./4104. * k4 - 1./5. * k5);
        x5 = x + dt * (16./135. * k1 + 6656./12825. * k3 + 28561./56430. * k4 - 9./50. * k5 + 2./55. * k6);
       
        auto temp  = (25./216. * k1 + 1408./2565. * k3 + 2197./4104. * k4 - 1./5. * k5) - (16./135. * k1 + 6656./12825. * k3 + 28561./56430. * k4 - 9./50. * k5 + 2./55. * k6);
        auto R = std::sqrt(temp(0, 0) * temp(0, 0) + temp(1, 0) * temp(1, 0) + temp(2, 0) * temp(2, 0) + temp(3, 0) * temp(3, 0));

        //std::cout << t << " " << dt << " " << R << " " << std::endl;
        
        auto delta = 0.86 * std::pow(e_tol * dt / R, 0.25);

        //この辺は適当にチューニング．最適な刻み幅制御は自分もよくわかっていない．経験的に決定する．
        if(delta < 0.1){
            dt = 0.1 * dt;
            if(dt < t_min){
                dt = t_min;
            }
        }else if(delta > 4){
            dt = 4. * dt;
            if(dt > 0.1){
                dt = 0.1;
            }
        }else{
            dt = delta * dt;
            if(dt > 0.1){
                dt = 0.1;
            }
            if(dt < t_min){
                dt = t_min;
            }
        }

        x  = x4;
        t += dt;
    }
    pclose(gp);
}

