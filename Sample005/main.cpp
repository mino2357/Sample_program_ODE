/*
 * RK4法で二階の常微分方程式を解く．
 * テーマ：単振り子．
 * 線形代数ライブラリのEigenを使用．
 *
 * by みーくん．
 */

#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

// パラメータ
constexpr double g      = 1.;//9.8
constexpr double l      = 1.;
//初期位置
constexpr double theta_init = 1.;
//初速度
constexpr double omega_init = 0.;

//時刻に関するパラメータ
constexpr double dt      = 0.001;
constexpr double t_limit = 20.0;

//R^2からR^2への関数．
Eigen::Matrix<double, 2, 1> func(const Eigen::Matrix<double, 2, 1>& x){
	Eigen::Matrix<double, 2, 1> a;
	a(0, 0) = x(1, 0);
	a(1, 0) = - g / l * std::sin(x(0, 0));
	return a;
}

int main(){
    Eigen::Matrix<double, 2, 1> x(theta_init, omega_init);
	std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);

	Eigen::Matrix<double ,2, 1> k1;
	Eigen::Matrix<double ,2, 1> k2;
	Eigen::Matrix<double ,2, 1> k3;
	Eigen::Matrix<double ,2, 1> k4;

    double t {};

    for(std::size_t i {}; t<t_limit; ++i){
        std::cout << t << " " << x(0, 0) << " " << x(1, 0) << std::endl;
        //RK4法で常微分方程式を解く．
        k1 = func(x);
        k2 = func(x + dt / 2. * k1);
        k3 = func(x + dt / 2. * k2);
        k4 = func(x + dt * k3);
        x = x + dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
        t = i * dt;
    }
}
