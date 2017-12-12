/*
 * RK4法で二階の常微分方程式を解く．
 * テーマ：二重振り子．
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
constexpr double g       = 1.;//9.8
constexpr double m1      = 1.;
constexpr double m2      = 1.;
constexpr double L1      = 1.;
constexpr double L2      = 1.;
//初期角
constexpr double theta1_init = 0.1;
constexpr double theta2_init = 0.1;
//初角速度
constexpr double omega1_init = 0.;
constexpr double omega2_init = 0.;

//時刻に関するパラメータ
constexpr double dt      =  0.001;
constexpr double t_limit = 20.0;

//R^2からR^2への関数．
Eigen::Matrix<double, 4, 1> func(const Eigen::Matrix<double, 4, 1>& x){
	double theta1 = x(0, 0);
	double eta1   = x(1, 0);
	double theta2 = x(2, 0);
	double eta2   = x(3, 0);

	return Eigen::Matrix<double, 4, 1> {
		eta1,
		(- m1 * g * std::sin(theta1) - m2 * (g * std::sin(theta1) + L2 * eta2 * eta2 * std::sin(theta1 - theta2) + (L1 * eta1 * eta1 * std::sin(theta1 - theta2) - g * std::sin(theta1 - theta2)) * std::cos(theta1 - theta2))) / (L1 * (m1 + m2 * (std::sin(theta1 - theta2) * (std::sin(theta1 - theta2))))),
		eta2,
		((m1 + m2) * (L1 * eta1 * eta1 * std::sin(theta1 - theta2) - g * std::sin(theta2) + g * std::sin(theta1) * std::cos(theta1 - theta2)) + m2 * L2 * eta2 * eta2 * std::cos(theta1 - theta2) * std::sin(theta1 - theta2)) / (L2 * (m1 + m2 * std::sin(theta1 - theta2) * std::sin(theta1 - theta2)))
	};
}

int main(){
    Eigen::Matrix<double, 4, 1> x(theta1_init, omega1_init, theta2_init, omega2_init);
	std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);

	Eigen::Matrix<double ,4, 1> k1;
	Eigen::Matrix<double ,4, 1> k2;
	Eigen::Matrix<double ,4, 1> k3;
	Eigen::Matrix<double ,4, 1> k4;

    double t{};

    for(std::size_t i{}; t<t_limit; ++i){
        std::cout << t << " " << x(0, 0) << " " << x(1, 0) << " " << x(2, 0) << " " << x(3, 0) << std::endl;
        //RK4法で常微分方程式を解く．
        k1 = func(x);
        k2 = func(x + dt / 2. * k1);
        k3 = func(x + dt / 2. * k2);
        k4 = func(x + dt * k3);
        x = x + dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
        t = i * dt;
    }
}
