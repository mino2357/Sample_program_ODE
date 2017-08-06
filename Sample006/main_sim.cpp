/*
 * RK4法で二階の常微分方程式を解く．
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
constexpr double L2          = 0.6;
//初期角
constexpr double theta1_init = 2.5;
constexpr double theta2_init = 2.0;
//初角速度
constexpr double omega1_init = 0.;
constexpr double omega2_init = 0.;

//時刻に関するパラメータ
constexpr double dt          =     0.001;
constexpr double t_limit     = 20000.0;

//インターバル
constexpr int INTV = 40;

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
	std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);

	Eigen::Matrix<double ,4, 1> k1, k2, k3, k4;

    double t{};

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
	//system
	if(system(com_dir)){
		std::cout << dir_name + "was not created." << std::endl;
	}

	// gnuplot
	FILE *gp = popen("gnuplot -persist","w");
    fprintf(gp, "set xr [-2.1:2.1]\n");
    fprintf(gp, "set yr [%f:%f]\n", -2.1, 2.1);
    fprintf(gp, "set size square\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set term png\n");

    for(std::size_t i{}; t<t_limit; ++i){
        //std::cout << t << " " << x(0, 0) << " " << x(1, 0) << " " << x(2, 0) << " " << x(3, 0) << " " << energy<>(x) << std::endl;

		//描画
        if(i%INTV == 0){
        	std::cout << t << " " << kineticEnergy<>(x) << " " << potentialEnergy<>(x) << " " << kineticEnergy<>(x) + potentialEnergy<>(x) << std::endl;
            fprintf(gp, "set output '%s/%06d.png'\n", c_dir_name, static_cast<int>(i/INTV));
            fprintf(gp, "plot '-' w lp lw 3 pt 7 ps 3\n");
            fprintf(gp, "0.0 0.0\n");
            fprintf(gp, "%f %f\n", L1 * std::sin(x(0, 0)), - L1 * std::cos(x(0,0)));
            fprintf(gp, "%f %f\n", L1 * std::sin(x(0, 0)) + L2 * std::sin(x(2,0)), - L1 * std::cos(x(0, 0)) - L2 * std::cos(x(2, 0)));
            fprintf(gp, "e\n");
            fprintf(gp, "set output\n");
            fflush(gp);
        }

        //RK4法で常微分方程式を解く．
        k1 = func<>(x);
        k2 = func<>(x + dt / 2. * k1);
        k3 = func<>(x + dt / 2. * k2);
        k4 = func<>(x + dt * k3);
        x = x + dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
        t = i * dt;
    }
    pclose(gp);
}
