/*
 * RKF法で二階の常微分方程式を解く．
 * テーマ：多体問題．
 * 線形代数ライブラリのEigenを使用．
 * Boostも使用．
 *
 * by みーくん．
 *
 * memo
 * 1m18s
 * 
 * 
 */


#define EIGEN_NO_DEBUG

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

int main(){

    Eigen::Matrix<multiFloat, N, 1> x{};
    // initial values
	x(0, 0) = 0.0;
	x(1, 0) = 0.0;

	x(2, 0) = 3.0;
	x(3, 0) = 0.0;

	x(4, 0) = 3.0;
	x(5, 0) = 4.0;
    
	//std::cout << std::fixed << std::setprecision(std::numeric_limits<multiFloat>::digits10 + 1);
    std::cout << std::fixed << std::setprecision(12);
    std::cerr << std::fixed << std::setprecision(12);

    multiFloat t{};

    mino2357::RKF78<multiFloat> rkf78(ATol, RTol);

    // gnuplot
    FILE *gp = popen("gnuplot -persist","w");
    fprintf(gp, "set term png size 1280, 720\n");
    fprintf(gp, "set xr [-4.:8.]\n");
    fprintf(gp, "set yr [-4.:8.]\n");
    //fprintf(gp, "set zr [-2.:2.]\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set tics font 'Times New Roman,20'\n");
    //fprintf(gp, "set term png\n");

	int time{};

    for(int i=0; t<t_limit; i++){
			
		if(i%INTV == 0){
            std::cerr << t << std::endl;
            //std::cout << x(0,0) << " " << x(1,0) << " " << x(2,0) << " " << x(3,0) << " "<< x(4,0) << " " << x(5,0) << std::endl;
			time++;
            
            fprintf(gp, "set title 't = %lf      log_{10}dt = %lf'\n", static_cast<double>(t),std::log10(static_cast<double>(dt)));
			fprintf(gp, "plot '-' with points pt 6 ps variable lc rgb var\n");
            fprintf(gp, "%f %f %d 0x0000FF\n", static_cast<double>(x(0,0)), static_cast<double>(x(1,0)), 4);
            fprintf(gp, "%f %f %d 0x00FF00\n", static_cast<double>(x(2,0)), static_cast<double>(x(3,0)), 5);
            fprintf(gp, "%f %f %d 0xFF0000\n", static_cast<double>(x(4,0)), static_cast<double>(x(5,0)), 3);
            fprintf(gp, "e\n");
            fprintf(gp, "set output '%06d.png'\n", time);
            fflush(gp);
       	}
        rkf78.Integrate(t, dt, x);
    }
    pclose(gp);
}
