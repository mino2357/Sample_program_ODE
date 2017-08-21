/*
 * RKF45法で二階の常微分方程式を解く．
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
#include <random>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <boost/multiprecision/cpp_dec_float.hpp>

int main(){

    std::mt19937 mt;
    std::random_device rnd;
    mt.seed(rnd());
    std::uniform_real_distribution<> rand10(-1, 1);

    Eigen::Matrix<multiFloat, N, 1> x{};
    // initial values
    for(int i=0; i<N/2; ++i){
        x(i, 0) = mp::pow(i%2 ==0? 2: 8, static_cast<multiFloat>(0.01*i));//rand10(mt);
        //x(N/2 + i, 0) = rand10(mt);
    }
	//x(0, 0) = 0;
	//x(1, 0) = 0;
    
    //std::cout << std::fixed << std::setprecision(std::numeric_limits<multiFloat>::digits10 + 1);
    std::cout << std::fixed << std::setprecision(12);
    std::cerr << std::fixed << std::setprecision(12);

    multiFloat t{};

    mino2357::RKF78<multiFloat> rkf78(ATol, RTol);

    // gnuplot
    FILE *gp = popen("gnuplot -persist","w");
    fprintf(gp, "set xr [-2.:2.]\n");
    fprintf(gp, "set yr [-2.:2.]\n");
    //fprintf(gp, "set zr [-2.:2.]\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "unset key\n");
    //fprintf(gp, "set term png\n");


    for(int i=0; t<t_limit; i++){
        if(i%INTV == 0){
            std::cerr << t << " " << mp::log10(dt) << std::endl;
            
			for(int j=0; j<N; ++j){
                std::cout << x(j,0) << " ";
            }
            std::cout << std::endl;
            
            fprintf(gp, "plot '-' pt 7 ps 1 lt rgb 'blue'\n");
            //fprintf(gp, "splot '-' pt 7 ps 1 lt rgb 'blue'\n");
            for(int j=0; j<N/2; j+=dim){
                fprintf(gp, "%f %f\n", static_cast<double>(x(j,0)), static_cast<double>(x(j+1,0)));
                //fprintf(gp, "%f %f %f\n", static_cast<double>(x(j,0)), static_cast<double>(x(j+1,0)), static_cast<double>(x(j+2,0)));
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        rkf78.Integrate(t, dt, x);
    }
    pclose(gp);
}
