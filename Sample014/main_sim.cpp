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
namespace mp = boost::multiprecision;
using multiFloat = mp::cpp_dec_float_100;

int main(){

    std::mt19937 mt;
    std::random_device rnd;
    mt.seed(rnd());
    std::uniform_real_distribution<> rand10(0, 10);

    Eigen::Matrix<multiFloat, N, 1> x{};
    // initial values
    for(int i=0; i<N/2; ++i){
        x(i, 0) = rand10(mt);
    }
    
    //std::cout << std::fixed << std::setprecision(std::numeric_limits<multiFloat>::digits10 + 1);
    std::cout << std::fixed << std::setprecision(20);
    std::cerr << std::fixed << std::setprecision(20);

    multiFloat t{};

    mino2357::RKF78<multiFloat> rkf78(ATol, RTol);

    // gnuplot
    FILE *gp = popen("gnuplot -persist","w");
    fprintf(gp, "set xr [-10.:20.]\n");
    fprintf(gp, "set yr [%f:%f]\n", -10., 20.);
    fprintf(gp, "set size square\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "unset key\n");
    //fprintf(gp, "set term png\n");


    for(int i=0; t<t_limit; i++){
        if(i%INTV == 0){
            //std::cerr << t << " " << mp::log10(dt) << std::endl;
            
            fprintf(gp, "plot '-' pt 7 ps 2 lt rgb 'blue'\n");
            for(int j=0; j<N/2; j+=2){
                fprintf(gp, "%f %f\n", static_cast<double>(x(j,0)), static_cast<double>(x(j+1,0)));
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        rkf78.Integrate(t, dt, x);
    }
    pclose(gp);
}
