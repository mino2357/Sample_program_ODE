/*
 * Euler法で二階の常微分方程式を解く．
 * テーマ：単振動．調和振動子．
 *
 * by みーくん．
 */

#include <iostream>
#include "myvector.hpp"

// パラメータ
constexpr double k      = 1.;
constexpr double m      = 1.;
constexpr double x_init = 1.;
constexpr double v_init = 0.;

//時刻に関するパラメータ
constexpr double dt      = 0.001;
constexpr double t_limit = 20.0;

//R^2からR^2への関数．
mino2357::vector<> func(mino2357::vector<>& x){
    return mino2357::vector<>{x.getComponentY(), - k / m * x.getComponentX()};
}

int main(){
    mino2357::vector<> x(x_init, v_init);

    double t {};

    for(size_t i {}; t<t_limit; ++i){
        std::cout << t << " " << x.getComponentX() << " " << x.getComponentY() << std::endl;
        //Euler法で常微分方程式を解く．       
        x = x + dt * func(x);
        t = i * dt;
    }
}
