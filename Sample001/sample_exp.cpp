/*
 * dy/dt = y
 * y(0) = y_init
 *
 * をオイラー法で解くよ．
 *
 * みーくん氏
 */


#include <iostream>

//パラメータ
const double y_init  = 1.0;
const double dt      = 0.001; //時間の刻み幅だよ．
const double T_limit = 10;

double func(double y){
    return y;
}

int main(){
    //t=0からはじめるよ．
    double t = 0.;
    //常微分方程式の初期条件を設定するよ．
    double y = y_init;

    //はじめの値を表示するよ.
    std::cout << t << " " << y_init << std::endl;

    //漸化式を解くよ．T_limit秒まで解くよ．
    for(int i=1; t<T_limit; i++){
        //dt秒後のyの値を求めるよ．
        y = y + dt * func(y);
        //時間tを進めるよ．
        t = i * dt;
        //t秒後のときのyの値を表示するよ．
        std::cout << t << " " << y << std::endl;
    }
}
