/*
 * dy/dt = alpha * (1 - y / beta) y
 * y(0) = y_init
 *
 * を古典的ルンゲ・クッタ法で解く．
 *
 * みーくん氏
 */


#include <iostream>
#include <iomanip>
#include <cmath>

//パラメータ
const double y_init  =  0.1;
const double dt      =  0.001; //時間の刻み幅．
const double T_limit = 10.0;
const double alpha   =  1.0;
const double beta    =  1.0;

double func(double y){
    return alpha * (1.0 - y / beta) * y;
}

int main(){
    std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        
     //t=0からはじめる．
    double t = 0.0;
    //常微分方程式の初期条件を設定．
    double y = y_init;

    //RK法で使うkたちを宣言して初期化．
    double k1, k2, k3, k4;
    k1 = k2 = k3 = k4 = 0.0;
    
    double C2 = 0.0;

    //はじめの値を表示. 「時刻t シミュレーション結果 相対誤差の対数」の順で出力．
    C2 = 1.0 - beta / y_init;
    std::cout << t << " " << y_init << " " << std::log10(std::abs(y - beta / (1.0 - C2 * std::exp(- alpha * t))) / std::abs(beta / (1.0 - C2 * std::exp(- alpha * t)))) << std::endl;

    //漸化式を解くよ．T_limit秒まで解く．
    for(int i=1; t<T_limit; i++){
        //k_i(i=1,2,3,4)の値を求める．
        k1 = func(y);
        k2 = func(y + dt * k1 / 2.0);
        k3 = func(y + dt * k2 / 2.0);
        k4 = func(y + dt * k3);
        //dt秒後のyの値を求める．
        y = y + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        //時間tをdt秒進める．
        t = i * dt;
        //t秒後のときのyの値を表示する．
        std::cout << t << " " << y << " "<< std::log10(std::abs(y - beta / (1.0 - C2 * std::exp(- alpha * t))) / std::abs(beta / (1.0 - C2 * std::exp(- alpha * t)))) << std::endl;
    }
}
