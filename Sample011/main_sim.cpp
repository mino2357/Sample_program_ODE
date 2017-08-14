/*
 * RKF45法で二階の常微分方程式を解く．
 * テーマ：二体問題．
 * 線形代数ライブラリのEigenを使用．
 * Boostも使用．
 *
 * by みーくん．
 */

#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <string>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace mp = boost::multiprecision;

using multiFloat = mp::cpp_dec_float_100;

// パラメータ
const multiFloat e("0.9");

//時刻に関するパラメータ
multiFloat dt("1.0e-5");
const multiFloat t_limit("100.0");

const multiFloat RTol("10e-10");
const multiFloat ATol("10e-10");
const multiFloat t_min("10e-50");
const multiFloat t_max("0.1");

const multiFloat alpha("0.8");

//インターバル
constexpr int INTV = 1;

//movie
constexpr int sim = 0;

namespace mino2357{

    template <typename T>
    constexpr T ratio(int a, int b){
        return static_cast<T>(a) / static_cast<T>(b);
    }

    template <typename T>
    class ButcherRKF45{
    public:
        T table[6][6];
        T order5[6];
        T order6[6];
        T minus[6];

        constexpr ButcherRKF45();

        constexpr T operator()(int i, int j){
            return table[i][j];
        };

        constexpr T o5(int i){
            return order5[i];
        };
        
        constexpr T o6(int i){
            return order6[i];
        };

        constexpr T R(int i){
            return minus[i];
        }
    };

    template <typename T>
    constexpr ButcherRKF45<T>::ButcherRKF45(){
        table[0][0] =   ratio<T>(0, 1);
        table[1][0] =   ratio<T>(1, 4);
        table[2][0] =   ratio<T>(3, 32);      table[2][1] =   ratio<T>(9, 32);
        table[3][0] =   ratio<T>(1932, 2197); table[3][1] = - ratio<T>(7200, 2197);  table[3][2] =   ratio<T>(7296, 2197);
        table[4][0] =   ratio<T>(439, 216);   table[4][1] = - ratio<T>(8, 1);        table[4][2] =   ratio<T>(3680, 513);   table[4][3] = - ratio<T>(845, 4104);
        table[5][0] = - ratio<T>(8, 27);      table[5][1] =   ratio<T>(2, 1);        table[5][2] = - ratio<T>(3544, 2565);  table[5][3] =   ratio<T>(1859, 4104);  table[5][4] = - ratio<T>(11, 40);
        
        order6[0] = ratio<T>(16, 135);      order6[1] =   ratio<T>(0, 1);   order6[2] = ratio<T>(6656, 12825);
        order6[3] = ratio<T>(28561, 56430); order6[4] = - ratio<T>(9, 50);  order6[5] = ratio<T>(2, 55);
    
        order5[0] = ratio<T>(25, 216);      order5[1] =   ratio<T>(0, 1);   order5[2] = ratio<T>(1408, 2565);
        order5[3] = ratio<T>(2197, 4104);   order5[4] = - ratio<T>(1, 5);   order5[5] = ratio<T>(0, 1);

        minus[0] =   ratio<T>(1, 360);         minus[1] = ratio<T>(0, 1);      minus[2] = - ratio<T>(128, 4275);
        minus[3] = - ratio<T>(2197, 75240);    minus[4] = ratio<T>(1, 50);     minus[5] =   ratio<T>(2, 55);
    }

    //R^4からR^4への関数．
    template <typename T = multiFloat>
    Eigen::Matrix<T, 4, 1> func(const Eigen::Matrix<multiFloat, 4, 1>& u){
        T x       = u(0, 0);
        T x_dot   = u(1, 0);
        T y       = u(2, 0);
        T y_dot   = u(3, 0);

        T R = x * x + y * y;

        return Eigen::Matrix<T, 4, 1> {
            x_dot,
            - x / mp::sqrt(R * R * R),
            y_dot,
            - y / mp::sqrt(R * R * R)
        };
    }

    template <typename T>
    class RKF45{
    public:
        T crt_h;
        T next_h;
        T A_Tol;
        T R_Tol;

        constexpr RKF45(T, T);

        inline constexpr void Integrate(T& ,T&, Eigen::Matrix<T, 4, 1>&) noexcept;
    };

    template <typename T>
    inline constexpr void RKF45<T>::Integrate(T& t, T& dt, Eigen::Matrix<T, 4, 1>& x) noexcept{
        crt_h = dt;

        Eigen::Matrix<T, 4, 1> x5, x6, temp; //Eigen::Matrix<T, 4, 1> t;
        T delta; //T tDelta;

        Eigen::Matrix<T ,4, 1> k0, k1, k2, k3, k4, k5;

        ButcherRKF45<T> bf45;

        k0 = func<T>(x);
        k1 = func<T>(x + crt_h * bf45(1, 0) * k0);
        k2 = func<T>(x + crt_h * bf45(2, 0) * k0 + crt_h * bf45(2, 1) * k1);
        k3 = func<T>(x + crt_h * bf45(3, 0) * k0 + crt_h * bf45(3, 1) * k1 + crt_h * bf45(3, 2) * k2);
        k4 = func<T>(x + crt_h * bf45(4, 0) * k0 + crt_h * bf45(4, 1) * k1 + crt_h * bf45(4, 2) * k2 + crt_h * bf45(4, 3) * k3);
        k5 = func<T>(x + crt_h * bf45(5, 0) * k0 + crt_h * bf45(5, 1) * k1 + crt_h * bf45(5, 2) * k2 + crt_h * bf45(5, 3) * k3 + crt_h * bf45(5, 4) * k4);

        x5 = x + crt_h * (bf45.o5(0) * k0 + bf45.o5(1) * k1 + bf45.o5(2) * k2 + bf45.o5(3) * k3 + bf45.o5(4) * k4 + bf45.o5(5) * k5);
        //x6 = x + crt_h * (bf45.o6(0) * k0 + bf45.o6(1) * k1 + bf45.o6(2) * k2 + bf45.o6(3) * k3 + bf45.o6(4) * k4 + bf45.o6(5) * k5);

        temp  = (bf45.R(0) * k0 + bf45.R(2) * k2 + bf45.R(3) * k3 + bf45.R(4) * k4 + bf45.R(5) * k5);
        delta = sqrt(temp(0, 0) * temp(0, 0) + temp(1, 0) * temp(1, 0) + temp(2, 0) * temp(2, 0) + temp(3, 0) * temp(3, 0));

        //tt = x5 - x6;

        //tDelta = sqrt(tt(0, 0) * tt(0, 0) + tt(1, 0) * tt(1, 0) + tt(2, 0) * tt(2, 0) + tt(3, 0) * tt(3, 0));

/*
        //RK4
        k1 = func<T>(x);
        k2 = func<T>(x + crt_h * ratio<T>(1, 2) * k1);
        k3 = func<T>(x + crt_h * ratio<T>(1, 2) * k2);
        k4 = func<T>(x + crt_h * ratio<T>(1, 1) * k3);

        x = x + ratio<T>(1, 6) * crt_h * (k1 + 2 * k2 + 2 * k3 + k4);
*/

        
        if(delta > A_Tol){
            dt = crt_h * mp::pow(alpha * A_Tol / delta, ratio<T>(1, 5));
            return;
        }

        x = x5;
        t += crt_h;
        
        std::cout << t << " " << dt << " " << delta << " " << mp::pow(alpha * A_Tol / delta, ratio<T>(1, 5)) << std::endl;

        next_h = crt_h * mp::pow(alpha * A_Tol / delta, ratio<T>(1, 5));

        dt = next_h;
    }

    template <typename T>
    constexpr RKF45<T>::RKF45(T at,T rt){
        A_Tol = at;
        R_Tol = rt;
    }
}


int main(){
    Eigen::Matrix<multiFloat, 4, 1> x(1 - e, 0, 0, mp::sqrt((1 + e) / (1 - e)));
    Eigen::Matrix<multiFloat, 4, 1> x4, x5;
    //std::cout << std::fixed << std::setprecision(std::numeric_limits<multiFloat>::digits10 + 1);
    std::cout << std::fixed << std::setprecision(35);

    multiFloat t{};

    mino2357::RKF45<multiFloat> rkf45(ATol, RTol);

    for(int i=0; t<t_limit; i++){
        //std::cout << t << " " << dt << " " << x(0,0) << " " << x(1,0) << " " << x(2,0) << " " << x(3,0) << std::endl;
        rkf45.Integrate(t, dt, x);
    }
    
    return 0;
/*
    for(std::size_t i{}; t<t_limit; ++i){
        
        if(i%INTV == 0){
        }
        

        //RK4法で常微分方程式を解く．
        k1 = func<>(x);
        k2 = func<>(x + dt / 4 * k1);
        k3 = func<>(x + dt / 32 * (3 * k1 + 9 * k2));
        k4 = func<>(x + dt / 2197 * (1932 * k1 - 7200 * k2 + 7296 * k3));
        k5 = func<>(x + dt * (static_cast<multiFloat>(439)/216 * k1 - 8 * k2 + static_cast<multiFloat>(3680)/513 * k3 - static_cast<multiFloat>(845)/4104 * k4));
        k6 = func<>(x + dt * (- static_cast<multiFloat>(8)/27 * k1 + 2 * k2 - static_cast<multiFloat>(3544)/2565 * k3 + static_cast<multiFloat>(1859)/4104 * k4 - static_cast<multiFloat>(11)/40 * k5));
        
        x4 = x + dt * (static_cast<multiFloat>(25)/216 * k1 + static_cast<multiFloat>(1408)/2565 * k3 + static_cast<multiFloat>(2197)/4104 * k4 - static_cast<multiFloat>(1)/5 * k5);
        x5 = x + dt * (static_cast<multiFloat>(16)/135 * k1 + static_cast<multiFloat>(6656)/12825 * k3 + static_cast<multiFloat>(28561)/56430 * k4 - static_cast<multiFloat>(9)/50 * k5 + static_cast<multiFloat>(2)/55 * k6);
       
        auto temp  = (static_cast<multiFloat>(25)/216 * k1 + static_cast<multiFloat>(1408)/2565 * k3 + static_cast<multiFloat>(2197)/4104 * k4 - static_cast<multiFloat>(1)/5 * k5) - (static_cast<multiFloat>(16)/135 * k1 + static_cast<multiFloat>(6656)/12825 * k3 + static_cast<multiFloat>(28561)/56430 * k4 - static_cast<multiFloat>(9)/50 * k5 + static_cast<multiFloat>(2)/55 * k6);
        multiFloat R = mp::sqrt((temp(0, 0) * temp(0, 0) + temp(1, 0) * temp(1, 0) + temp(2, 0) * temp(2, 0) + temp(3, 0) * temp(3, 0)));

        multiFloat x4_abs = mp::sqrt(static_cast<multiFloat>(x4(0, 0) * x4(0,0) + x4(1, 0) * x4(1,0) + x4(2, 0) * x4(2,0) + x4(3, 0) * x4(3,0)));

        multiFloat E = static_cast<multiFloat>(1)/static_cast<multiFloat>(2) * (ATol + x4_abs * RTol);



        if(mp::pow(E/R, static_cast<multiFloat>(1)/5) < 1){
            dt = dt * 0.6 * mp::pow(E/R, static_cast<multiFloat>(1)/5);
            continue;
        }

        x  = x4;
        //std::cout << t << " " << mp::log10(dt) << " " << x(0, 0) << " " << x(1, 0) << " " << x(2, 0) << " " << x(3, 0) << std::endl;
        std::cout << t << " " << mp::log10(dt) << " " << R << " " << mp::pow(E/R, static_cast<multiFloat>(1)/5) << " " << x(0, 0) << " " << x(2, 0) << " " <<mp::sqrt(static_cast<multiFloat>(x4(1, 0) * x4(1,0) + x4(3, 0) * x4(3,0))) << std::endl;

        t += dt;

        //std::cout << x4 << std::endl;
        
        dt = 10.0 * dt;
        if(dt > t_max){
            dt = t_max;
        }
    }
    */  
}
