#pragma once

#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <string>
#include <random>

//#include <Eigen/Core>
//#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <boost/multiprecision/cpp_dec_float.hpp>

const multiFloat alpha("0.8");

const multiFloat order("-3.0");

namespace mino2357{

    template <typename T>
    constexpr T ratio(int a, int b){
        return static_cast<T>(a) / static_cast<T>(b);
    }

    /*
     * position x[i]
     *
     *
     */


    //R^NからR^Nへの関数．
    template <typename T = multiFloat>
    Eigen::Matrix<T, N, 1> func(const Eigen::Matrix<multiFloat, N, 1>& u){
        T x[planets][dim];
        T x_dot[planets][dim];
        for(int i=0; i<N/2; ++i){
            x[i/dim][i%dim]     = u(i, 0);
            x_dot[i/dim][i%dim] = u(N/2+i, 0);
        }

        T R[planets][planets] {};
        T F[planets][planets][dim] {};

        // distance
        for(int i=0; i<planets; ++i){
			for(int j=0; j<planets; ++j){
                for(int k=0; k<dim; ++k){
                    R[i][j] += (x[i][k] - x[j][k]) * (x[i][k] - x[j][k]);
                }
            }
        }

        for(int i=0; i<planets; ++i){
			for(int j=0; j<planets; ++j){
                R[i][j] = mp::sqrt(R[i][j]);
            }
        }

        Eigen::Matrix<multiFloat, N, 1> ans{};

        for(int i=0; i<N/2; ++i){
            ans(i, 0) = x_dot[i/dim][i%dim];
        }

		T m[planets];
        m[0] = 4.0;
        m[1] = 5.0;
        m[2] = 3.0;
/*
        m[3] = 1.0;
        m[4] = 1.0;
        m[5] = 1.0;
        m[6] = 1.0;
*/	
        for(int i=0; i<planets; ++i){
            for(int j=0; j<planets; ++j){
                if(i != j){
                    for(int k=0; k<dim; ++k){
                        F[i][j][k] += - m[j] * (x[i][k] - x[j][k])/mp::pow(R[i][j], 1.0 + order);
                        //F[i][j][k] += - (x[i][k] - x[j][k])/(R[i][j] * R[i][j]* R[i][j]);
                    }
                }
            }
        }

        T F_vector[planets][dim] {};

        for(int i=0; i<planets; ++i){
            for(int j=0; j<planets; ++j){
                for(int k=0; k<dim; ++k){
                    F_vector[i][k] += F[i][j][k];
                }
            }
        }
        
        for(int j=0; j<N/2;){
            for(int i=0; i<planets; ++i){
                for(int k=0; k<dim; ++k){
                    ans(N/2 + j, 0) = F_vector[i][k];
                    j++;
                }
            }
        }
        return ans;
    }
    
    /*************************************************************************/    
    /*************************************************************************/    
    /*************************************************************************/    

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

    template <typename T>
    class RKF45{
    public:
        T crt_h;
        T next_h;
        T A_Tol;
        T R_Tol;//今は使わない

        constexpr RKF45(T, T);

        inline constexpr void Integrate(T& ,T&, Eigen::Matrix<T, N, 1>&) noexcept;
    };

    template <typename T>
    inline constexpr void RKF45<T>::Integrate(T& t, T& dt, Eigen::Matrix<T, N, 1>& x) noexcept{
        crt_h = dt;

        Eigen::Matrix<T, N, 1> x5, x6, temp;
        T delta, R;

        R = 0;

        Eigen::Matrix<T, N, 1> k0, k1, k2, k3, k4, k5;

        ButcherRKF45<T> bf45;

        k0 = func<T>(x);
        k1 = func<T>(x + crt_h * bf45(1, 0) * k0);
        k2 = func<T>(x + crt_h * bf45(2, 0) * k0 + crt_h * bf45(2, 1) * k1);
        k3 = func<T>(x + crt_h * bf45(3, 0) * k0 + crt_h * bf45(3, 1) * k1 + crt_h * bf45(3, 2) * k2);
        k4 = func<T>(x + crt_h * bf45(4, 0) * k0 + crt_h * bf45(4, 1) * k1 + crt_h * bf45(4, 2) * k2 + crt_h * bf45(4, 3) * k3);
        k5 = func<T>(x + crt_h * bf45(5, 0) * k0 + crt_h * bf45(5, 1) * k1 + crt_h * bf45(5, 2) * k2 + crt_h * bf45(5, 3) * k3 + crt_h * bf45(5, 4) * k4);

        x5 = x + crt_h * (bf45.o5(0) * k0 + bf45.o5(1) * k1 + bf45.o5(2) * k2 + bf45.o5(3) * k3 + bf45.o5(4) * k4 + bf45.o5(5) * k5);
        x6 = x + crt_h * (bf45.o6(0) * k0 + bf45.o6(1) * k1 + bf45.o6(2) * k2 + bf45.o6(3) * k3 + bf45.o6(4) * k4 + bf45.o6(5) * k5);

        temp = x5 - x6;

        //temp  = (bf45.R(0) * k0 + bf45.R(2) * k2 + bf45.R(3) * k3 + bf45.R(4) * k4 + bf45.R(5) * k5) / crt_h;
        
        //delta = sqrt(temp(0, 0) * temp(0, 0) + temp(1, 0) * temp(1, 0) + temp(2, 0) * temp(2, 0) + temp(3, 0) * temp(3, 0));
        
        for(int i=0; i<N; ++i){
            R += temp(i, 0) * temp(i, 0);
        }
        R = R / crt_h;
        
        delta = mp::sqrt(R);

        if(delta > A_Tol){
            std::cerr << "Retry " << t << " " << dt << std::endl;
            dt = crt_h * mp::pow(alpha * A_Tol / delta, ratio<T>(1, 5));
            return;
        }

        x = x5;
        //x = x6;
        t += crt_h;

        next_h = crt_h * mp::pow(alpha * A_Tol / delta, ratio<T>(1, 5));

        dt = next_h;
    }

    template <typename T>
    constexpr RKF45<T>::RKF45(T at,T rt){
        A_Tol = at;
        R_Tol = rt;
    }
    
    /************************************************************************/
    
    template <typename T>
    class ButcherRKF78{
    public:
        T table[13][13];
        T order8[13];
        T order9[13];

        constexpr ButcherRKF78();
        
        constexpr T operator()(int i, int j){
            return table[i][j];
        };

        constexpr T o8(int i){
            return order8[i];
        };
        
        constexpr T o9(int i){
            return order9[i];
        };
    };

    template <typename T>
    constexpr ButcherRKF78<T>::ButcherRKF78(){
        T zero = 0;
        
        table[0][0] = zero;

        table[1][0] = ratio<T>(2, 27);
        table[0][1] = zero;
        
        table[2][0] = ratio<T>(1, 36);
        table[2][1] = ratio<T>(1, 12);
        table[2][2] = zero;
        
        table[3][0] = ratio<T>(1, 24);
        table[3][1] = zero;
        table[3][2] = ratio<T>(1, 8);
        table[3][3] = zero;
        
        table[4][0] = ratio<T>(5, 12);
        table[4][1] = zero;
        table[4][2] = ratio<T>(-25, 16);
        table[4][3] = ratio<T>(25, 16);
        table[4][4] = zero;
        
        table[5][0] = ratio<T>(1, 20);
        table[5][1] = zero;
        table[5][2] = zero;
        table[5][3] = ratio<T>(1, 4);
        table[5][4] = ratio<T>(1, 5);
        table[5][5] = zero;
        
        table[6][0] = ratio<T>(-25, 108);
        table[6][1] = zero;
        table[6][2] = zero;
        table[6][3] = ratio<T>(125, 108);
        table[6][4] = ratio<T>(-65, 27);
        table[6][5] = ratio<T>(125, 54);
        table[6][6] = zero;
        
        table[7][0] = ratio<T>(31, 300);
        table[7][1] = zero;
        table[7][2] = zero;
        table[7][3] = zero;
        table[7][4] = ratio<T>(61, 225);
        table[7][5] = ratio<T>(-2, 9);
        table[7][6] = ratio<T>(13, 900);
        table[7][7] = zero;
        
        table[8][0] = ratio<T>(2, 1);
        table[8][1] = zero;
        table[8][2] = zero;
        table[8][3] = ratio<T>(-53, 6);
        table[8][4] = ratio<T>(704, 45);
        table[8][5] = ratio<T>(-107, 9);
        table[8][6] = ratio<T>(67, 90);
        table[8][7] = ratio<T>(3, 1);
        table[8][8] = zero;
        
        table[9][0] = ratio<T>(-91, 108);
        table[9][1] = zero;
        table[9][2] = zero;
        table[9][3] = ratio<T>(23, 108);
        table[9][4] = ratio<T>(-976, 135);
        table[9][5] = ratio<T>(311, 54);
        table[9][6] = ratio<T>(-19, 60);
        table[9][7] = ratio<T>(17, 6);
        table[9][8] = ratio<T>(-1, 12);
        table[9][9] = zero;
        
        table[10][0]  = ratio<T>(2383, 4100);
        table[10][1]  = zero;
        table[10][2]  = zero;
        table[10][3]  = ratio<T>(-341, 164);
        table[10][4]  = ratio<T>(4496, 1025);
        table[10][5]  = ratio<T>(-301, 82);
        table[10][6]  = ratio<T>(2133, 4100);
        table[10][7]  = ratio<T>(45, 82);
        table[10][8]  = ratio<T>(45, 164);
        table[10][9]  = ratio<T>(18, 41);
        table[10][10] = zero;
        
        table[11][0]  = ratio<T>(3, 205);
        table[11][1]  = zero;
        table[11][2]  = zero;
        table[11][3]  = zero;
        table[11][4]  = zero;
        table[11][5]  = ratio<T>(-6, 41);
        table[11][6]  = ratio<T>(-3, 205);
        table[11][7]  = ratio<T>(-3, 41);
        table[11][8]  = ratio<T>(3, 41);
        table[11][9]  = ratio<T>(6, 41);
        table[11][10] = zero;
        table[11][11] = zero;
        
        table[12][0]  = ratio<T>(-1777, 4100);
        table[12][1]  = zero;
        table[12][2]  = zero;
        table[12][3]  = ratio<T>(-341, 164);
        table[12][4]  = ratio<T>(4496, 1025);
        table[12][5]  = ratio<T>(-289, 82);
        table[12][6]  = ratio<T>(2193, 4100);
        table[12][7]  = ratio<T>(51, 82);
        table[12][8]  = ratio<T>(33, 164);
        table[12][9]  = ratio<T>(12, 41);
        table[12][10] = zero;
        table[12][11] = ratio<T>(1, 1);
        table[12][12] = zero;

        order8[0]  = ratio<T>(41, 840);
        order8[1]  = zero;
        order8[2]  = zero;
        order8[3]  = zero;
        order8[4]  = zero;
        order8[5]  = ratio<T>(34, 105);
        order8[6]  = ratio<T>(9, 35);
        order8[7]  = ratio<T>(9, 35);
        order8[8]  = ratio<T>(9, 280);
        order8[9]  = ratio<T>(9, 280);
        order8[10] = ratio<T>(41, 840);
        order8[11] = zero;
        order8[12] = zero;
        
        order9[0]  = zero;
        order9[1]  = zero;
        order9[2]  = zero;
        order9[3]  = zero;
        order9[4]  = zero;
        order9[5]  = ratio<T>(34, 105);
        order9[6]  = ratio<T>(9, 35);
        order9[7]  = ratio<T>(9, 35);
        order9[8]  = ratio<T>(9, 280);
        order9[9]  = ratio<T>(9, 280);
        order9[10] = zero;
        order9[11] = ratio<T>(41, 840);
        order9[12] = ratio<T>(41, 840);
    }

    template <typename T>
    class RKF78{
    public:
        T crt_h;
        T next_h;
        T A_Tol;
        T R_Tol;//今は使わない

        constexpr RKF78(T, T);
        
        inline constexpr void Integrate(T& t, T& dt, Eigen::Matrix<T, N, 1>& x) noexcept;
    };

    template <typename T>
    constexpr RKF78<T>::RKF78(T at,T rt){
        A_Tol = at;
        R_Tol = rt;
    }

    template <typename T>
    inline constexpr void RKF78<T>::Integrate(T& t, T& dt, Eigen::Matrix<T, N, 1>& x) noexcept{
        crt_h = dt;

        Eigen::Matrix<T, N, 1> x9, x8, temp;
        T R, delta;
        R = 0;

        Eigen::Matrix<T, N, 1> k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12;

        ButcherRKF78<T> f;

        k0 = func<T>(x);
        k1 = func<T>(x + crt_h * f(1, 0) * k0);
        k2 = func<T>(x + crt_h * f(2, 0) * k0 + crt_h * f(2, 1) * k1);
        k3 = func<T>(x + crt_h * f(3, 0) * k0 + crt_h * f(3, 1) * k1 + crt_h * f(3, 2) * k2);
        k4 = func<T>(x + crt_h * f(4, 0) * k0 + crt_h * f(4, 1) * k1 + crt_h * f(4, 2) * k2 + crt_h * f(4, 3) * k3);
        k5 = func<T>(x + crt_h * f(5, 0) * k0 + crt_h * f(5, 1) * k1 + crt_h * f(5, 2) * k2 + crt_h * f(5, 3) * k3 + crt_h * f(5, 4) * k4);
        k6 = func<T>(x + crt_h * f(6, 0) * k0 + crt_h * f(6, 1) * k1 + crt_h * f(6, 2) * k2 + crt_h * f(6, 3) * k3 + crt_h * f(6, 4) * k4
                       + crt_h * f(6, 5) * k5);
        k7 = func<T>(x + crt_h * f(7, 0) * k0 + crt_h * f(7, 1) * k1 + crt_h * f(7, 2) * k2 + crt_h * f(7, 3) * k3 + crt_h * f(7, 4) * k4
                       + crt_h * f(7, 5) * k5 + crt_h * f(7, 6) * k6);
        k8 = func<T>(x + crt_h * f(8, 0) * k0 + crt_h * f(8, 1) * k1 + crt_h * f(8, 2) * k2 + crt_h * f(8, 3) * k3 + crt_h * f(8, 4) * k4
                       + crt_h * f(8, 5) * k5 + crt_h * f(8, 6) * k6 + crt_h * f(8, 7) * k7);
        k9 = func<T>(x + crt_h * f(9, 0) * k0 + crt_h * f(9, 1) * k1 + crt_h * f(9, 2) * k2 + crt_h * f(9, 3) * k3 + crt_h * f(9, 4) * k4
                       + crt_h * f(9, 5) * k5 + crt_h * f(9, 6) * k6 + crt_h * f(9, 7) * k7 + crt_h * f(9, 8) * k8);
        k10= func<T>(x + crt_h * f(10, 0) * k0 + crt_h * f(10, 1) * k1 + crt_h * f(10, 2) * k2 + crt_h * f(10, 3) * k3 + crt_h * f(10, 4) * k4
                       + crt_h * f(10, 5) * k5 + crt_h * f(10, 6) * k6 + crt_h * f(10, 7) * k7 + crt_h * f(10, 8) * k8 + crt_h * f(10, 9) * k9);
        k11= func<T>(x + crt_h * f(11, 0) * k0 + crt_h * f(11, 1) * k1 + crt_h * f(11, 2) * k2 + crt_h * f(11, 3) * k3 + crt_h * f(11, 4) * k4
                       + crt_h * f(11, 5) * k5 + crt_h * f(11, 6) * k6 + crt_h * f(11, 7) * k7 + crt_h * f(11, 8) * k8 + crt_h * f(11, 9) * k9
                       + crt_h * f(11, 10) * k10);
        k12= func<T>(x + crt_h * f(12, 0) * k0 + crt_h * f(12, 1) * k1 + crt_h * f(12, 2) * k2 + crt_h * f(12, 3) * k3 + crt_h * f(12, 4) * k4
                       + crt_h * f(12, 5) * k5 + crt_h * f(12, 6) * k6 + crt_h * f(12, 7) * k7 + crt_h * f(12, 8) * k8 + crt_h * f(12, 9) * k9
                       + crt_h * f(12, 10) * k10 + crt_h * f(12, 11) * k11);
        
        x8 = x + crt_h * (f.o8(0) * k0
                        + f.o8(1) * k1
                        + f.o8(2) * k2
                        + f.o8(3) * k3
                        + f.o8(4) * k4
                        + f.o8(5) * k5
                        + f.o8(6) * k6
                        + f.o8(7) * k7
                        + f.o8(8) * k8
                        + f.o8(9) * k9
                        + f.o8(10) * k10
                        + f.o8(11) * k11
                        + f.o8(12) * k12);
        
        x9 = x + crt_h * (f.o9(0) * k0
                        + f.o9(1) * k1
                        + f.o9(2) * k2
                        + f.o9(3) * k3
                        + f.o9(4) * k4
                        + f.o9(5) * k5
                        + f.o9(6) * k6
                        + f.o9(7) * k7
                        + f.o9(8) * k8
                        + f.o9(9) * k9
                        + f.o9(10) * k10
                        + f.o9(11) * k11
                        + f.o9(12) * k12);

        temp = x8 - x9;

        for(int i=0; i<N; ++i){
            R += temp(i, 0) * temp(i, 0);
        }
        R = R / crt_h;

        delta = mp::sqrt(R);
       
       
        if(delta > A_Tol){
            std::cerr << "Retry " << t << " " << dt << std::endl;
            dt = crt_h * mp::pow(alpha * A_Tol / delta, ratio<T>(1, 5));
            return;
        }if(10.0 * delta < A_Tol){
			dt = 10.0 * dt;
		}

        x = x9;

        t += crt_h;
        
        next_h = crt_h * mp::pow(alpha * A_Tol / delta, ratio<T>(1, 8));

		if(next_h > t_max){
			next_h = t_max;
		}

        dt = next_h;
    }


}
