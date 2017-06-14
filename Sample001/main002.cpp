/*
 * calculationg ODE by Euler method.
 * dx/dt = f(x).
 * x is real number.
 *
 * @Takaaki MINOMO
 */


#include <iostream>
#include <iomanip>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <exception>

namespace mp = boost::multiprecision;
using f100 = mp::cpp_dec_float_100;

namespace mino2357{
    template <typename T = double>
    //constexpr T ratio(const T &a, const T &b){
    T ratio(const T &a, const T &b){
        if(b == 0){
            throw std::range_error("divided by zero!.");
        }
        return static_cast<T>(a / b);
    }       
}

/* parameter */
template <typename T = double>
const T alpha = mino2357::ratio<T>(1, 1);

template <typename T = double>
const T beta = mino2357::ratio<T>(1, 1);

template <typename T = double>
const T t_limit = mino2357::ratio<T>(10, 1);

template <typename T = double>
const T dt = mino2357::ratio<T>(1, 1000);

/* initial value */
template <typename T = double>
const T x_init = mino2357::ratio<T>(1, 0); // divided by zero!!!!

/* f(x) */
template<typename T = double>
constexpr T func(T x){
    return alpha<T> * (1. - x / beta<T>) * x;
}


int main(){
    std::cout << std::fixed << std::setprecision(std::numeric_limits<f100>::digits10 + 1);
    
    try{
        x_init<f100>;
        t_limit<f100>;
        dt<f100>;
        func(x_init<f100>);
    }catch(std::range_error& exception){
        std::cout << "catch!!!!!!!!!!!" << std::endl;
    }

    auto x = x_init<f100>;
    auto t = static_cast<f100>(0.0);

    // initial value
    std::cout << t << " " << x << std::endl;

    // time loop
    for(auto i = 0; t<t_limit<f100>; i++){
        t = i*dt<f100>;
        // Euler method 
        x = x + dt<f100> * func(x);
        std::cout << t << " " << x << std::endl;
    }
}
