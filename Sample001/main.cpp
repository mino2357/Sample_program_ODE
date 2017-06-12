/*
 * calculationg ODE by Euler method.
 * dx/dt = f(x).
 * x is real number.
 *
 * @Takaaki MINOMO
 */


#include <iostream>
#include <iomanip>

/* parameter */
template <typename T>
constexpr T alpha = static_cast<T>(1.);

template <typename T>
constexpr T beta = static_cast<T>(1.);

template <typename T>
constexpr T t_limit = static_cast<T>(10.);

template <typename T>
constexpr T dt = static_cast<T>(1./1000.);

/* initial value */
template <typename T>
constexpr T x_init = static_cast<T>(1./10.);

/* f(x) */
template<typename T>
constexpr T func(T x){
    return alpha<T> * (1. - x / beta<T>) * x;
}

int main(){
    std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    auto x = x_init<double>;
    auto t = 0.;

    /* initial value */
    std::cout << t << " " << x << std::endl;

    /* time loop */
    for(auto i = 0; t<t_limit<double>; i++){
        t = i*dt<double>;
        /* Euler method */
        x = x + dt<double> * func(x);
        std::cout << t << " " << x << std::endl;
    }
}
