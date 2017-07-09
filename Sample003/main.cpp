#include <iostream>
#include "myvector.hpp"

constexpr double k      = 1.;
constexpr double m      = 1.;
constexpr double x_init = 1.;
constexpr double v_init = 0.;

constexpr double dt      = 0.001;
constexpr double t_limit = 20.0;

mino2357::vector<> func(mino2357::vector<>& x){
    return mino2357::vector<>{x.getComponentY(), - k / m * x.getComponentX()};
}

int main(){
    mino2357::vector<> x(x_init, v_init);

    double t {};

    for(size_t i {}; t<t_limit; ++i){
        std::cout << t << " " << x.getComponentX() << " " << x.getComponentY() << std::endl;
        
        x = x + dt * func(x);
        t = i * dt;
    }
}
