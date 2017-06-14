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

int main(){
    std::cout << std::fixed << std::setprecision(std::numeric_limits<f100>::digits10 + 1);
    
}
