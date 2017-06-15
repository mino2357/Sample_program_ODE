/*
 * written by rittai3d
 *
 */

#include <iostream>
#include <string>
#include <stdexcept>
#include <type_traits>

/*
template <bool B, class T = void>
using std::enable_if_t = typename std::enable_if<B, T>::type;

template <typename T>
inline constexpr bool std::is_integral_v = std::is_integral<T>::value;
*/

template <
    typename T,
    typename std::enable_if_t<
        std::is_integral_v<T>, T
    >*  = nullptr
>
constexpr  T ratio(T const& a, T const& b){
    return b == static_cast<T>(0.) ? throw std::runtime_error{"devided by zero!"} : a / b;
}

int main(){
    // OK
    static_assert(ratio(  1, 1), "");
    // NG
    //static_assert(ratio(1.0, 1), "");
}
