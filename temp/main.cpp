#include <iostream>

int main(){

    double a = 100;

    for(double i=1; i<10; i += 1.){
        std::cout << i << " " << a << std::endl;
        a = (i - 6) / i * a + 2 / i;
    }
}
