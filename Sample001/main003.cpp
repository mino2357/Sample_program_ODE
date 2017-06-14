/*
 * 例外処理のテストプログラム
 * written by rittai3d
 * 
 * @mino2357
 */

#include <iostream>
#include <string>
#include <stdexcept>

constexpr int ratio(int f, int s){
    return s == 0 ? throw std::runtime_error {"ratio() : divided by zero."} : f / s;
}

int main(){
    for(std::string first{}, second{} ; ; ){
        std::cin >> first >> second;
        if(first == "q" || second == "Q"){
            break;
        }

        int result {};

        try{
            result = ratio(std::stoi(first), std::stoi(second));
        }catch(std::runtime_error& e){
            std::cerr << e.what() << std::endl;
            break;
        }

        std::cout << result << std::endl;
    }
}
