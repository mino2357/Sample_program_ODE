#include <iostream>
#include <limits>
#include <iomanip>
#include <cmath>

namespace mino2357{
    template<typename T = double>
    class vector{
    private:
        T componentX;
        T componentY;
    public:
        vector(T x, T y) noexcept : componentX(x), componentY(y) {}
        vector() noexcept : vector{T{}, T{}} {}
        
        inline void setComponetX(T) noexcept;
        inline void setComponetY(T) noexcept;
        
        inline T getComponentX() const noexcept;
        inline T getComponentY() const noexcept;
        
        inline vector& operator=( const vector<T>&) noexcept;
        inline vector& operator+=(const vector<T>&) noexcept;
        inline vector& operator-=(const vector<T>&) noexcept;
        
        inline T norm() noexcept;
        
        inline void printPosition() noexcept;
    };
    
    template <typename T>
    inline void vector<T>::setComponetX(T x) noexcept {
        componentX = x;
    }
    
    template <typename T>
    inline void vector<T>::setComponetY(T y) noexcept {
        componentY = y;
    }
    
    template <typename T>
    inline T vector<T>::getComponentX() const noexcept {
        return this->componentX;
    }
    
    template <typename T>
    inline T vector<T>::getComponentY() const noexcept {
        return this->componentY;
    }
    
    template <typename T>
    inline void vector<T>::printPosition() noexcept{
        std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        std::cout << componentX << " " << componentY << std::endl;
    }
        
    template <typename T>
    inline vector<T>& vector<T>::operator=(const vector<T>& v) noexcept {
        this->componentX = v.getComponentX();
        this->componentY = v.getComponentY();
        return *this;
    }
    
    template <typename T>
    inline vector<T>& vector<T>::operator+=(const vector<T>& v) noexcept {
        this->componentX += v.getComponentX();
        this->componentY += v.getComponentY();
        return *this;
    }
        
    template <typename T>
    inline vector<T>& vector<T>::operator-=(const vector<T>& v) noexcept {
        this->componentX -= v.getComponentX();
        this->componentY -= v.getComponentY();
        return *this;
    }
    
    template <typename T>
    inline T vector<T>::norm() noexcept {
        return std::sqrt(componentX * componentX + componentY * componentY);
    }
    
    template <typename T>
    inline vector<T> operator+(const vector<T>& a, const vector<T>& b) noexcept {
        return vector<T>{
            a.getComponentX() + b.getComponentX(),
            a.getComponentY() + b.getComponentY(),
        };
    }

    template <typename T>
    inline vector<T> operator-(const vector<T>& a, const vector<T>& b) noexcept {
        return vector<T>{
            a.getComponentX() - b.getComponentX(),
            a.getComponentY() - b.getComponentY(),
        };
    }

    template <typename T>
    inline T operator*(const vector<T>& a, const vector<T>& b) noexcept {
        return
            a.getComponentX() * b.getComponentX()
          + a.getComponentY() * b.getComponentY();
    }

    template <typename T>
    inline vector<T> operator*(const T& a, const vector<T> v) noexcept {
        return vector<T>{a * v.getComponentX(), a * v.getComponentY()};
    }
}
