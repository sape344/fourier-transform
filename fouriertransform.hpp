#if !defined(FOURIERTRANSFPRM_HPP)
#define FOURIERTRANSFPRM_HPP

#define _USE_MATH_DEFINES
#include<fstream>
#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include <complex>
#include <array>
#include<memory>
#include<chrono>
#include<bitset>

namespace FourierTransform
{

    

template <typename T, int size>
inline void DFT(std::array<std::complex< T>, size> &X, std::array<std::complex< T>, size> &Y) {

    constexpr std::complex<double> Con(0.0, (-2.0 / size) * M_PI);
    for (double xx = 0; xx <size; xx++)
    {
        Y[xx].imag(0);
        Y[xx].real(0);

        for (double yy = 0; yy < size; yy++)
        {
            Y[xx] += X[yy]*std::exp(xx * yy * Con);
        }

    }
}



template<typename T>
inline void Binary_Read(const char* filename, T* buffer, const int size) {
    std::ifstream myFile(filename, std::ios::in | std::ios::binary);
    if (!myFile) {
        std::cout << "Cannot open file!\n";

    }
    else {

        myFile.read(reinterpret_cast<char*>(buffer), size * sizeof(T));

    }
}

template<typename T>
inline void Binary_Write(const char* filename, T* buffer, const int size) {
    std::ofstream myFile(filename, std::ios::out | std::ios::binary);
    if (!myFile) {
        std::cout << "Cannot open file!\n";

    }
    else {

        myFile.write(reinterpret_cast<char*>(buffer), size * sizeof(T));

    }
}

template<size_t X, size_t Y=(size_t)std::log2(X)>
inline void bit_reverse(std::array<int, X> &Out) {  
    //size_t Y = (size_t)std::log2(X);
    size_t amount_step = 0;
    if (Y % 2 == 0) {
        amount_step = Y / 2;
    }
    else{
        amount_step = (Y - 1) / 2;
    }
    std::bitset<Y> I_value = 0;
    for (size_t i = 0; i < X; i++)
        {
        I_value = i;
            for (size_t i1 = 0; i1 < amount_step; i1++)
            {
                bool temp = (bool)I_value[Y - i1 - 1];
                I_value[Y - i1 - 1] = (bool)I_value[i1];
                I_value[i1] = (bool)temp;                
            }
            Out[i] = (int)(I_value.to_ulong());
        }   

}



//template<typename T, int size>
//using SPCA = std::unique_ptr < std::array< std::complex<T>, size>>; //Smart Pointer Complex Array
template<typename T, int size, int LoS = static_cast<int>(std::log2(size))>
inline void FFT(const std::array<std::complex< T>, size>& in, std::array<std::complex< T>, size>& out) {
    //constexpr int LoS =(int)std::log2(size); //logarithm of Size
    std::array<int, size> bit_place{};
  //  std::cout << "size of : " << sizeof(bit_place) << "\n";
 //   auto start = std::chrono::steady_clock::now();
    bit_reverse<size>(bit_place);
  //  std::cout << "bit_reverse processing time " << std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - start).count() << " ms\n";
   // start = std::chrono::steady_clock::now();
    std::array< std::complex<T>, size> temp{};
    std::array< std::complex<T>, size> W_arr{};
    W_arr[0] = std::complex<T>(1, 0);
    W_arr[1] = std::exp(std::complex<T>(0.0, (-2.0 / size) * M_PI));
    for (size_t i = 2; i < size; i++)
    {
        W_arr[i] = W_arr[i - 1] * W_arr[1];
    }
    //std::cout << "W_arr processing time " << std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - start).count() << " ms\n";

    int k = 0;
    int h = 0;
    int j = 0;
    for (int y = 0; y < LoS; y++)
    {
        int Pow2_y1 = std::pow(2, y + 1);
        int Pow2_y = std::pow(2, y);
        if (LoS % 2 == 0 && y == 0) {
            for (int x = 0; x < size; x++) {
                k = std::floor((x / Pow2_y1)) * Pow2_y1 + x % Pow2_y;
                h = k + Pow2_y;
                j = std::floor(x % Pow2_y1) * (size / Pow2_y1);
                temp[x] = in[bit_place[k]] + W_arr[j] * in[bit_place[h]];

            }
        }
        else if (LoS % 2 == 1 && y == 0) {
            for (int x = 0; x < size; x++) {
                k = std::floor((x / Pow2_y1)) * Pow2_y1 + x % Pow2_y;
                h = k + Pow2_y;
                j = std::floor(x % Pow2_y1) * (size / Pow2_y1);
                out[x] = in[bit_place[k]] + W_arr[j] * in[bit_place[h]];

            }
        }
        else if((LoS % 2 == 0 && y % 2 == 1) || (LoS % 2 == 1 && y % 2 == 0) ) {
            for (int x = 0; x < size; x++) {
                k = std::floor((x / Pow2_y1)) * Pow2_y1 + x % Pow2_y;
                h = k + Pow2_y;
                j = std::floor(x % Pow2_y1) * (size / Pow2_y1);
                out[x] = temp[k] + W_arr[j] * temp[h];
            }}
        else {
            for (int x = 0; x < size; x++) {
                k = std::floor((x / Pow2_y1)) * Pow2_y1 + x % Pow2_y;
                h = k + Pow2_y;
                j = std::floor(x % Pow2_y1) * (size / Pow2_y1);
                temp[x] = out[k] + W_arr[j] * out[h];

            }
        }            
    }
}

}



#endif // FOURIERTRANSFPRM_HPP
