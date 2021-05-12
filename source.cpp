// Own_FFT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//source: http://www.cs.cmu.edu/afs/andrew/scs/cs/15-463/2001/pub/www/notes/fourier/fourier.pdf

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







template <typename T, int size>
void DFT(std::unique_ptr<std::array<std::complex< T>, size>> &X, std::unique_ptr<std::array<std::complex< T>, size>> &Y) {

    const std::complex<double> Con(0.0, (-2.0 / size) * M_PI);
    

    for (double xx = 0; xx <size; xx++)
    {
        
        std::complex<double> sum(0, 0);
        for (double yy = 0; yy < size; yy++)
        {
            sum += X->at(yy)*std::exp(xx * yy * Con);
        }
        Y->at(xx) = sum;

    }
}



template<typename T>
void Binary_Read(const char* filename, T* buffer, const int size) {
    std::ifstream myFile(filename, std::ios::in | std::ios::binary);
    if (!myFile) {
        std::cout << "Cannot open file!\n";

    }
    else {

        myFile.read(reinterpret_cast<char*>(buffer), size * sizeof(T));

    }
}

template<typename T>
void Binary_Write(const char* filename, T* buffer, const int size) {
    std::ofstream myFile(filename, std::ios::out | std::ios::binary);
    if (!myFile) {
        std::cout << "Cannot open file!\n";

    }
    else {

        myFile.write(reinterpret_cast<char*>(buffer), size * sizeof(T));

    }
}

template<size_t X, size_t Y=(size_t)std::log2(X)>
void bit_reverse(std::array<int, X> &Out) {  
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

template<size_t X, size_t Y = (size_t)std::log2(X)>
void bit_reverse(std::unique_ptr<std::array<int, X>>& Out) {

    size_t amount_step = 0;
    if (Y % 2 == 0) {
        amount_step = Y / 2;
    }
    else {
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
            I_value[i1] = (bool)temp;        }
        //Out->at(i) = static_cast<int>(I_value.to_ulong());
        Out->at(i) = (int)(I_value.to_ulong());

        
    }

}



//template<typename T, int size>
//using SPCA = std::unique_ptr < std::array< std::complex<T>, size>>; //Smart Pointer Complex Array
template<typename T, int size, int LoS = static_cast<int>(std::log2(size))>
void FFT(const std::unique_ptr<std::array<std::complex< T>, size>>& in, std::unique_ptr<std::array<std::complex< T>, size>>& out) {
    //constexpr int LoS =(int)std::log2(size); //logarithm of Size
    std::unique_ptr<std::array<int, size>> bit_place{ new std::array<int, size> };
    std::cout << "size of : " << sizeof(bit_place) << "\n";
    auto start = std::chrono::steady_clock::now();
    bit_reverse<size,LoS>(bit_place);
    std::cout << "bit_reverse processing time " << std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - start).count() << " ms\n";
    start = std::chrono::steady_clock::now();
    std::unique_ptr < std::array< std::complex<T>, size>> temp{ new std::array< std::complex<T>, size> };
    std::unique_ptr < std::array< std::complex<T>, size>> W_arr{ new std::array< std::complex<T>, size> };
    W_arr->at(0) = std::complex<T>(1, 0);
    W_arr->at(1) = std::exp(std::complex<T>(0.0, (-2.0 / size) * M_PI));
    for (size_t i = 2; i < size; i++)
    {
        W_arr->at(i) = W_arr->at(i - 1) * W_arr->at(1);
    }
    std::cout << "W_arr processing time " << std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - start).count() << " ms\n";

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
                temp->at(x) = in->at(bit_place->at(k)) + W_arr->at(j) * in->at(bit_place->at(h));

            }
        }
        else if (LoS % 2 == 1 && y == 0) {
            for (int x = 0; x < size; x++) {
                k = std::floor((x / Pow2_y1)) * Pow2_y1 + x % Pow2_y;
                h = k + Pow2_y;
                j = std::floor(x % Pow2_y1) * (size / Pow2_y1);
                out->at(x) = in->at(bit_place->at(k)) + W_arr->at(j) * in->at(bit_place->at(h));

            }
        }
        else if((LoS % 2 == 0 && y % 2 == 1) || (LoS % 2 == 1 && y % 2 == 0) ) {
            for (int x = 0; x < size; x++) {
                k = std::floor((x / Pow2_y1)) * Pow2_y1 + x % Pow2_y;
                h = k + Pow2_y;
                j = std::floor(x % Pow2_y1) * (size / Pow2_y1);
                out->at(x) = temp->at(k) + W_arr->at(j) * temp->at(h);
            }}
        else {
            for (int x = 0; x < size; x++) {
                k = std::floor((x / Pow2_y1)) * Pow2_y1 + x % Pow2_y;
                h = k + Pow2_y;
                j = std::floor(x % Pow2_y1) * (size / Pow2_y1);
                temp->at(x) = out->at(k) + W_arr->at(j) * out->at(h);

            }
        }            
    }
}








int main()
{
    const int size_ = 1024*1; // FFT size
    const int LoS = 10; // User should change here according to size, LoS=log2(1024)
    double* buffer = new double[size_];
    Binary_Read("Signal.bin", buffer, size_);
    //std::array<std::complex< double>, size_> X{};
    std::unique_ptr<std::array<std::complex< double>, size_>> X{ new std::array<std::complex< double>, size_> }; //input
    std::unique_ptr<std::array<std::complex< double>, size_>> Y{ new std::array<std::complex< double>, size_> }; // output

    X.get()->fill(0);
    for (size_t i = 0; i < size_; i++) // moving data from buffer to input array
    {
        X.get()->at(i).real(buffer[i]);
        
    }
    
    
    for (size_t i = 0; i < 10; i++) //first ten samples
    {
        std::cout << X->at(i).real() << " ," << X->at(i).imag() << "\n";
    }

    auto start = std::chrono::steady_clock::now();
    DFT<double, size_>(X, Y);
    std::cout << "DFT processing time:  " << std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - start).count() << " ms\n";
    for (size_t i = 0; i < 10; i++) output first 10 sample
    {
        std::cout << Y->at(i).real() << " ," << Y->at(i).imag() << "\n";
    }


    start = std::chrono::steady_clock::now();
    FFT<double, size_, LoS>(X, Y);
    std::cout<<"FFT processing time: " << std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now()- start).count() << " ms\n";

    for (size_t i = 0; i < 10; i++) output first 10 sample
    {
        std::cout << Y->at(i).real() << " ," << Y->at(i).imag() << "\n";
    }

    delete[] buffer;
}

