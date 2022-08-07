// Own_FFT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//source: http://www.cs.cmu.edu/afs/andrew/scs/cs/15-463/2001/pub/www/notes/fourier/fourier.pdf

#include "fouriertransform.hpp"





using namespace FourierTransform;








int main()
{

    const int size_ = 1024
    ; // FFT size
    double* buffer = new double[size_];

    for (int i = 0; i < size_; i++)
    {
        buffer[i]=i;
    }

    //std::array<std::complex< double>, size_> X{};
    std::array<std::complex< double>, size_> X{}; //input
    std::array<std::complex< double>, size_> Y{}; // output

    X.fill(0);
    for (size_t i = 0; i < size_; i++) // moving data from buffer to input array
    {
        X[i].real(buffer[i]);
        
    }
    
    
    for (size_t i = 0; i < 10; i++) //first ten samples
    {
        std::cout << X[i].real() << " ," << X[i].imag() << "\n";
    }

    auto start = std::chrono::steady_clock::now();
    DFT<double, size_>(X, Y);
    std::cout << "DFT processing time:  " << std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - start).count() << " ms\n";
    /*for (size_t i = 0; i < 10; i++) //output first 10 sample
    {
        std::cout << Y[i].real() << " ," << Y[i].imag() << "\n";
    }*/


    start = std::chrono::steady_clock::now();
    FFT<double, size_>(X, Y);
    std::cout<<"FFT processing time: " << std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now()- start).count() << " ms\n";

   /* for (size_t i = 0; i < 10; i++) //output first 10 sample
    {
        std::cout << Y[i].real() << " ," << Y[i].imag() << "\n";
    }*/

    delete[] buffer;
}

