#include <iostream>
#include <complex>

using Complex = std::complex<double>;

#define ARC(i, a, b) for(int i = (a); i < (b); ++i)
#define RC(i, n) ARC(i, 0, n)

template<int M, int N>
Complex zeta(){
	return Complex(cos(2 * M_PI * N / M), sin(2 * M_PI * N / M));
}

template<int N, int I = 0>
struct array_butterfly{
	static void butterfly(Complex x[]){
		Complex tmp0 = x[I], tmp1 = x[I + N];
		x[I] = tmp0 + tmp1;
		x[I + N] = (tmp0 - tmp1) * zeta<N * 2, I>();
		array_butterfly<N, I + 1>::butterfly(x + 1);
	}
};
template<int N>
struct array_butterfly<N, N>{
	static void butterfly(Complex x[]){
	}
};

template<int N>
struct array_FFT{
	static void FFT(Complex x[]){
		array_butterfly<N / 2>::butterfly(x);
		array_FFT<N / 2>::FFT(x);
		array_FFT<N / 2>::FFT(x + N / 2);
	}
};
template<>
struct array_FFT<1>{
	static void FFT(Complex x[]){}
};

int main(){
	constexpr int n = 4;
	Complex x[n], X[n];
	RC(i, n) x[i] = X[i] = 1 << (n - i + 1);
	RC(i, n) std::cout << x[i] << std::endl;
	std::cout << std::endl;
	array_FFT<n>::FFT(X);
	/*
	RC(i, n) std::cout << X[i] << std::endl;
	std::cout << std::endl;
	RC(i, n){
		Complex tmp;
		RC(j, n) tmp += x[j] * Complex(cos(2 * M_PI * i * j / n), sin(2 * M_PI * i * j / n));
		std::cout << tmp << std::endl;
	}*/
}
