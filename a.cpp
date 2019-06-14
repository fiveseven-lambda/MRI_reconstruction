#include <iostream>
#include <complex>
#include <array>
#include <random>

using Real = double;
using Complex = std::complex<Real>;

template<int M, int N>
constexpr Complex zeta(){
	return Complex(cos(2 * M_PI * N / M), sin(2 * M_PI * N / M));
}

template<int N, int SIGN, int I = 0>
struct array_mulzeta{
	static void mulzeta(Complex x[]){
		*x *= zeta<N * 2, SIGN * I>();
		array_mulzeta<N, SIGN, I + 1>::mulzeta(x + 1);
	}
};
template<int N, int SIGN>
struct array_mulzeta<N, SIGN, N>{
	static void mulzeta(Complex []){}
};
template<int SIZE, int I = 0>
struct array{
	static void add(Complex a[], Complex b[], Complex c[]){
		*c = *a + *b;
		array<SIZE, I + 1>::add(a + 1, b + 1, c + 1);
	}
	static void subst(Complex a[], Complex b[], Complex c[]){
		*c = *a - *b;
		array<SIZE, I + 1>::subst(a + 1, b + 1, c + 1);
	}
	static void div(Complex x[]){
		*x /= SIZE;
		array<SIZE, I + 1>::div(x + 1);
	}
};
template<int SIZE>
struct array<SIZE, SIZE>{
	static void add(Complex [], Complex [], Complex []){}
	static void subst(Complex [], Complex [], Complex []){}
	static void div(Complex []){};
};

template<int SIZE, int SIGN, int GAP = 1>
struct array_FFT{
	static void FFT(Complex x[], Complex X[]){
		Complex even[SIZE / 2], odd[SIZE / 2];
		array_FFT<SIZE / 2, SIGN, GAP * 2>::FFT(x, even);
		array_FFT<SIZE / 2, SIGN, GAP * 2>::FFT(x + GAP, odd);
		array_mulzeta<SIZE / 2, SIGN>::mulzeta(odd);
		array<SIZE / 2>::add(even, odd, X);
		array<SIZE / 2>::subst(even, odd, X + SIZE / 2);
	}
};
template<int SIGN, int GAP>
struct array_FFT<1, SIGN, GAP>{
	static void FFT(Complex x[], Complex X[]){
		*X = *x;
	}
};

template<int N>
void DFT(Complex x[], Complex X[]){
	array_FFT<N, 1>::FFT(x, X);
}
template<int N>
void IDFT(Complex X[], Complex x[]){
	array_FFT<N, -1>::FFT(X, x);
	array<N>::div(x);
}

int main(){
	constexpr int size = 16;
	Complex x[size], X[size];

	// substitute random complex numbers to x
	std::random_device gen;
	std::uniform_real_distribution<Real> rnd(-1., 1.);
	for(int i = 0; i < size; ++i){
		x[i].real(rnd(gen));
		x[i].imag(rnd(gen));
	}

	// copy x to X
	for(int i = 0; i < size; ++i) X[i] = x[i];


	// DFT ( FFT, time complexity : O(size * log(size)) )
	DFT<size>(X, X);

	for(int i = 0; i < size; ++i){

		// calculate DFT by definition ( time complexity: O(size * size) )
		Complex tmp = 0;
		for(int j = 0; j < size; ++j) tmp += x[j] * Complex(cos(2 * M_PI * i * j / size), sin(2 * M_PI * i * j / size));

		// check the answer
		std::cout << X[i] << '=' << tmp << std::endl;
	}

	std::cout << std::endl;

	// IDFT ( FFT )
	IDFT<size>(X, X);

	for(int i = 0; i < size; ++i){

		// check the answer (IDFT is the inverse transform of DFT, so X = IDFT(DFT(x)) should be equal to x)
		std::cout << X[i] << '=' << x[i] << std::endl;
	}
}
