#include <iostream>
#include <complex>
#include <array>
#include <cmath>

using num = std::complex<double>;
constexpr int size = 256;

#define ARC(i, a, b) for(int i = (a); i < (b); ++i)
#define RC(i, n) ARC(i, 0, n)

template<int N, int S>
std::array<num, N> FFT(const std::array<num, N> &x){
	std::array<num, N / 2> half[2];
	for(int i = 0; i < N; ++i) half[i & 1][i >> 1] = x[i];

	half[0] = FFT<N / 2, S>(half[0]);
	half[1] = FFT<N / 2, S>(half[1]);

	num z = 1;
	constexpr num a(cos(2 * M_PI * S / N), sin(2 * M_PI * S / N));
	for(auto &i : half[1]){
		i *= z;
		z *= a;
	}
	
	std::array<num, N> ret;
	for(int i = 0; i < N / 2; ++i){
		ret[i] = half[0][i] + half[1][i];
		ret[i + N / 2] = half[0][i] - half[1][i];
	}

	return ret;
}

template<>
std::array<num, 1> FFT<1, 1>(const std::array<num, 1> &x){
	return x;
}

template<>
std::array<num, 1> FFT<1, -1>(const std::array<num, 1> &x){
	return x;
}

template<int N>
std::array<num, N> DFT(const std::array<num, N> &x){
	return FFT<N, +1>(x);
}

template<int N>
std::array<num, N> IDFT(const std::array<num, N> &x){
	std::array<num, N> ret = FFT<N, -1>(x);
	for(auto &i : ret) i /= N;
	return ret;
}


int main(){
	std::array<num, size> x, X;
	for(auto &i : x){
		double a, b;
		std::cin >> a >> b;
		i.real(a);
		i.imag(b);
	}
	X = DFT<size>(x);
	std::array<num, size> x_ = IDFT<size>(X);
}
