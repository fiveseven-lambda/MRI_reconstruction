#include <iostream>
#include <fstream>
#include <complex>
#include <array>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <unistd.h>
#include <fcntl.h>

template<int N>
class Block{
	char data[N];
	int cursor = 0;
public:
	constexpr uint32_t size(){return N;};
	Block &operator<<(const char *);
	template<class T>
	Block &operator<<(T);
	int write(int);
};


template<int N>
Block<N> &Block<N>::operator<<(const char *a){
	int i;
	for(i = 0; a[i]; ++i) data[cursor + i] = a[i];
	cursor += i;
	return *this;
}

template<int N>
template<class T>
Block<N> &Block<N>::operator<<(T a){
	memcpy(data + cursor, &a, sizeof a);
	cursor += sizeof a;
	return *this;
}

template<int N>
int Block<N>::write(int fd){
	::write(fd, data, N);
	return N;
}

using num = std::complex<double>;

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
	constexpr int size = 256;
	std::array<std::array<num, size>, size> x;
	static unsigned char bitmap_data[size * size * 3];
	std::fstream real("2DFSE_real.csv"), imag("2DFSE_imag.csv");
	for(auto &i : x) for(auto &j : i){
		double tmp;
		real >> tmp;
		j.real(tmp);
		imag >> tmp;
		j.imag(tmp);
	}
	for(int i = 0; i < size; ++i) x[i] = IDFT<size>(x[i]);
	for(int i = 0; i < size; ++i){
		std::array<num, size> tmp;
		for(int j = 0; j < size; ++j) tmp[j] = x[j][i];
		tmp = IDFT<size>(tmp);
		for(int j = 0; j < size; ++j) x[j][i] = tmp[j];
	}
	double max = DBL_MIN, min = DBL_MAX;
	for(auto i : x){
		for(auto j : i){
			if(max < abs(j)) max = abs(j);
			if(min > abs(j)) min = abs(j);
		}
	}
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			bitmap_data[(i * size + j) * 3] =
			bitmap_data[(i * size + j) * 3 + 1] =
			bitmap_data[(i * size + j) * 3 + 2] = abs(x[i][(j + size / 2) % size]) / max * 256;
		}
	}

	uint32_t width = 256, height = 256;
	Block<14> bitmapfileheader;
	Block<40> bitmapinfoheader;
	int fd = creat("out.bmp", S_IWUSR | S_IRUSR);
	bitmapfileheader << "BM";
	bitmapfileheader << bitmapfileheader.size() + bitmapinfoheader.size() + width * height * 3;
	bitmapfileheader << uint16_t(0) << uint16_t(0);
	bitmapfileheader << bitmapfileheader.size() + bitmapinfoheader.size();
	bitmapinfoheader << bitmapinfoheader.size() << width << height << uint16_t(1) << uint16_t(24) << uint32_t(0) << width * height * 3 << uint32_t(0) << uint32_t(0) << uint32_t(0) << uint32_t(0);

	bitmapfileheader.write(fd);
	bitmapinfoheader.write(fd);
	write(fd, bitmap_data, width * height * 3);
}
