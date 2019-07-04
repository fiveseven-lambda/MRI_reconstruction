#include <iostream>
#include <complex>
#include <fstream>
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

using Complex = std::complex<double>;

template<int M, int N>
constexpr Complex zeta(){
	return Complex(cos(2 * M_PI * N / M), sin(2 * M_PI * N / M));
}

template<int N, int GAP, int I>
struct array_butterfly{
	static void butterfly(Complex x[]){
		Complex tmp0 = x[0], tmp1 = x[GAP * N / 2];
		x[0] = tmp0 + tmp1;
		x[GAP * N / 2] = (tmp0 - tmp1) * zeta<N, I / 2>();
		array_butterfly<N, GAP, I + 2>::butterfly(x + GAP);
	}
};
template<int N, int GAP>
struct array_butterfly<N, GAP, N>{
	static void butterfly(Complex x[]){
	}
};

template<int N, int GAP>
struct array_FFT{
	static void FFT(Complex x[]){
		array_butterfly<N, GAP, 0>::butterfly(x);
		array_FFT<N / 2, GAP>::FFT(x);
		array_FFT<N / 2, GAP>::FFT(x + GAP * N / 2);
	}
};
template<int GAP>
struct array_FFT<1, GAP>{
	static void FFT(Complex x[]){
	}
};

int rev8bit(int a){
	a = ((a & 0x55) << 1) | ((a & 0xAA) >> 1);
	a = ((a & 0x33) << 2) | ((a & 0xCC) >> 2);
	return ((a & 0x0F) << 4) | (a >> 4);
}

int main(){
	/*
	constexpr int n = 8;
	Complex x[n * 2], X[n * 2];
	RC(i, n) x[i * 2] = X[i * 2] = 1 << (n - i + 1);
	array_FFT<n, 2>::FFT(X);
	RC(i, n) std::cout << X[i * 2] << std::endl;
	std::cout << std::endl;
	RC(i, n){
		Complex tmp;
		RC(j, n) tmp += x[j * 2] * Complex(cos(2 * M_PI * i * j / n), sin(2 * M_PI * i * j / n));
		std::cout << tmp << std::endl;
	}
	*/
	constexpr int size = 256;
	Complex x[size][size];
	Complex X[size];
	std::fstream real("2DFSE_real.csv"), imag("2DFSE_imag.csv");
	static unsigned char bitmap_data[size * size * 3];
	for(int i = 0; i < size; ++i) for(int j = 0; j < size; ++j){
		double tmp;
		real >> tmp;
		x[i][j].real(tmp);
		imag >> tmp;
		x[i][j].imag(tmp);
	}
	for(int i = 0; i < size; ++i) array_FFT<size, 1>::FFT(&x[i][0]);
	for(int i = 0; i < size; ++i) array_FFT<size, size>::FFT(&x[0][i]);
	double max = DBL_MIN;
	for(int i = 0; i < size; ++i) for(int j = 0; j < size; ++j) if(max < abs(j)) max = abs(j);
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			bitmap_data[(i * size + j) * 3] =
			bitmap_data[(i * size + j) * 3 + 1] =
			bitmap_data[(i * size + j) * 3 + 2] = abs(x[rev8bit(i)][rev8bit(j)]) / max * 256;
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
