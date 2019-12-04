#include <iostream>
#include <fftw3.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <math.h>
#include <unistd.h>

class Index{
	const size_t n;
	size_t val;
public:
	constexpr Index(const size_t n, size_t val = 0): n(n), val(val){}
	bool operator()(size_t r = 0){
		return val < n - r;
	}
	operator size_t() const {
		return val;
	}
	operator size_t &(){
		return val;
	}
	Index operator,(const Index &i) const {
		return Index(n, val * i.n + i.val);
	}
private:
	template<int I>
	Index shift_() const {
		long tmp = val - (n + I) / 2;
		if(tmp < 0) tmp += n;
		return Index(n, tmp);
	}
public:
	Index shift() const { return shift_<0>(); }
	Index ishift() const { return shift_<1>(); }
};

struct fftwf_space{
	fftwf_complex *in, *out;
	fftwf_plan plan;
	int height, width;
	fftwf_space() = delete;
	fftwf_space(int height, int width, int sign, unsigned flags = FFTW_ESTIMATE):
		in(reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * height * width))),
		out(reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * height * width))),
		plan(fftwf_plan_dft_2d(height, width, in, out, sign, flags)) {}
	~fftwf_space(){
		fftwf_destroy_plan(plan);
		fftwf_free(in);
		fftwf_free(out);
	}
	void execute(){
		fftwf_execute(plan);
	}
};

void bitmap(const char filename[], float data[], int height, int width){
	int fd = open(filename, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	int filesize = height * width * 3 + 26;
	ftruncate(fd, filesize);
	auto map = reinterpret_cast<unsigned char *>(mmap(NULL, filesize, PROT_WRITE, MAP_SHARED, fd, 0));
	map[0] = 'B';
	map[1] = 'M';
	*reinterpret_cast<int32_t *>(map + 2) = filesize;
	*reinterpret_cast<int16_t *>(map + 6) = 0;
	*reinterpret_cast<int16_t *>(map + 8) = 0;
	*reinterpret_cast<int32_t *>(map + 10) = 26;
	*reinterpret_cast<int32_t *>(map + 14) = 12;
	*reinterpret_cast<int16_t *>(map + 18) = width;
	*reinterpret_cast<int16_t *>(map + 20) = height;
	*reinterpret_cast<int16_t *>(map + 22) = 1;
	*reinterpret_cast<int16_t *>(map + 24) = 24;
	for(Index i(height); i(); ++i){
		for(Index j(width); j(); ++j){
			unsigned char tmp = data[i, Index(width, width - j)] * ((1 << 8) - 1);
			for(Index k(3); k(); ++k){
				map[(i, j, k) + 26] = tmp;
			}
		}
	}
	close(fd);
	munmap(map, filesize);
}

int main(int argc, char *argv[]){
	if(argc < 4){
		printf("usage: %s <path to k-space> <height> <width>\n", argv[0]);
		return -1;
	}
	char *kspace_filename = argv[1];
	unsigned height, width;
	sscanf(argv[2], "%u", &height);
	if(height == 0) fprintf(stderr, "error: invalid height\n");
	sscanf(argv[3], "%u", &width);
	if(width == 0) fprintf(stderr, "error: invalid width\n");
	unsigned kspace_filesize = sizeof(float) * height * width * 2;
	int kspace_fd = open(kspace_filename, O_RDONLY);
	if(kspace_fd == -1) fprintf(stderr, "error: cannot open file \'%s\'\n", kspace_filename);
	float *kspace = reinterpret_cast<float *>(mmap(NULL, kspace_filesize, PROT_READ, MAP_PRIVATE, kspace_fd, 0));

	fftwf_space space(height, width, FFTW_FORWARD);
	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j) for(Index k(2); k(); ++k) space.in[i.shift(), j.shift()][k] = kspace[i, j, k];

	close(kspace_fd);
	munmap(kspace, kspace_filesize);

	space.execute();

	float abs[height][width], abs_max = 0;
	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j){
		abs[i][j] = hypot(space.out[i.shift(), j.shift()][0], space.out[i.shift(), j.shift()][1]);
		if(abs_max < abs[i][j]) abs_max = abs[i][j];
	}
	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j){
		abs[i][j] /= abs_max;
	}

	bitmap("out.bmp", reinterpret_cast<float *>(abs), height, width);
}
