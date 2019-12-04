#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <fftw3.h>
#include <cmath>
#include <ctype.h>

class Index{
	size_t val;
	const size_t n;
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
	Index shift(){
		if(val < n / 2) return Index(n, val + (n + 1) / 2);
		else return Index(n, val - n / 2);
	}
};

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

	auto in = reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * height * width));

	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j) for(Index k(2); k(); ++k) in[i, j][k] = kspace[i, j, k];

	close(kspace_fd);
	munmap(kspace, kspace_filesize);

	auto out = reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * height * width));
	fftwf_plan plan = fftwf_plan_dft_2d(height, width, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	fftwf_free(in);

	float abs[height][width], abs_max = 0;
	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j){
		abs[i][j] = hypot(out[i.shift(), j.shift()][0], out[i.shift(), j.shift()][1]);
		if(abs_max < abs[i][j]) abs_max = abs[i][j];
	}
	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j){
		abs[i][j] /= abs_max;
	}

	static char bmp_filename[256];
	sprintf(bmp_filename, "%s.bmp", kspace_filename);
	int bmp_filesize = height * width * 3 + 26;
	int bmp_fd = open(bmp_filename, O_CREAT | O_RDWR, S_IWUSR | S_IRUSR | S_IRGRP | S_IROTH);
	ftruncate(bmp_fd, bmp_filesize);
	auto bmp_map = reinterpret_cast<char *>(mmap(NULL, bmp_filesize, PROT_WRITE, MAP_SHARED, bmp_fd, 0));
	bmp_map[0] = 'B';
	bmp_map[1] = 'M';
	*reinterpret_cast<int32_t *>(bmp_map + 2) = bmp_filesize;
	*reinterpret_cast<int16_t *>(bmp_map + 6) = 0;
	*reinterpret_cast<int16_t *>(bmp_map + 8) = 0;
	*reinterpret_cast<int32_t *>(bmp_map + 10) = 26;
	*reinterpret_cast<int32_t *>(bmp_map + 14) = 12;
	*reinterpret_cast<int16_t *>(bmp_map + 18) = width;
	*reinterpret_cast<int16_t *>(bmp_map + 20) = height;
	*reinterpret_cast<int16_t *>(bmp_map + 22) = 1;
	*reinterpret_cast<int16_t *>(bmp_map + 24) = 24;
	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j) for(Index k(3); k(); ++k) bmp_map[(i, j, k) + 26] = abs[i][j] * (1u << 8);

	munmap(bmp_map, bmp_filesize);
	close(bmp_fd);
}
