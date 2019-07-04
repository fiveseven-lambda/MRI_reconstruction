#include <iostream>
#include <fstream>
#include <fftw3.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <cstring>
#include <cmath>
#include <climits>

int main(){
	int16_t size = 256;

	auto in = reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * size * size));

	std::ifstream ifs("P56320.cmp", std::ios::binary);
	for(size_t i = 0; i < size * size; ++i) ifs.read(reinterpret_cast<char *>(in[i]), sizeof in[i]);

	auto out = reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * size * size));
	auto plan = fftwf_plan_dft_2d(size, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	fftwf_free(in);
	
	int fd = open("out.bmp", O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	int filesize = size * size * 3 + 26;
	ftruncate(fd, filesize);
	auto map = reinterpret_cast<char *>(mmap(NULL, filesize, PROT_WRITE, MAP_SHARED, fd, 0));

	map[0] = 'B';
	map[1] = 'M';
	struct{
		int a:32, b:16, c:16, d:32, e:32, f:16, g:16, h:16, i:16;
	} header = {
		filesize,
		0, 0,
		26, 12,
		size, size, 1, 24
	};
	memcpy(map + 2, &header, sizeof header);

	float max = 0;
	for (size_t i = 0; i < size * size; ++i){
		float tmp = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
		if(max < tmp) max = tmp;
	}
	for(size_t i = 0; i < size * size; ++i){
		unsigned char tmp = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]) / max * UCHAR_MAX;
		for(size_t j = 0; j < 3; ++j) map[i * 3 + j + 26] = tmp;
	}

	fftwf_free(out);
	close(fd);
	munmap(map, filesize);
}
