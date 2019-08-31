#include <iostream>
#include <fstream>
#include <fftw3.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <cstring>
#include <cmath>
#include <climits>
#include <sstream>
#include <ios>
#include <iomanip>

int main(){
	uint16_t size = 256;
	unsigned int slice = 24;

	auto in = reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * size * size * size));
	auto out = reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * size * size * size));
	auto plan = fftwf_plan_dft_3d(slice, size, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	std::ifstream ifs("k-space", std::ios::binary);

	for(size_t i = 0; i < size * size * slice; ++i) ifs.read(reinterpret_cast<char *>(in[i]), sizeof in[i]);

	fftwf_execute(plan);

	fftwf_destroy_plan(plan);

	float max_r = 0;
	for (size_t i = 0; i < size * size * slice; ++i){
		float tmp;
		tmp = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
		if(max_r < tmp) max_r = tmp;
	}

	for(size_t i = 0; i < slice; ++i){
		std::stringstream filename;
		filename << "out" << std::setw(2) << std::setfill('0') << i << ".bmp";
		int fd = open(filename.str().c_str(), O_CREAT | O_RDWR, S_IWUSR | S_IRUSR | S_IRGRP | S_IROTH);
		int filesize = size * size * 3 + 26;
		ftruncate(fd, filesize);
		auto map = reinterpret_cast<char *>(mmap(NULL, filesize, PROT_WRITE, MAP_SHARED, fd, 0));
		for(int j = 0; j < size; ++j){
			for(int k = 0; k < size; ++k){
				char tmp_r = sqrt(
					out[i * size * size + j * size + k][0]
					* out[i * size * size + j * size + k][0]
					+ out[i * size * size + j * size + k][1]
					* out[i * size * size + j * size + k][1]
				) / max_r * UCHAR_MAX;

				// swap the left and right half
				if(j < size / 2) j += size / 2;
				else j -= size / 2;
				if(k < size / 2) k += size / 2;
				else k -= size / 2;

				map[(j * size + k) * 3] = tmp_r;
				map[(j * size + k) * 3 + 1] = tmp_r;
				map[(j * size + k) * 3 + 2] = tmp_r;
			}
		}
		map[0] = 'B';
		map[1] = 'M';
		*reinterpret_cast<int32_t *>(map + 2) = filesize;
		*reinterpret_cast<int16_t *>(map + 6) = 0;
		*reinterpret_cast<int16_t *>(map + 8) = 0;
		*reinterpret_cast<int32_t *>(map + 10) = 26;
		*reinterpret_cast<int32_t *>(map + 14) = 12;
		*reinterpret_cast<int16_t *>(map + 18) = size;
		*reinterpret_cast<int16_t *>(map + 20) = size;
		*reinterpret_cast<int16_t *>(map + 22) = 1;
		*reinterpret_cast<int16_t *>(map + 24) = 24;
	}

	fftwf_free(in);
	fftwf_free(out);
}
