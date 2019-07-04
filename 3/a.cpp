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
	int16_t size = 256;
	int slice = 24;

	auto in = reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * size * size));
	auto out = reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * size * size));
	auto plan = fftwf_plan_dft_2d(size, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	std::ifstream ifs("P56832.cmp", std::ios::binary);

	for(int i = 0; i < slice; ++i){
		for(size_t i = 0; i < size * size; ++i) ifs.read(reinterpret_cast<char *>(in[i]), sizeof in[i]);
		fftwf_execute(plan);

		std::stringstream filename;
		filename << "out" << std::setw(2) << std::setfill('0') << i << ".bmp";
		int fd = open(filename.str().c_str(), O_CREAT | O_RDWR, S_IWUSR | S_IRUSR | S_IRGRP | S_IROTH);
		int filesize = size * size * 3 + 26;
		ftruncate(fd, filesize);
		auto map = reinterpret_cast<char *>(mmap(NULL, filesize, PROT_WRITE, MAP_SHARED, fd, 0));
		float max_r = 0;
		for (size_t i = 0; i < size * size; ++i){
			float tmp;
		       	tmp = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
			if(max_r < tmp) max_r = tmp;
		}
		for(int i = 0; i < size; ++i){
			for(int j = 0; j < size; ++j){
				char tmp_r, tmp_x, tmp_y, tmp_theta;
				tmp_x = out[i * size + j][0] / max_r * CHAR_MAX + CHAR_MAX;
				tmp_y = out[i * size + j][1] / max_r * CHAR_MAX + CHAR_MAX;
				tmp_r = sqrt(out[i * size + j][0] * out[i * size + j][0] + out[i * size + j][1] * out[i * size + j][1]) / max_r * UCHAR_MAX;
				tmp_theta = atan2(out[i * size + j][1], out[i * size + j][0]) / M_PI * CHAR_MAX + CHAR_MAX;

				// swap the left and right half
				if(j < size / 2) j += size / 2;
				else j -= size / 2;

				map[(i * size + j) * 3] = tmp_theta;
				map[(i * size + j) * 3 + 1] = tmp_theta;
				map[(i * size + j) * 3 + 2] = tmp_theta;
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
	fftwf_destroy_plan(plan);
	fftwf_free(in);
	fftwf_free(out);
}
