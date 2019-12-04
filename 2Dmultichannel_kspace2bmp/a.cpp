#include <stdio.h>
#include <gdcm-3.0/gdcmImageWriter.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <fftw3.h>
#include <cmath>
#include <ctype.h>
#include <string.h>

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
	char kspace_filename[256];
	strcpy(kspace_filename, argv[1]);
	int kspace_filename_offset = strlen(argv[1]);
	unsigned height, width;
	sscanf(argv[2], "%u", &height);
	if(height == 0) fprintf(stderr, "error: invalid height\n");
	sscanf(argv[3], "%u", &width);
	if(width == 0) fprintf(stderr, "error: invalid width\n");
	unsigned kspace_filesize = sizeof(float) * height * width * 2;

	float ans[height][width] = {}, ans_max = 0;

	for(int channel = 0; channel < 12; ++channel){
		sprintf(kspace_filename + kspace_filename_offset, "%d", channel);

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

		for(Index i(height); i(); ++i) for(Index j(width); j(); ++j){
			ans[i][j] += out[i.shift(), j.shift()][0] * out[i.shift(), j.shift()][0] + out[i.shift(), j.shift()][1] * out[i.shift(), j.shift()][1];
		}
	}
	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j) ans[i][j] = sqrt(ans[i][j]);
	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j) if(ans_max < ans[i][j]) ans_max = ans[i][j];
	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j) ans[i][j] /= ans_max;

	static uint16_t dcm_buf[512][308];
	for(Index i(height); i(); ++i) for(Index j(width); j(); ++j) dcm_buf[i][j] = ans[i][j] * UINT16_MAX;
	gdcm::ImageWriter writer;
	gdcm::Image &image = writer.GetImage();
	image.SetNumberOfDimensions(2);
	unsigned int dims[2] = {width, height};
	image.SetDimensions(dims);
	gdcm::PixelFormat pf = gdcm::PixelFormat::UINT16;
	pf.SetSamplesPerPixel(4);
	image.SetPixelFormat(pf);
	image.SetPhotometricInterpretation(gdcm::PhotometricInterpretation::ARGB);
	image.SetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);
	static char filename[] = "out.dcm";
	writer.SetFileName(filename);
	gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));
	pixeldata.SetByteValue((const char *)(dcm_buf), (uint32_t)2 * height * width);
	image.SetDataElement(pixeldata);
	writer.Write();


	static char bmp_filename[256] = "out.bmp";
	int bmp_filesize = height * width * 3 * 3 + 26;
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
	*reinterpret_cast<int16_t *>(bmp_map + 18) = width * 3;
	*reinterpret_cast<int16_t *>(bmp_map + 20) = height;
	*reinterpret_cast<int16_t *>(bmp_map + 22) = 1;
	*reinterpret_cast<int16_t *>(bmp_map + 24) = 24;
	for(Index i(height); i(); ++i) for(Index j(width * 3); j(); ++j) for(Index k(3); k(); ++k) bmp_map[(i, j, k) + 26] = ans[i][width - j % width] * (1u << 8);

	munmap(bmp_map, bmp_filesize);
	close(bmp_fd);
}
