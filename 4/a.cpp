#include <fftw3.h>
#include <cstdint>
#include <cmath>
#include <stdio.h>
#include <gdcm-3.0/gdcmImageWriter.h>

int main(){
	constexpr uint16_t size = 256;
	constexpr unsigned int slice = 24;

	fftwf_complex
		*in = reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * size * size * slice)),
		*out = reinterpret_cast<fftwf_complex *>(fftwf_malloc(sizeof(fftwf_complex) * size * size * slice));
	fftwf_plan plan = fftwf_plan_dft_3d(slice, size, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	std::ifstream ifs("k-space", std::ios::binary);
	for(size_t i = 0; i < size * size * slice; ++i) ifs.read(reinterpret_cast<char *>(in[i]), sizeof in[i]);

	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	fftwf_free(in);

	static float fbuf[slice][size][size], fbuf_max = -1;
	for(size_t i = 0; i < size * size * slice; ++i){
		reinterpret_cast<float *>(fbuf)[i] = sqrtf(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
		if(reinterpret_cast<float *>(fbuf)[i] > fbuf_max) fbuf_max = reinterpret_cast<float *>(fbuf)[i];
	}

	static uint32_t buf[slice][size][size];
	for(size_t i = 0; i < slice; ++i) for(size_t j = 0; j < size; ++j) for(size_t k = 0; k < size; ++k){
		buf[i][j < size / 2 ? j + size / 2 : j - size / 2][k < size / 2 ? k + size / 2 : k - size / 2] = fbuf[i][j][k] / fbuf_max * UINT32_MAX;
	}
	fftwf_free(out);

	gdcm::ImageWriter writer;

	gdcm::Image &image = writer.GetImage();

	image.SetNumberOfDimensions(2);

	unsigned int dims[2] = {size, size};
	image.SetDimensions(dims);

	gdcm::PixelFormat pf = gdcm::PixelFormat::UINT32;
	pf.SetSamplesPerPixel(4);
	image.SetPixelFormat(pf);

	image.SetPhotometricInterpretation(gdcm::PhotometricInterpretation::ARGB);

	image.SetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);

	for(size_t i = 0; i < slice; ++i){
		static char filename[10];
		sprintf(filename, "out%02lu.dcm", i);
		writer.SetFileName(filename);

		gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));

		pixeldata.SetByteValue((const char *)(buf + i), (uint32_t)4 * size * size);

		image.SetDataElement(pixeldata);

		if(!writer.Write()) return -1;
	}

	return 0;
}
