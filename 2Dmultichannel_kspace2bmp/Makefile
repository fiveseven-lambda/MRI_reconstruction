build: a.cpp
	g++ a.cpp -lfftw3f -lgdcmcharls -lgdcmCommon -lgdcmDICT -lgdcmDSED -lgdcmexpat -lgdcmIOD -lgdcmjpeg12 -lgdcmjpeg16 -lgdcmjpeg8 -lgdcmMEXD -lgdcmMSFF -lgdcmuuid -lgdcmzlib
test: a.out
	./a.out 75776-kspace 320 136
	./a.out 76800-kspace 512 308
