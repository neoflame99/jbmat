#include "jbBmp.h"
#include <stdio.h>


bool read_bmp(const std::string& fname, Mat& desIm){


    char *memblock = nullptr;

    std::ifstream file (fname.c_str(), std::ios::in | std::ios::binary );
    if (!file.is_open()) {
        fprintf(stderr,"file \'%s\' is not able to open\n",fname.c_str());
        return false;
    }
    BmpFileheader fileheader;
    BmpInfoheader Infoheader;
    std::streampos size;
    file.read((char *)&fileheader, sizeof(BmpFileheader));
    // check file signature 'BM'.
    if (fileheader.bfType != (unsigned short)0x424D){ // 0x42='B', 0x4D='M'
        fprintf(stderr,"file \'%s\' is not a BMP file\n",fname.c_str());
        return false;
    }

    file.seekg(0, std::ios::end);
    size = file.tellg(); // get file size
    file.seekg(0, std::ios::beg);

    memblock = new char [size];
    file.read(memblock, size);
    file.close();
    memcpy(&Infoheader, memblock+sizeof(BmpFileheader), sizeof(BmpInfoheader));


    if (Infoheader.biBitCount == 32){

    }else if (Infoheader.biBitCount == 24){

    }else if (Infoheader.biBitCount == 16){

    }else if (Infoheader.biBitCount <= 8 ){
        char *palette_ptr = memblock+sizeof(BmpFileheader)+Infoheader.biSize;
    }

    delete[] memblock;
    return true;
}
