#include "jbBmp.h"
#include <stdio.h>


Mat read_bmp(const std::string& fname){


    char *memblock = nullptr;

    std::ifstream file (fname.c_str(), std::ios::in | std::ios::binary );
    if (!file.is_open()) {
        fprintf(stderr,"file \'%s\' is not able to open\n",fname.c_str());
        return Mat();
    }
    BmpFileheader fileheader;
    BmpInfoheader Infoheader;
    std::streampos size;
    file.read((char *)&fileheader, sizeof(BmpFileheader));
    // check file signature 'BM'.
    if (fileheader.bfType != (unsigned short)0x424D){ // 0x42='B', 0x4D='M'
        fprintf(stderr,"file \'%s\' is not a BMP file\n",fname.c_str());
        return Mat();
    }

    file.seekg(0, std::ios::end);
    size = file.tellg(); // get file size
    file.seekg(0, std::ios::beg);

    memblock = new char [size];
    file.read(memblock, size);
    file.close();
    memcpy(&Infoheader, memblock+sizeof(BmpFileheader), sizeof(BmpInfoheader));

    Mat desIm(DTYP::UCHAR, abs(Infoheader.biHeight), abs(Infoheader.biWidth), 3);

    int32 mat_cur_r, mat_r_dir, biHt, biWd, rowSz, rowPadSz;
    if( Infoheader.biHeight< 0){
        mat_cur_r = 0;
        mat_r_dir = 1;
        biHt      = -Infoheader.biHeight;
    }else{
        mat_cur_r = Infoheader.biHeight-1;
        mat_r_dir = -1;
        biHt      = Infoheader.biHeight;
    }
    biWd  = abs(Infoheader.biWidth);
    rowSz = ((biWd*Infoheader.biBitCount+31)/32)*4;  // row size is to be multiples of 4 bytes.
    rowPadSz = rowSz - biWd*Infoheader.biBitCount/8; // padding bytes

    rgbT *mat_row_pt = nullptr;

    int m, n, k, j, step;
    if (Infoheader.biBitCount == 32){
        step = 4;
        rgbQ *mem = (rgbQ *)(memblock+sizeof(BmpFileheader)+fileheader.bfOffBits);
        for(m=0, k=0, j=0; m < biHt; ++m, mat_cur_r += mat_r_dir){
            mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
            for(n=0 ; n < rowSz; n+=step, ++k, ++j){
                mat_row_pt[j].blu = mem[k].blu;
                mat_row_pt[j].grn = mem[k].grn;
                mat_row_pt[j].red = mem[k].red;
            }
        }
    }else if (Infoheader.biBitCount == 24){
        step = 3;
        rgbT *mem = (rgbT *)(memblock+sizeof(BmpFileheader)+fileheader.bfOffBits);
        for(m=0, k=0, j=0 ; m < biHt; ++m, mat_cur_r += mat_r_dir){
            mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
            for(n=0 ; n < rowSz-rowPadSz; n+=step, ++k, ++j){
                mat_row_pt[j] = mem[k];
            }
            mem = (rgbT *)((char *)mem+rowPadSz);
        }
    }else if (Infoheader.biBitCount == 16){

    }else if (Infoheader.biBitCount <= 8 ){
        //char *palette_ptr = memblock+sizeof(BmpFileheader)+Infoheader.biSize;
    }

    delete[] memblock;
    return desIm;
}
