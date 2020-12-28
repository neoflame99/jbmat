#include "jbBmp.h"
#include <stdio.h>


Mat read_bmp(const std::string& fname){


    char *memblock = nullptr;

    std::ifstream file (fname.c_str(), std::ios::in | std::ios::binary );
    if (!file.is_open()) {
        fprintf(stderr,"file \'%s\' is not able to open\n",fname.c_str());
        return Mat();
    }
    BmpFileheader   fileheader;
    BmpV3Infoheader v3Infoheader;
    std::streampos size;
    file.read((char *)&fileheader, sizeof(BmpFileheader));
    // check file signature 'BM'.
    //if (fileheader.bfType != (unsigned short)0x424D){ // 0x42='B', 0x4D='M'
    if (fileheader.bfType != (unsigned short)0x4d42){ // 0x42='B', 0x4D='M'
        fprintf(stderr,"file \'%s\' is not a BMP file\n",fname.c_str());
        return Mat();
    }

    file.seekg(0, std::ios::end);
    size = file.tellg(); // get file size
    file.seekg(0, std::ios::beg);

    memblock = new char [size];
    file.read(memblock, size);
    file.close();
    memcpy(&v3Infoheader.bmpInfoHeader, memblock+sizeof(BmpFileheader), sizeof(BmpInfoheader));
    if(v3Infoheader.bmpInfoHeader.biSize >= sizeof(BmpV3Infoheader)){
        memcpy(&v3Infoheader, memblock+sizeof(BmpFileheader), sizeof(BmpV3Infoheader));
    }

    Mat desIm(DTYP::UCHAR, abs(v3Infoheader.bmpInfoHeader.biHeight), abs(v3Infoheader.bmpInfoHeader.biWidth), 3);

    int32 mat_cur_r, mat_r_dir, biHt, biWd, rowSz, rowPadSz;
    if( v3Infoheader.bmpInfoHeader.biHeight< 0){
        mat_cur_r = 0;
        mat_r_dir = 1;
        biHt      = -v3Infoheader.bmpInfoHeader.biHeight;
    }else{
        mat_cur_r = v3Infoheader.bmpInfoHeader.biHeight-1;
        mat_r_dir = -1;
        biHt      = v3Infoheader.bmpInfoHeader.biHeight;
    }
    biWd  = abs(v3Infoheader.bmpInfoHeader.biWidth);
    rowSz = ((biWd*v3Infoheader.bmpInfoHeader.biBitCount+31)/32)*4;  // row size is to be multiples of 4 bytes.
    rowPadSz = rowSz - biWd*v3Infoheader.bmpInfoHeader.biBitCount/8; // padding bytes

    rgbT *mat_row_pt = nullptr;
    unsigned short bitCount = v3Infoheader.bmpInfoHeader.biBitCount;
    int m, n, k, j, step;
    if (bitCount == 32){
        step = 4;
        if(v3Infoheader.bmpInfoHeader.biCompression== BI_RGB){
            rgbQ *mem = (rgbQ *)(memblock+fileheader.bfOffset);
            for(m=0, k=0; m < biHt; ++m, mat_cur_r += mat_r_dir){
                mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
                for(n=0, j=0 ; n < rowSz; n+=step, ++k, ++j){
                    mat_row_pt[j] = mem[k].rgbt;
                }
            }
        }else if(v3Infoheader.bmpInfoHeader.biCompression== BI_BITFIELDS){
            unsigned int *mem = (unsigned int*)(memblock+fileheader.bfOffset);
            unsigned int bitfield;
            unsigned int dat;
            for(m=0, k=0; m < biHt; ++m, mat_cur_r += mat_r_dir){
                mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
                for(n=0, j=0 ; n < rowSz; n+=step, ++k, ++j){
                    dat = mem[k] & v3Infoheader.biBluBitField;
                    for(bitfield = v3Infoheader.biBluBitField; (bitfield & 0x01) == 0; bitfield = bitfield >> 1){
                        dat = dat >> 1;
                    }
                    mat_row_pt[j].blu = (uchar) dat;
                    dat = mem[k] & v3Infoheader.biGrnBitField;
                    for(bitfield = v3Infoheader.biGrnBitField; (bitfield & 0x01) == 0; bitfield = bitfield >> 1){
                        dat = dat >> 1;
                    }
                    mat_row_pt[j].grn = (uchar) dat;
                    dat = mem[k] & v3Infoheader.biRedBitField;
                    for(bitfield = v3Infoheader.biRedBitField; (bitfield & 0x01) == 0; bitfield = bitfield >> 1){
                        dat = dat >> 1;
                    }
                    mat_row_pt[j].red = (uchar) dat;
                }
            }
        }
    }else if (bitCount == 24){
        step = 3;
        rgbT *mem = (rgbT *)(memblock+fileheader.bfOffset);
        for(m=0, k=0 ; m < biHt; ++m, mat_cur_r += mat_r_dir){
            mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
            for(n=0, j=0; n < rowSz-rowPadSz; n+=step, ++k, ++j){
                mat_row_pt[j] = mem[k];
            }
            mem = (rgbT *)((char *)mem+rowPadSz);
        }
    }else if (bitCount == 16){
        step = 2;
        unsigned short* mem = (unsigned short *)(memblock+fileheader.bfOffset);
        for(m=0, k=0 ; m < biHt; ++m, mat_cur_r += mat_r_dir){
            mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
            for(n=0, j=0 ; n < rowSz-rowPadSz; n+=step, ++k, ++j){
                mat_row_pt[j].blu = (mem[k] & 0x001f) << 2;
                mat_row_pt[j].grn = (mem[k] & 0x03e0) >> 3;
                mat_row_pt[j].red = (mem[k] & 0x7c00) >> 8;
            }
            mem = (unsigned short *)mem+rowPadSz;
        }
    }else if (bitCount == 8 ){
        step = 1;
        rgbQ *palette_ptr = (rgbQ *)(memblock+sizeof(BmpFileheader)+v3Infoheader.bmpInfoHeader.biSize);
        uchar* mem = (uchar *)(memblock+fileheader.bfOffset);
        for(m=0, k=0 ; m < biHt; ++m, mat_cur_r += mat_r_dir){
            mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
            for(n=0, j=0 ; n < rowSz-rowPadSz; n+=step, ++k, ++j){
                mat_row_pt[j] = palette_ptr[mem[k]].rgbt;
            }
            mem = (uchar *)mem+rowPadSz;
        }
    }else if (bitCount == 4 ){
        step = 1;
        int mat_step = 2;
        uchar bitDat;
        rgbQ *palette_ptr = (rgbQ *)(memblock+sizeof(BmpFileheader)+v3Infoheader.bmpInfoHeader.biSize);
        uchar* mem = (uchar *)(memblock+fileheader.bfOffset);
        for(m=0, k=0 ; m < biHt; ++m, mat_cur_r += mat_r_dir){
            mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
            for(n=0, j=0 ; n < rowSz-rowPadSz; n+=step, ++k, j+=mat_step){
                bitDat = mem[k] & 0xf0;
                mat_row_pt[j  ] = palette_ptr[bitDat].rgbt;
                bitDat = mem[k] & 0x0f;
                mat_row_pt[j+1] = palette_ptr[bitDat].rgbt;
            }
            mem = (uchar *)mem+rowPadSz;
        }
    }else if (bitCount == 2 ){
        step = 1;
        int mat_step = 4;
        uchar bitDat;
        rgbQ *palette_ptr = (rgbQ *)(memblock+sizeof(BmpFileheader)+v3Infoheader.bmpInfoHeader.biSize);
        uchar* mem = (uchar *)(memblock+fileheader.bfOffset);
        for(m=0, k=0 ; m < biHt; ++m, mat_cur_r += mat_r_dir){
            mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
            for(n=0, j=0; n < rowSz-rowPadSz; n+=step, ++k, j+=mat_step){
                bitDat = mem[k] & 0xc0;
                mat_row_pt[j  ] = palette_ptr[bitDat].rgbt;
                bitDat = mem[k] & 0x30;
                mat_row_pt[j+1] = palette_ptr[bitDat].rgbt;
                bitDat = mem[k] & 0x0c;
                mat_row_pt[j+2] = palette_ptr[bitDat].rgbt;
                bitDat = mem[k] & 0x03;
                mat_row_pt[j+3] = palette_ptr[bitDat].rgbt;
            }
            mem = (uchar *)mem+rowPadSz;
        }
    }else if (bitCount == 1 ){
        step = 1;
        int mat_step = 8;
        uchar bitDat;
        rgbQ *palette_ptr = (rgbQ *)(memblock+sizeof(BmpFileheader)+v3Infoheader.bmpInfoHeader.biSize);
        uchar* mem = (uchar *)(memblock+fileheader.bfOffset);
        for(m=0, k=0 ; m < biHt; ++m, mat_cur_r += mat_r_dir){
            mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
            for(n=0, j=0 ; n < rowSz-rowPadSz; n+=step, ++k, j+=mat_step){
                for(int i=0; i < mat_step; ++i){
                    bitDat = mem[k] & (0x80 >> i);
                    mat_row_pt[j+i] = palette_ptr[bitDat].rgbt;
                }
            }
            mem = (uchar *)mem+rowPadSz;
        }
    }

    delete[] memblock;
    return desIm;
}
