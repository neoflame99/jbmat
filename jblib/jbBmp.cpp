#include "jbBmp.h"
#include <stdio.h>


Mat read_bmp(const std::string& fname){

    std::ifstream file (fname.c_str(), std::ios::in | std::ios::binary );
    if (!file.is_open()) {
        fprintf(stderr,"file \'%s\' is not able to open\n",fname.c_str());
        return Mat();
    }
    char *memblock = nullptr;
    BmpFileheader   fileheader;
    BmpV3Infoheader v3Infoheader;
    std::streampos  size;
    file.read((char *)&fileheader, sizeof(BmpFileheader));
    // check file signature 'BM'.
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
        if(v3Infoheader.bmpInfoHeader.biCompression== BI_CMPRESS::bi_rgb){
            rgbQ *mem = (rgbQ *)(memblock+fileheader.bfOffset);
            for(m=0, k=0; m < biHt; ++m, mat_cur_r += mat_r_dir){
                mat_row_pt = (rgbT *)desIm.getRowPtr(mat_cur_r);
                for(n=0, j=0 ; n < rowSz; n+=step, ++k, ++j){
                    mat_row_pt[j] = mem[k].rgbt;
                }
            }
        }else if(v3Infoheader.bmpInfoHeader.biCompression== BI_CMPRESS::bi_bitfields){
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

bool write_bmp(const std::string& fname, const Mat& Im, const int biBits){
    // bitsCount: only support rgb32, rgb24, rgb8(gray)
    // Compression: only support BI_RGB
    //

    switch(biBits){
    case 32: break;
    case 24: break;
    case  8: break;
    default:
        fprintf(stderr,"for biBits of write_bmp() func. It can only supports 32, 24 and 8 bits\n");
        return false;
    }
    switch(Im.getChannel()){
    case 4: break;
    case 3: break;
    case 1: break;
    default:
        fprintf(stderr,"for Mat of write_bmp() func. It can only supports 4, 3 and 1 channels\n");
        return false;
    }

    std::ofstream file (fname.c_str(), std::ios::out | std::ios::binary );
    if (!file.is_open()) {
        fprintf(stderr,"file \'%s\' is not able to open\n",fname.c_str());
        return false;
    }
    char *memblock = nullptr;
    BmpFileheader   fileheader;
    BmpInfoheader   infoheader;

    uint32 palette_sz ;
    size_t bmpsize;
    size_t data_sz;
    int32  mat_cur_r, biHt, rowSz, rowPadSz;

    rowSz      = ((Im.getCol()*biBits+31)/32)*4;  // row size is to be multiples of 4 bytes.
    rowPadSz   = rowSz - Im.getCol()*biBits/8; // padding bytes
    data_sz    = rowSz * Im.getRow();
    palette_sz = biBits > 8 ? 0 : 0x00000004 << biBits;
    bmpsize    = sizeof(BmpFileheader)+sizeof(BmpInfoheader)+ palette_sz + data_sz;
    memblock   = (char*)calloc(bmpsize, sizeof(char));

    // bmp file header
    fileheader.bfType       = (unsigned short)0x4d42; // 0x42='B', 0x4D='M'
    fileheader.bfSize       = bmpsize;
    fileheader.bfReserved1  = 0;
    fileheader.bfReserved2  = 0;
    fileheader.bfOffset     = sizeof(BmpFileheader)+sizeof(BmpInfoheader)+palette_sz ;

    // bmp info header
    infoheader.biSize         = sizeof(BmpInfoheader);
    infoheader.biWidth        = Im.getCol();
    infoheader.biHeight       = Im.getRow();
    infoheader.biPlanes       = 1;
    infoheader.biBitCount     = biBits;
    infoheader.biCompression  = BI_CMPRESS::bi_rgb;
    infoheader.biSizeImage    = data_sz;
    infoheader.biXPelsPerMeter= 0;
    infoheader.biYPelsPerMeter= 0;
    infoheader.biClrUsed      = palette_sz >> 2;
    infoheader.biClrImportant = 0;


    // writing bmp file
    memcpy(memblock, &fileheader, sizeof(BmpFileheader));
    memcpy(memblock+sizeof(BmpFileheader), &infoheader, sizeof(BmpInfoheader));

    mat_cur_r = Im.getRow()-1;
    biHt      = Im.getRow();

    int m, n, k, j, step;
    if (biBits == 32){
        step = 4;
        rgbQ *mem = (rgbQ *)(memblock+fileheader.bfOffset);
        if( Im.getChannel()==4){
            rgbQ *mat_row_pt = nullptr;
            for(m=0, k=0; m < biHt; ++m, mat_cur_r-- ){
                mat_row_pt = (rgbQ *)Im.getRowPtr(mat_cur_r);
                for(n=0, j=0 ; n < rowSz; n+=step, ++k, ++j){
                    mem[k] = mat_row_pt[j];
                }
            }
        }else if(Im.getChannel()==3){
            rgbT *mat_row_pt = nullptr;
            for(m=0, k=0; m < biHt; ++m, mat_cur_r-- ){
                mat_row_pt = (rgbT *)Im.getRowPtr(mat_cur_r);
                for(n=0, j=0 ; n < rowSz; n+=step, ++k, ++j){
                    mem[k].alp = 0xff;
                    mem[k].rgbt = mat_row_pt[j];
                }
            }
        }
    }else if (biBits == 24){
        step = 3;
        rgbT *mem = (rgbT *)(memblock+fileheader.bfOffset);
        if( Im.getChannel()==4){
            rgbQ *mat_row_pt = nullptr;
            for(m=0, k=0; m < biHt; ++m, mat_cur_r-- ){
                mat_row_pt = (rgbQ *)Im.getRowPtr(mat_cur_r);
                for(n=0, j=0 ; n < rowSz; n+=step, ++k, ++j){
                    mem[k] = mat_row_pt[j].rgbt;
                }
                mem = (rgbT *)((char *)mem+rowPadSz);
            }
        }else if(Im.getChannel()==3){
            rgbT *mat_row_pt = nullptr;
            for(m=0, k=0; m < biHt; ++m, mat_cur_r-- ){
                mat_row_pt = (rgbT *)Im.getRowPtr(mat_cur_r);
                for(n=0, j=0 ; n < rowSz; n+=step, ++k, ++j){
                    mem[k] = mat_row_pt[j];
                }
                mem = (rgbT *)((char *)mem+rowPadSz);
            }
        }
    }else if (biBits == 8 ){ // gray8 only supported
        step = 1;
        uchar *palette_ptr = (uchar *)(memblock+sizeof(BmpFileheader)+infoheader.biSize);
        uchar* mem = (uchar *)(memblock+fileheader.bfOffset);
        for(uint32 mm=0, mn=0; mm < palette_sz; mm+=4, mn++){
            palette_ptr[mm  ] = mn;
            palette_ptr[mm+1] = mn;
            palette_ptr[mm+2] = mn;
            palette_ptr[mm+3] = 00;
        }
        if( Im.getChannel()==4){
            rgbQ *mat_row_pt = nullptr;
            for(m=0, k=0 ; m < biHt; ++m, mat_cur_r--){
                mat_row_pt = (rgbQ *)Im.getRowPtr(mat_cur_r);
                for(n=0, j=0 ; n < rowSz-rowPadSz; n+=step, ++k, ++j){
                    mem[k] = (uchar)(mat_row_pt[j].rgbt.red*0.299+mat_row_pt[j].rgbt.grn*0.587+mat_row_pt[j].rgbt.blu*0.114);
                }
                mem = (uchar *)mem+rowPadSz;
            }
        }else if( Im.getChannel()==3){
            rgbT *mat_row_pt = nullptr;
            for(m=0, k=0 ; m < biHt; ++m, mat_cur_r--){
                mat_row_pt = (rgbT *)Im.getRowPtr(mat_cur_r);
                for(n=0, j=0 ; n < rowSz-rowPadSz; n+=step, ++k, ++j){
                    mem[k] = (uchar)(mat_row_pt[j].red*0.299+mat_row_pt[j].grn*0.587+mat_row_pt[j].blu*0.114);
                }
                mem = (uchar *)mem+rowPadSz;
            }
        }else if( Im.getChannel()==1){
            uchar *mat_row_pt = nullptr;
            for(m=0, k=0 ; m < biHt; ++m, mat_cur_r--){
                mat_row_pt = Im.getRowPtr(mat_cur_r);
                for(n=0, j=0 ; n < rowSz-rowPadSz; n+=step, ++k, ++j){
                    mem[k] = mat_row_pt[j];
                }
                mem = (uchar *)mem+rowPadSz;
            }
        }
    }

    file.write(memblock, bmpsize);
    file.close();

    free(memblock);
    return true;
}
