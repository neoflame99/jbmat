#ifndef JBBMP_H
#define JBBMP_H

#include <iostream>
#include <fstream>
#include <string>
#include "types.h"
#include "jbMat.h"


using namespace jmat;
enum BI_CMPRESS {bi_rgb=0, bi_bitfields=3};

#pragma pack (push, 1)               // 구조체를 1바이트 크기로 정렬
typedef struct _BitmapFileHeader     // BMP file header
{
    unsigned short bfType;           // BMP file id
    unsigned int   bfSize;           // file size
    unsigned short bfReserved1;      //
    unsigned short bfReserved2;      //
    unsigned int   bfOffset;         // bitmap data starting point
} BmpFileheader;

typedef struct _BitmapInfoHeader     // BMP info header, basic version
{
    unsigned int   biSize;           // size of info header
    int            biWidth;          // width of bitmap image
    int            biHeight;         // height of bitmap image
    unsigned short biPlanes;         // number of color planes
    unsigned short biBitCount;       // bits per pixel
    unsigned int   biCompression;    // compression type
    unsigned int   biSizeImage;      // storage size of image(including padding)
    int            biXPelsPerMeter;  // resolution of horizontal(meter per pixel)
    int            biYPelsPerMeter;  // resolution of vertical(meter per pixel)
    unsigned int   biClrUsed;        // number of colors in color palette
    unsigned int   biClrImportant;   // number of color indices for bitmap
} BmpInfoheader;

typedef struct _BitmapV3InfoHeader     // info header version 3
{
    BmpInfoheader  bmpInfoHeader;
    unsigned int   biRedBitField;
    unsigned int   biGrnBitField;
    unsigned int   biBluBitField;
    unsigned int   biAlpBitField;
} BmpV3Infoheader;

typedef struct _rgbT      // 24 bits per pixel, 24 color depth
{
    uchar blu;
    uchar grn;
    uchar red;
} rgbT;
typedef struct _rgbQ      // 32 bits per pixel, 32 color depth
{
    rgbT  rgbt;
    uchar alp;
} rgbQ;
#pragma pack(pop)


Mat read_bmp(const std::string& fname);
bool write_bmp(const std::string& fname, const Mat& Im, const int biBits=24);

#endif // JBBMP_H
