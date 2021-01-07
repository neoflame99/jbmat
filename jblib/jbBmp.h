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
typedef struct _BitmapFileHeader     // BMP 비트맵 파일 헤더 구조체
{
    unsigned short bfType;           // BMP 파일 매직 넘버
    unsigned int   bfSize;           // 파일 크기
    unsigned short bfReserved1;      // 예약
    unsigned short bfReserved2;      // 예약
    unsigned int   bfOffset;        // 비트맵 데이터의 시작 위치
} BmpFileheader;

typedef struct _BitmapInfoHeader     // BMP 비트맵 정보 헤더 구조체(DIB 헤더)
{
    unsigned int   biSize;           // 현재 구조체의 크기
    int            biWidth;          // 비트맵 이미지의 가로 크기
    int            biHeight;         // 비트맵 이미지의 세로 크기
    unsigned short biPlanes;         // 사용하는 색상판의 수
    unsigned short biBitCount;       // 픽셀 하나를 표현하는 비트 수
    unsigned int   biCompression;    // 압축 방식
    unsigned int   biSizeImage;      // 비트맵 이미지의 픽셀 데이터 크기(including padding)
    int            biXPelsPerMeter;  // 그림의 가로 해상도(미터당 픽셀)
    int            biYPelsPerMeter;  // 그림의 세로 해상도(미터당 픽셀)
    unsigned int   biClrUsed;        // 색상 테이블에서 실제 사용되는 색상 수
    unsigned int   biClrImportant;   // 비트맵을 표현하기 위해 필요한 색상 인덱스 수
} BmpInfoheader;

typedef struct _BitmapV3InfoHeader     // BMP 비트맵 정보 헤더 구조체(DIB 헤더)
{
    BmpInfoheader  bmpInfoHeader;
    unsigned int   biRedBitField;
    unsigned int   biGrnBitField;
    unsigned int   biBluBitField;
    unsigned int   biAlpBitField;
} BmpV3Infoheader;

typedef struct _rgbT      // 24비트 비트맵 이미지의 픽셀 구조체
{
    uchar blu;
    uchar grn;
    uchar red;
} rgbT;
typedef struct _rgbQ      // 32비트 비트맵 이미지의 픽셀 구조체
{
    rgbT  rgbt;
    uchar alp;
} rgbQ;
#pragma pack(pop)


Mat read_bmp(const std::string& fname);
bool write_bmp(const std::string& fname, const Mat& Im, const int biBits=24);

#endif // JBBMP_H