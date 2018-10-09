#ifndef JBIMGPROC_H
#define JBIMGPROC_H
#include "jbMat.h"


class jbImgproc {
public:
    jbImgproc();

    static jbMat rgb2ycc(const jbMat& rgbIm, const int sel_eq = 0);
    static jbMat ycc2rgb(const jbMat& yccIm, const int sel_eq = 0);
    static jbMat rgb2gray(const jbMat& rgbIm, const int HowToGray = 0);
    static jbMat histoPmf(const jbMat& src);
    static jbMat histoCmf(const jbMat& src);
    static jbMat clip_HistoPmf(const jbMat& src, const int clipVal);
    static jbMat clip_HistoCmf(const jbMat& src, const int clipVal);
};

#endif // JBIMGPROC_H
