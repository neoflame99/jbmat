#include "qimmat.h"
#include <stdio.h>

namespace qimmat {

    Mat qim2mat( QImage& src, DTYP matDtype){
        QImage::Format fmt = src.format();
        uint32 ch = 0;
        if(fmt == QImage::Format_Grayscale8 ){
            ch = 1;
            fprintf(stdout,"image format: grayscale\n");
        }else if(fmt == QImage::Format_RGB888){
            ch = 3;
            fprintf(stdout,"image format: RGB888\n");
        }else if(fmt == QImage::Format_ARGB32 || fmt == QImage::Format_RGB32){
            ch = 3;
            fprintf(stdout,"image format: RGB32\n");
        }else{
            fprintf(stdout,"qim2jbmat is not supporting format\n");
            return Mat();
        }
        uint32 row = static_cast<uint32>(src.height());
        uint32 col = static_cast<uint32>(src.width());

        Mat mat(matDtype, row, col, ch);

        switch(matDtype){
        case DTYP::DOUBLE : _datqim2mat<double>(mat, src); break;
        case DTYP::FLOAT  : _datqim2mat<float >(mat, src); break;
        case DTYP::INT    : _datqim2mat<int32 >(mat, src); break;
        case DTYP::UCHAR  : _datqim2mat<uchar >(mat, src); break;
        case DTYP::CMPLX  : _datqim2mat<cmplx >(mat, src);
        }

        return mat;
    }

    QImage mat2qim(const Mat& src){
    QImage::Format fmt;
    if(src.getChannel()==1)
        fmt = QImage::Format_Grayscale8;
    else if(src.getChannel()==3 || src.getChannel()==4)
        fmt = QImage::Format_RGB32;
    else{
        fprintf(stdout,"jbmat2qim is not supporting format\n");
        return QImage();
    }

    int32 col = static_cast<int32>(src.getCol());
    int32 row = static_cast<int32>(src.getRow());
    QImage qim(col, row, fmt);
    qim.fill(0);

    switch(src.getDatType()){
    case DTYP::DOUBLE : _datmat2qim<double>(qim, src); break;
    case DTYP::FLOAT  : _datmat2qim<float >(qim, src); break;
    case DTYP::INT    : _datmat2qim<int32 >(qim, src); break;
    case DTYP::UCHAR  : _datmat2qim<uchar >(qim, src); break;
    case DTYP::CMPLX  : _datmat2qim<cmplx >(qim, src);
    }
    return qim;
}

} // qimmat namespace
