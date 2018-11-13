#include "qimmat.h"
#include <stdio.h>

namespace qimmat {

    jbMat qim2jbmat( QImage& src, DTYP matDtype){
        QImage::Format fmt = src.format();
        int ch = 0;
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
            return jbMat();
        }
        uint row = src.height();
        uint col = src.width();

        jbMat mat(matDtype, row, col, ch);

        switch(matDtype){
        case DTYP::DOUBLE : _datqim2jbmat<double>(mat, src); break;
        case DTYP::FLOAT  : _datqim2jbmat<float >(mat, src); break;
        case DTYP::INT    : _datqim2jbmat<int   >(mat, src); break;
        case DTYP::UCHAR  : _datqim2jbmat<uchar >(mat, src); break;
        }

        return mat;
    }

    QImage jbmat2qim(const jbMat& src){
    QImage::Format fmt;
    if(src.getChannel()==1)
        fmt = QImage::Format_Grayscale8;
    else if(src.getChannel()==3 || src.getChannel()==4)
        fmt = QImage::Format_RGB32;
    else{
        fprintf(stdout,"jbmat2qim is not supporting format\n");
        return QImage();
    }

    int col = src.getCol();
    int row = src.getRow();
    QImage qim(col, row, fmt);
    qim.fill(0);

    switch(src.getDatType()){
    case DTYP::DOUBLE : _datjbmat2qim<double>(qim, src); break;
    case DTYP::FLOAT  : _datjbmat2qim<float >(qim, src); break;
    case DTYP::INT    : _datjbmat2qim<int   >(qim, src); break;
    case DTYP::UCHAR  : _datjbmat2qim<uchar >(qim, src); break;
    }
    return qim;
}

} // qimmat namespace
