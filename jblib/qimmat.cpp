#include "qimmat.h"
#include <stdio.h>

QimMat::QimMat()
{
}

jbMat QimMat::qim2jbmat(const QImage& src){
    QImage::Format fmt = src.format();
    int ch;
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
    int row = src.height();
    int col = src.width();
    int len = row*col*ch;
    jbMat mat(row, col, ch);

    const unsigned char* qim_dat = src.bits();
    double *mat_dat = mat.getMat().get();
    int y, k;
    if( fmt == QImage::Format_ARGB32 || fmt == QImage::Format_RGB32){
        for( y=0,k=0; y < len; y+=3, k+=4){
            mat_dat[y  ] = qim_dat[k+2]; // r
            mat_dat[y+1] = qim_dat[k+1]; // g
            mat_dat[y+2] = qim_dat[k  ]; // b
        }
    }else if( fmt == QImage::Format_RGB888){
        int rowcol = row*col ;
        int rowcol2= rowcol*2;
        for( y=0, k=0; y < rowcol; y++, k+=3){
            mat(y        ) = qim_dat[k  ];
            mat(y+rowcol ) = qim_dat[k+1];
            mat(y+rowcol2) = qim_dat[k+2];
        }
    }else {
        for( y=0; y < len; y++)
            mat[y] = qim_dat[y];
    }

    return mat;
}

QImage QimMat::jbmat2qim(const jbMat& src){
    QImage::Format fmt;
    if(src.getChannel()==1)
        fmt = QImage::Format_Grayscale8;
    else if(src.getChannel()==3 || src.getChannel()==4)
        fmt = QImage::Format_RGB32;
    else{
        fprintf(stdout,"jbmat2qim is not supporting format\n");
        return QImage();
    }

    QImage qim(src.getCol(), src.getRow(), fmt);
    qim.fill(0);
    unsigned char* qim_dat = qim.bits();
    double* mat_dat = src.getMat().get();
    int y, k;
    int len = src.getLength();

    double a,b,c;

    if(fmt == QImage::Format_RGB32){
        for( y=0, k=0 ; y < len; y+=3, k+=4){
            a = mat_dat[y  ]; // r
            b = mat_dat[y+1]; // g
            c = mat_dat[y+2]; // b
            if(a > 255) a = 255;
            else if( a < 0) a = 0;
            if(b > 255) b = 255;
            else if( b < 0) b = 0;
            if(c > 255) c = 255;
            else if( c < 0) c = 0;

            qim_dat[k  ] = static_cast<unsigned char>(c);
            qim_dat[k+1] = static_cast<unsigned char>(b);
            qim_dat[k+2] = static_cast<unsigned char>(a);
            qim_dat[k+3] = 255;
        }
    }else if(fmt==QImage::Format_RGB888){
        for( y=0 ; y < len; y+=3){
            a = mat_dat[y  ];
            b = mat_dat[y+1];
            c = mat_dat[y+2];
            if(a > 255) a = 255;
            else if( a < 0) a = 0;
            if(b > 255) b = 255;
            else if( b < 0) b = 0;
            if(c > 255) c = 255;
            else if( c < 0) c = 0;

            qim_dat[y  ] = static_cast<unsigned char>(a);
            qim_dat[y+1] = static_cast<unsigned char>(b);
            qim_dat[y+2] = static_cast<unsigned char>(c);
        }
    }else if(fmt==QImage::Format_Grayscale8){
        for( y=0 ; y < len ; y++){
            a = mat_dat[y];
            if(a > 255) a = 255;
            else if( a < 0) a = 0;
            qim_dat[y] = static_cast<unsigned char>(a);
        }
    }

    return qim;
}
