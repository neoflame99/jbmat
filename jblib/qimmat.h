#ifndef QIMMAT_H
#define QIMMAT_H

#include <QImage>
#include <QFile>
#include "types.h"
#include "jbMat.h"
#include "satcast.h"


using namespace jmat;

namespace qimmat {
    Mat  qim2mat( QImage& src, const DTYP matDtype=DTYP::DOUBLE);
    QImage mat2qim(const Mat& src);

    Mat read_hdr(QString& filename, bool normalizing = false);

    template <typename _T> void _datqim2mat(Mat& tar, QImage& src);
    template <typename _T> void _datmat2qim(QImage& tar, const Mat& src);




    template <typename _T> void _datqim2mat(Mat& tar, QImage& src){
        if(tar.isEmpty()) return ;

        uchar* qim_dat = src.bits();
        _T*    mat_dat = tar.getDataPtr<_T>();

        uint32 row = static_cast<uint32>(src.height());
        uint32 col = static_cast<uint32>(src.width());
        uint32 lenRowCol  = row*col;
        uint32 lenRowCol2 = lenRowCol << 1;

        uint32 y, k;
        QImage::Format fmt =  src.format();
        if( fmt == QImage::Format_ARGB32 || fmt == QImage::Format_RGB32){
            for( y=0,k=0; y < lenRowCol; y++, k+=4){
                mat_dat[y           ] = qim_dat[k+2]; // r
                mat_dat[y+lenRowCol ] = qim_dat[k+1]; // g
                mat_dat[y+lenRowCol2] = qim_dat[k  ]; // b
            }
        }else if( fmt == QImage::Format_RGB888){
            for( y=0, k=0; y < lenRowCol; y++, k+=3){
                mat_dat[y           ] = qim_dat[k+2];  // r
                mat_dat[y+lenRowCol ] = qim_dat[k+1];  // g
                mat_dat[y+lenRowCol2] = qim_dat[k  ];  // b
            }
        }else {
            for( y=0; y < lenRowCol; y++)
                mat_dat[y] = qim_dat[y];
        }
    }

    template <typename _T> void _datmat2qim(QImage& tar, const Mat& src){
        if(src.isEmpty()) return;

        uint32 col = src.getCol();
        uint32 row = src.getRow();

        QImage::Format fmt = tar.format();

        uchar* qim_dat = tar.bits();
        _T*    mat_dat = src.getDataPtr<_T>();
        //_T a,b,c;
        uint32 y, k;
        uint32 lenRowCol  = row*col;
        uint32 lenRowCol2 = lenRowCol << 1;

        if(fmt == QImage::Format_RGB32){
            for( y=0, k=0 ; y < lenRowCol; y++, k+=4){
                //a = mat_dat[y           ]; // r
                //b = mat_dat[y+lenRowCol ]; // g
                //c = mat_dat[y+lenRowCol2]; // b
                //if(a > 255) a = 255; else if( a < 0) a = 0;
                //if(b > 255) b = 255; else if( b < 0) b = 0;
                //if(c > 255) c = 255; else if( c < 0) c = 0;
                //qim_dat[k  ] = static_cast<unsigned char>(c);
                //qim_dat[k+1] = static_cast<unsigned char>(b);
                //qim_dat[k+2] = static_cast<unsigned char>(a);
                //qim_dat[k+3] = 255;
                qim_dat[k  ] = sat_cast<uchar>(mat_dat[y+lenRowCol2]);
                qim_dat[k+1] = sat_cast<uchar>(mat_dat[y+lenRowCol ]);
                qim_dat[k+2] = sat_cast<uchar>(mat_dat[y           ]);
                qim_dat[k+3] = 255;
            }
        }else if(fmt==QImage::Format_RGB888){
            for( y=0, k=0 ; y < lenRowCol; y++, k+=3 ){
                //a = mat_dat[y           ]; // r
                //b = mat_dat[y+lenRowCol ]; // g
                //c = mat_dat[y+lenRowCol2]; // b
                //if(a > 255) a = 255; else if( a < 0) a = 0;
                //if(b > 255) b = 255; else if( b < 0) b = 0;
                //if(c > 255) c = 255; else if( c < 0) c = 0;
                //qim_dat[k  ] = static_cast<unsigned char>(c);
                //qim_dat[k+1] = static_cast<unsigned char>(b);
                //qim_dat[k+2] = static_cast<unsigned char>(a);
                qim_dat[k  ] = sat_cast<uchar>(mat_dat[y+lenRowCol2]);
                qim_dat[k+1] = sat_cast<uchar>(mat_dat[y+lenRowCol ]);
                qim_dat[k+2] = sat_cast<uchar>(mat_dat[y           ]);
            }
        }else if(fmt==QImage::Format_Grayscale8){
            for( y=0 ; y < lenRowCol ; y++){
                //a = mat_dat[y];
                //if(a > 255) a = 255; else if( a < 0) a = 0;
                //qim_dat[y] = static_cast<unsigned char>(a);
                qim_dat[y] = sat_cast<uchar>(mat_dat[y]);
            }
        }
    }

    template <> inline void _datmat2qim<cmplx>(QImage& tar, const Mat& src){
        if(src.isEmpty()) return;

        uint32 col = src.getCol();
        uint32 row = src.getRow();

        QImage::Format fmt = tar.format();

        uchar* qim_dat = tar.bits();
        cmplx* mat_dat = src.getDataPtr<cmplx>();
        //cmplx a,b,c;
        uint32 y, k;
        uint32 lenRowCol  = row*col;
        uint32 lenRowCol2 = lenRowCol << 1;

        if(fmt == QImage::Format_RGB32){
            for( y=0, k=0 ; y < lenRowCol; y++, k+=4){
                qim_dat[k  ] = sat_cast<uchar>(mat_dat[y+lenRowCol2].re);
                qim_dat[k+1] = sat_cast<uchar>(mat_dat[y+lenRowCol ].re);
                qim_dat[k+2] = sat_cast<uchar>(mat_dat[y           ].re);
                qim_dat[k+3] = 255;
            }
        }else if(fmt==QImage::Format_RGB888){
            for( y=0, k=0 ; y < lenRowCol; y++, k+=3 ){
                qim_dat[k  ] = sat_cast<uchar>(mat_dat[y+lenRowCol2].re);
                qim_dat[k+1] = sat_cast<uchar>(mat_dat[y+lenRowCol ].re);
                qim_dat[k+2] = sat_cast<uchar>(mat_dat[y           ].re);
            }
        }else if(fmt==QImage::Format_Grayscale8){
            for( y=0 ; y < lenRowCol ; y++){
                qim_dat[y] = sat_cast<uchar>(mat_dat[y].re);
            }
        }
    }

} // end of qimmat namespace

#endif // QIMMAT_H
