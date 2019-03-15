#include "qimmat.h"
#include <stdio.h>

namespace qimmat {

    Mat qim2mat( QImage& src, DTYP matDtype){
        QImage::Format fmt = src.format();
        uint32 ch = 0;
        if(fmt == QImage::Format_Grayscale8 ){
            ch = 1;
            fprintf(stderr,"image format: grayscale\n");
        }else if(fmt == QImage::Format_RGB888){
            ch = 3;
            fprintf(stderr,"image format: RGB888\n");
        }else if(fmt == QImage::Format_ARGB32 || fmt == QImage::Format_RGB32){
            ch = 3;
            fprintf(stderr,"image format: RGB32\n");
        }else{
            fprintf(stderr,"qim2jbmat is not supporting format\n");
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
            fprintf(stderr,"jbmat2qim is not supporting format\n");
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

    Mat read_hdr(QFile& file, bool normalizing ){
        int32   ysize       = 0;
        int32   xsize       = 0;
        int32   rgbe_format = 0;
        bool    isHdrfile   = false;
        uchar   *rgbe       = nullptr;
        double  exposure    = 1.0;
        int32   colorformat = 0;     // 0 : rgb, 1: xyz, 2: Yxy ...

        int32   ypos, xpos;
        int32   t;
        bool    isOk;
        QStringList strlist;

        if( !file.open(QIODevice::ReadOnly | QIODevice::Text) ){
            fprintf(stderr, "Not open file!");
            return Mat();
        }
        while( !file.atEnd()) {
            QString line = file.readLine().trimmed();
            if(line.contains(QString("#?RADIANCE"),Qt::CaseInsensitive)){
                isHdrfile = true;
            }
            else if(line.contains(QString("format"),Qt::CaseInsensitive) ){
                t = line.indexOf(QString("="));
                if ( line.mid(t+1).compare(QString("32-bit_rle_rgbe"),Qt::CaseInsensitive) ==0){
                    rgbe_format = 0;
                    colorformat = 0;
                }else if( line.mid(t+1).compare(QString("32-bit_rle_xyze"),Qt::CaseInsensitive)==0){
                    rgbe_format = 1;
                    colorformat = 1;
                }
            }
            else if(line.contains(QString("exposure"),Qt::CaseInsensitive) ){
                t = line.indexOf(QString("="));
                exposure = line.mid(t+1).toDouble(&isOk);
                if( !isOk)
                    exposure = 1.0;
            }
            else if(line.contains(QString("pvalue"),Qt::CaseInsensitive) ){
                t = line.indexOf(QString("+e"));
                strlist = line.mid(t).trimmed().split(QRegExp("\\s"));
                exposure = strlist.at(1).toDouble(&isOk);
                if( !isOk)
                    exposure = 1.0;
            }
            else if( line.isEmpty() ) {
                QString headfl = file.readLine().simplified();
                ypos = headfl.indexOf(QString("Y"), 0 , Qt::CaseInsensitive);
                xpos = headfl.indexOf(QString("X"), 0 , Qt::CaseInsensitive);

                ysize = headfl.mid(ypos+1, xpos-ypos-2).toInt();
                xsize = headfl.mid(xpos+1).toInt();
                break;
            }
        }
        if( !isHdrfile) {
            fprintf(stderr,"The file is not Radiance file!\n");
            return Mat();
        }

        fprintf(stdout, "ysize : %d , xsize : %d \n ", ysize , xsize);
        // read binary from here
        file.setTextModeEnabled(false);
        QByteArray hdr_data = file.readAll();


        int32 dlen = (ysize * xsize) << 2;

        rgbe = new uchar [dlen];
        if (rgbe == nullptr ){
            fprintf(stderr,"Memory allocation error\n");
            return Mat();
        }

        uchar   a,b;
        uchar   c,d;
        uchar   rlc;
        uchar   lc, rdat;
        int32   xcnt = 0, cxcnt = 0;
        int32   ycnt = 0;
        uchar   *pdata = nullptr;
        bool    c0, c1, c2, c3;
        int32   hbyteX = xsize >> 8;
        int32   lbyteX = xsize - (hbyteX << 8);
        int32   y, ccnt;


        pdata = (uchar *)hdr_data.data();

        a = pdata[0];
        b = pdata[1];
        c = pdata[2];
        d = pdata[3];
        if( dlen == hdr_data.length()) { // uncompressed
            for ( y=0; y < dlen; y++)
                rgbe[y] = pdata[y];

        }else if( a==2 && b==2 && c==hbyteX && d==lbyteX) { // check scanline codes for New RLE
            for ( y=0; y < ysize; y++){
                ycnt = (y * xsize);

                a = *pdata++;
                b = *pdata++;
                c = *pdata++;
                d = *pdata++;
                c0 = a == 0x02;
                c1 = b == 0x02;
                c2 = c == hbyteX;
                c3 = d == lbyteX;
                if( c0 && c1 && c2 && c3){  // New RLE
                    // decoding each channel for the scanline
                    for( ccnt = 0; ccnt < 4; ccnt++){
                        xcnt = 0;
                        while(xcnt < xsize){
                            //rlc = (unsigned char)hdr_data.at(idx++);
                            rlc = *pdata++;
                            if(rlc <= 128){  // read literally values
                                for(lc = 0; lc < rlc; lc++ ){
                                    //rdat = (unsigned char)hdr_data.at(idx++);
                                    rdat  = *pdata++;
                                    cxcnt = ((ycnt + xcnt+lc) << 2) + ccnt;
                                    rgbe[cxcnt] = rdat;
                                }
                            }else{ // read run length value
                                rlc = rlc - 128;
                                rdat   = *pdata++;
                                for(lc = 0; lc < rlc; lc++){
                                    cxcnt = ((ycnt + xcnt+lc) << 2) + ccnt;
                                    rgbe[cxcnt] = rdat;
                                }
                            }
                            xcnt += rlc;
                            if(xcnt > xsize)
                                fprintf(stderr," run length coding error : xcnt = %d \n", xcnt);
                        }
                    }
                }

            }
        }else { // Old RLE
            int32 num = 0;
            int32 consecRpt = 0;
            int32 i,k;
            xcnt = 0;
            cxcnt = 0;
            for ( i=0; i < hdr_data.length(); i+=4){
                a = pdata[i  ];
                b = pdata[i+1];
                c = pdata[i+2];
                d = pdata[i+3];
                if( a==1 && b==1 && c==1){
                    num = d << (consecRpt * 8);
                    xcnt = cxcnt-4;  // store the last rgbe point
                    for( k=0; k < num; k++){
                        rgbe[cxcnt++] = rgbe[xcnt  ];
                        rgbe[cxcnt++] = rgbe[xcnt+1];
                        rgbe[cxcnt++] = rgbe[xcnt+2];
                        rgbe[cxcnt++] = rgbe[xcnt+3];
                    }
                    consecRpt++;
                }else{
                    rgbe[cxcnt++] = a;
                    rgbe[cxcnt++] = b;
                    rgbe[cxcnt++] = static_cast<uchar>(c);
                    rgbe[cxcnt++] = static_cast<uchar>(d);
                    consecRpt = 0;
                }
            }

        }

        Mat hdrMat(DTYP::DOUBLE, (uint32)ysize, (uint32)xsize, 3);
        double* dat64f = hdrMat.getDataPtr<double>();

        int32 xc, x;
        double rt, gt, bt, et;
        double max = -1;
        ypos = 0;
        for( y=0; y < ysize; y++){
            ypos  = y * xsize;
            for( x=0; x < xsize; x++){
                cxcnt = (ypos + x) << 2;
                rt = rgbe[cxcnt  ];
                gt = rgbe[cxcnt+1];
                bt = rgbe[cxcnt+2];
                et = rgbe[cxcnt+3];

                if(rgbe[cxcnt+3]==0){
                    rt = 0.0;
                    gt = 0.0;
                    bt = 0.0;
                }else{
                    rt = (rt+0.5) * pow(2,et-128-8) ;
                    gt = (gt+0.5) * pow(2,et-128-8) ;
                    bt = (bt+0.5) * pow(2,et-128-8) ;
                }
                xc = ( ypos + x );
                xc = ( xc << 1) + xc ; // xc = (ypos+x)*3;
                dat64f[xc  ] = rt;
                dat64f[xc+1] = gt;
                dat64f[xc+2] = bt;

                max = (max < rt) ? rt : max;
                max = (max < gt) ? gt : max;
                max = (max < bt) ? bt : max;
            }
        }

        // normalizing
        if(normalizing){
            for( uint32 j = 0; j < hdrMat.getLength(); ++j)
                dat64f[j] /= max;
        }


        return hdrMat;
    }

} // qimmat namespace
