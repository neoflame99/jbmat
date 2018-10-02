#include "jreadhdr.h"
#include <QByteArray>
#include <stdio.h>
#include <QRegExp>
#include <QStringList>
#include <math.h>
double rgb2xyz_bt709[3][3] = { {0.4124564, 0.3575761, 0.1804375},
                               {0.2126729, 0.7151522, 0.0721750},
                               {0.0193339, 0.1191920, 0.9503041} };
double xyz2rgb_bt709[3][3] = { { 3.2404542, -1.5371385, -0.4985314},
                               {-0.9692660,  1.8760108,  0.0415560},
                               { 0.0556434, -0.2040259,  1.0572252} };
void Jreadhdr::init(){
    ysize        = 0;
    xsize        = 0;
    rgbe_format  = 0;
    isHdrfile    = false;
    exposure     = 1.0;
    hasValidData = false;
    rgbe         = 0; // initial value is null
    dat64f       = 0; // inital null
    maxval       = 0.0;
    colorformat  = 0; // 0 : rgb, 1: xyz, 2: Yxy ...
}
Jreadhdr::Jreadhdr()
{
    init();
}
Jreadhdr::Jreadhdr(QFile& file)
{
    init();
    read_hdr(file);
}

Jreadhdr::~Jreadhdr()
{
    if( rgbe != 0)
        delete rgbe;
    if(dat64f != 0)
        delete dat64f;

    rgbe   = 0; // null pointer
    dat64f = 0; // null pointer
}
int Jreadhdr::read_hdr(QFile& file){

    int ypos, xpos;
    int t;
    bool isOk;
    QStringList strlist;
    if( !file.open(QIODevice::ReadOnly | QIODevice::Text) ){
        fprintf(stderr, "Not open file!");
        return -2;
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
                exposure = 1.0f;
        }
        else if(line.contains(QString("pvalue"),Qt::CaseInsensitive) ){
            t = line.indexOf(QString("+e"));
            strlist = line.mid(t).trimmed().split(QRegExp("\\s"));
            exposure = strlist.at(1).toDouble(&isOk);
            if( !isOk)
                exposure = 1.0f;
        }
        else if( line.isEmpty() ) {
            QString headfl = file.readLine().simplified();
            ypos = headfl.indexOf(QString("Y"), 0 , Qt::CaseInsensitive);
            xpos = headfl.indexOf(QString("X"), 0 , Qt::CaseInsensitive);

            ysize = headfl.mid(ypos+1,xpos-ypos-2).toInt();
            xsize = headfl.mid(xpos+1).toInt();
            break;
        }
    }
    if( !isHdrfile) {
        fprintf(stderr,"The file is not Radiance file!\n");
        return -2;
    }

    fprintf(stdout, "ysize : %zd , xsize : %zd \n ", ysize , xsize);
    // read binary from here
    file.setTextModeEnabled(false);
    QByteArray hdr_data = file.readAll();


    size_t dlen = (ysize * xsize) << 2;

    rgbe = new unsigned char [dlen];
    if (rgbe == 0){
        fprintf(stderr,"Memory allocation error\n");
        return -1;
    }


    unsigned char a,b;
    unsigned char c,d;
    unsigned char rlc;
    unsigned char lc, rdat;
    size_t        xcnt = 0, cxcnt = 0;
    size_t        ycnt = 0;
    unsigned char *pdata=0;
    bool          c0, c1, c2, c3;
    unsigned int  hbyteX, lbyteX;
    hbyteX = this->xsize >> 8;
    lbyteX = this->xsize - (hbyteX << 8);
    pdata = (unsigned char*)hdr_data.data();

    a = pdata[0];
    b = pdata[1];
    c = pdata[2];
    d = pdata[3];
    if( dlen == hdr_data.length()) { // uncompressed
        for ( size_t y=0; y < dlen; y++)
            rgbe[y] = pdata[y];

    }else if( a==2 && b==2 && c==hbyteX && d==lbyteX) { // check scanline codes for New RLE
        for (size_t y=0; y < ysize; y++){
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
                for( size_t ccnt = 0; ccnt < 4; ccnt++){
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
                            fprintf(stdout," run length coding error : xcnt = %zd \n", xcnt);
                    }
                }
            }

        }
    }else { // Old RLE
        int num=0;
        int consecRpt = 0;
        xcnt = 0;
        cxcnt = 0;
        for (size_t i=0; i < hdr_data.length(); i+=4){
            a = pdata[i  ];
            b = pdata[i+1];
            c = pdata[i+2];
            d = pdata[i+3];
            if( a==1 && b==1 && c==1){
                num = d << (consecRpt*8);
                xcnt = cxcnt-4;  // store the last rgbe point
                for( size_t k=0; k < num; k++){
                    rgbe[cxcnt++] = rgbe[xcnt  ];
                    rgbe[cxcnt++] = rgbe[xcnt+1];
                    rgbe[cxcnt++] = rgbe[xcnt+2];
                    rgbe[cxcnt++] = rgbe[xcnt+3];
                }
                consecRpt++;
            }else{
                rgbe[cxcnt++] = a;
                rgbe[cxcnt++] = b;
                rgbe[cxcnt++] = (unsigned char)c;
                rgbe[cxcnt++] = (unsigned char)d;
                consecRpt = 0;
            }
        }

    }

    hasValidData = true;
    return 0;
}

int Jreadhdr::conv_rgbe2rgb64f(bool normalizing){

   if( !hasValidData ) return -2;

    dat64f = new double [ysize*xsize*3];
    if(dat64f ==0 ) {
        fprintf(stderr,"In cov_rgbe2dat64f: Memory allocation error for dat64f!\n");
        return -1;
    }


    size_t xc, cxcnt;
    size_t ypos = 0;
    float rt, gt, bt, et;
    double max = -1;
    for(size_t y=0; y < ysize; y++){
        ypos  = y * xsize;
        for(size_t x=0; x < xsize; x++){
            cxcnt = (ypos + x) << 2;
            rt = rgbe[cxcnt  ];
            gt = rgbe[cxcnt+1];
            bt = rgbe[cxcnt+2];
            et = rgbe[cxcnt+3];

            if(rgbe[cxcnt+3]==0){
                rt = 0.0f;
                gt = 0.0f;
                bt = 0.0f;
            }else{
                rt = (rt+0.5) * pow(2,et-128-8) ;
                gt = (gt+0.5) * pow(2,et-128-8) ;
                bt = (bt+0.5) * pow(2,et-128-8) ;
            }
            xc = ( ypos + x );
            xc = (xc << 1) + xc ; // xc = (ypos+x)*3;
            dat64f[xc  ] = rt;
            dat64f[xc+1] = gt;
            dat64f[xc+2] = bt;

            max = (max < rt) ? rt : max;
            max = (max < gt) ? gt : max;
            max = (max < bt) ? bt : max;

        }
    }

    this->maxval = max;
    // normalizing
    if(normalizing){
        for(size_t y=0; y < ysize*xsize*3; y++){
            dat64f[y ] /= max;
        }
    }

    return 0;
}

int Jreadhdr::conv_rgb2xyz(){
    if( colorformat != 0){
        fprintf(stderr," The color format of the image is not RGB\n");
        return -1;
    }
    double x, y, z;
    for(size_t i=0; i < xsize*ysize*3; i+=3){
        x = rgb2xyz_bt709[0][0] * dat64f[i] + rgb2xyz_bt709[0][1] *dat64f[i+1] + rgb2xyz_bt709[0][2] *dat64f[i+3];
        y = rgb2xyz_bt709[1][0] * dat64f[i] + rgb2xyz_bt709[1][1] *dat64f[i+1] + rgb2xyz_bt709[1][2] *dat64f[i+3];
        z = rgb2xyz_bt709[2][0] * dat64f[i] + rgb2xyz_bt709[2][1] *dat64f[i+1] + rgb2xyz_bt709[2][2] *dat64f[i+3];
        dat64f[i]  = x;
        dat64f[i+1]= y;
        dat64f[i+2]= z;
    }
    colorformat = 1;

    return 0;
}

int Jreadhdr::conv_xyz2rgb(){
    if( colorformat != 1){
        fprintf(stderr," The color format of the image is not XYZ\n");
        return -1;
    }
    double r, g, b;
    for(size_t i=0; i < xsize*ysize*3; i+=3){
        r = xyz2rgb_bt709[0][0] * dat64f[i] + xyz2rgb_bt709[0][1] *dat64f[i+1] + xyz2rgb_bt709[0][2] *dat64f[i+3];
        g = xyz2rgb_bt709[1][0] * dat64f[i] + xyz2rgb_bt709[1][1] *dat64f[i+1] + xyz2rgb_bt709[1][2] *dat64f[i+3];
        b = xyz2rgb_bt709[2][0] * dat64f[i] + xyz2rgb_bt709[2][1] *dat64f[i+1] + xyz2rgb_bt709[2][2] *dat64f[i+3];
        dat64f[i]  = r;
        dat64f[i+1]= g;
        dat64f[i+2]= b;
    }
    colorformat = 0;
    return 0;
}

int Jreadhdr::conv_rgb2Yxy(){
    if( colorformat != 0){
        fprintf(stderr," The color format of the image is not RGB\n");
        return -1;
    }

    double X, Y, Z, W, x, y;
    double r,g,b;
    for(size_t i=0; i < xsize*ysize*3; i+=3){
        r = dat64f[i  ];
        g = dat64f[i+1];
        b = dat64f[i+2];
        X = rgb2xyz_bt709[0][0] * dat64f[i] + rgb2xyz_bt709[0][1] *dat64f[i+1] + rgb2xyz_bt709[0][2] *dat64f[i+2];
        Y = rgb2xyz_bt709[1][0] * dat64f[i] + rgb2xyz_bt709[1][1] *dat64f[i+1] + rgb2xyz_bt709[1][2] *dat64f[i+2];
        Z = rgb2xyz_bt709[2][0] * dat64f[i] + rgb2xyz_bt709[2][1] *dat64f[i+1] + rgb2xyz_bt709[2][2] *dat64f[i+2];
        W = X + Y + Z;
        if( W <= 0.0) {
            x = 0.0;
            y = 0.0;
        }else{
            x = X/W;
            y = Y/W;
        }

        dat64f[i  ]= Y;
        dat64f[i+1]= x;
        dat64f[i+2]= y;
    }
    colorformat = 2;

    return 0;
}

int Jreadhdr::conv_Yxy2rgb(){
    if( colorformat != 2){
        fprintf(stderr," The color format of the image is not XYZ\n");
        return -1;
    }
    double r, g, b;
    double X,Y,Z,x,y,W;
    for(size_t i=0; i < xsize*ysize*3; i+=3){
        Y = dat64f[i  ];
        x = dat64f[i+1];
        y = dat64f[i+2];
        W = Y/y;
        X = x * W;
        Z = W-Y-X;
        r = xyz2rgb_bt709[0][0] *X + xyz2rgb_bt709[0][1] *Y + xyz2rgb_bt709[0][2] *Z;
        g = xyz2rgb_bt709[1][0] *X + xyz2rgb_bt709[1][1] *Y + xyz2rgb_bt709[1][2] *Z;
        b = xyz2rgb_bt709[2][0] *X + xyz2rgb_bt709[2][1] *Y + xyz2rgb_bt709[2][2] *Z;
        dat64f[i  ]= r;
        dat64f[i+1]= g;
        dat64f[i+2]= b;
    }
    colorformat = 0;
    return 0;
}
