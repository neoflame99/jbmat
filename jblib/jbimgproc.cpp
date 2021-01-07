/*
 * Copyright (C) 2020. Jong B. Choi
 * License : MIT License
 * contact : neoflame99@naver.com
 */

#include "jbimgproc.h"
#include "jbimgproc.h"


namespace jmat {

namespace imgproc {

Mat rgb2ycc(const Mat& rgbIm, const int32 sel_eq){
/*
 *   if sel_eq = 0 (BT 601)
 *  Y = [  0.299    ,  0.587    ,  0.114    ]   [ r ]
 * Pb = [ -0.168736 , -0.331264 ,  0.5      ] * [ g ]
 * Pr = [  0.5      , -0.418688 , -0.081312 ]   [ b ]
 *
 *   if sel_eq = 1 (BT 709)
 *  Y = [  0.2126  ,  0.7152  ,  0.0722  ]   [ r ]
 * Pb = [ -0.11457 , -0.38543 ,	 0.5     ] * [ g ]
 * Pr = [  0.5     , -0.45415 ,	-0.04585 ]   [ b ]
 *
 */
    uint32 row = rgbIm.getRow();
    uint32 col = rgbIm.getCol();
    uint32 chsize = rgbIm.getChannel();
    assert(chsize==3);

    Mat  A(rgbIm.getDatType(), row, col, chsize);
    int32 coe1[]={    2990,    5870,    1140,
                   -168736, -331264,  500000,
                    500000, -418688, -81312 };
    int32 coe2[]={    2126,    7152,    722,
                   -114570, -385430, 500000,
                    500000, -454150, -45850};
    int32* coe = sel_eq == 1 ? coe2 : coe1;
    uint32 x;
    if(rgbIm.getDatType()==DTYP::DOUBLE){
        bgr_d* bgr = (bgr_d*)rgbIm.getElptr().f64_ptr;
        yuv_d* yuv = (yuv_d*)A.getElptr().f64_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            yuv[x].y = (coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/10000; // Y
            yuv[x].u = (coe[3]*bgr[x].r + coe[4]*bgr[x].g + coe[5]*bgr[x].b)/1000000; // U
            yuv[x].v = (coe[6]*bgr[x].r + coe[7]*bgr[x].g + coe[8]*bgr[x].b)/1000000; // V
        }
    }else if(rgbIm.getDatType()==DTYP::FLOAT){
        bgr_f* bgr = (bgr_f*)rgbIm.getElptr().f32_ptr;
        yuv_f* yuv = (yuv_f*)A.getElptr().f32_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            yuv[x].y = (coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/10000  ; // Y
            yuv[x].u = (coe[3]*bgr[x].r + coe[4]*bgr[x].g + coe[5]*bgr[x].b)/1000000; // U
            yuv[x].v = (coe[6]*bgr[x].r + coe[7]*bgr[x].g + coe[8]*bgr[x].b)/1000000; // V
        }
    }else if(rgbIm.getDatType()==DTYP::INT){
        bgr_i* bgr = (bgr_i*)rgbIm.getElptr().int_ptr;
        yuv_i* yuv = (yuv_i*)A.getElptr().int_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            yuv[x].y = sat_cast<int32>((coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/10000  ); // Y
            yuv[x].u = sat_cast<int32>((coe[3]*bgr[x].r + coe[4]*bgr[x].g + coe[5]*bgr[x].b)/1000000); // U
            yuv[x].v = sat_cast<int32>((coe[6]*bgr[x].r + coe[7]*bgr[x].g + coe[8]*bgr[x].b)/1000000); // V
        }
    }else if(rgbIm.getDatType()==DTYP::UCHAR){
        bgr_uc* bgr = (bgr_uc*)rgbIm.getElptr().uch_ptr;
        yuv_uc* yuv = (yuv_uc*)A.getElptr().uch_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            yuv[x].y = sat_cast<uchar>((coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/10000  ); // Y
            yuv[x].u = sat_cast<uchar>((coe[3]*bgr[x].r + coe[4]*bgr[x].g + coe[5]*bgr[x].b)/1000000); // U
            yuv[x].v = sat_cast<uchar>((coe[6]*bgr[x].r + coe[7]*bgr[x].g + coe[8]*bgr[x].b)/1000000); // V
        }
    }else{
        fprintf(stderr, " Unsupported DTYP in rgb2ycc func.\n");
        return Mat();
    }
    return A;
}


Mat ycc2rgb(const Mat& yuvIm, const int32 sel_eq){
/*
 *   if sel_eq = 0 (BT 601)
 *
 *  r = [ 1.0000 , -0.0000  ,  1.4020 ]   [ Y  ]
 *  g = [ 1.0000 , -0.3441  , -0.7141 ] * [ Pb ]
 *  b = [ 1.0000 ,  1.7720  ,  0.0000 ]   [ Pr ]
 *
 *   if sel_eq = 1 (BT 709)
 *  r = [ 1.0000 ,  0.0000  ,   1.5748  ]   [ Y  ]
 *  g = [ 1.0000 , -0.1873  ,  -0.4681  ] * [ Pb ]
 *  b = [ 1.0000 ,  1.8556  ,  -0.0000  ]   [ Pr ]
 *
 */
    uint32 row = yuvIm.getRow();
    uint32 col = yuvIm.getCol();
    uint32 chsize = yuvIm.getChannel();
    assert(chsize==3);

    Mat  A(yuvIm.getDatType(), row, col, chsize);
    int32 coe1[]={ 10000,     0, 14020,
                   10000, -3441, -7141,
                   10000, 17720,     0};
    int32 coe2[]={ 10000,     0, 15748,
                   10000, -1873, -4681,
                   10000, 18556,     0};
    int32* coe = sel_eq == 1 ? coe2 : coe1;
    uint32 x;
    if(yuvIm.getDatType()==DTYP::DOUBLE){
        yuv_d* yuv = (yuv_d*)yuvIm.getElptr().f64_ptr;
        bgr_d* bgr = (bgr_d*)A.getElptr().f64_ptr;
        for(x=0 ; x < yuvIm.getRowColSize(); ++x ){
            bgr[x].r = (coe[0]*yuv[x].y + coe[1]*yuv[x].u + coe[2]*yuv[x].v)/10000;
            bgr[x].g = (coe[3]*yuv[x].y + coe[4]*yuv[x].u + coe[5]*yuv[x].v)/10000;
            bgr[x].b = (coe[6]*yuv[x].y + coe[7]*yuv[x].u + coe[8]*yuv[x].v)/10000;
            bgr[x].r = (bgr[x].r < 0 ) ? 0 : bgr[x].r;
            bgr[x].g = (bgr[x].g < 0 ) ? 0 : bgr[x].g;
            bgr[x].b = (bgr[x].b < 0 ) ? 0 : bgr[x].b;
        }
    }else if(yuvIm.getDatType()==DTYP::FLOAT){
        yuv_f* yuv = (yuv_f*)yuvIm.getElptr().f32_ptr;
        bgr_f* bgr = (bgr_f*)A.getElptr().f32_ptr;
        for(x=0 ; x < yuvIm.getRowColSize(); ++x ){
            bgr[x].r = (coe[0]*yuv[x].y + coe[1]*yuv[x].u + coe[2]*yuv[x].v)/10000;
            bgr[x].g = (coe[3]*yuv[x].y + coe[4]*yuv[x].u + coe[5]*yuv[x].v)/10000;
            bgr[x].b = (coe[6]*yuv[x].y + coe[7]*yuv[x].u + coe[8]*yuv[x].v)/10000;
            bgr[x].r = (bgr[x].r < 0 ) ? 0 : bgr[x].r;
            bgr[x].g = (bgr[x].g < 0 ) ? 0 : bgr[x].g;
            bgr[x].b = (bgr[x].b < 0 ) ? 0 : bgr[x].b;
        }
    }else if(yuvIm.getDatType()==DTYP::INT){
        yuv_i* yuv = (yuv_i*)yuvIm.getElptr().int_ptr;
        bgr_i* bgr = (bgr_i*)A.getElptr().int_ptr;
        for(x=0 ; x < yuvIm.getRowColSize(); ++x ){
            bgr[x].r = (coe[0]*yuv[x].y + coe[1]*yuv[x].u + coe[2]*yuv[x].v)/10000;
            bgr[x].g = (coe[3]*yuv[x].y + coe[4]*yuv[x].u + coe[5]*yuv[x].v)/10000;
            bgr[x].b = (coe[6]*yuv[x].y + coe[7]*yuv[x].u + coe[8]*yuv[x].v)/10000;
            bgr[x].r = (bgr[x].r < 0 ) ? 0 : bgr[x].r;
            bgr[x].g = (bgr[x].g < 0 ) ? 0 : bgr[x].g;
            bgr[x].b = (bgr[x].b < 0 ) ? 0 : bgr[x].b;
        }
    }else if(yuvIm.getDatType()==DTYP::UCHAR){
        yuv_uc* yuv = (yuv_uc*)yuvIm.getElptr().uch_ptr;
        bgr_uc* bgr = (bgr_uc*)A.getElptr().uch_ptr;
        for(x=0 ; x < yuvIm.getRowColSize(); ++x ){
            bgr[x].r = sat_cast<uchar>((coe[0]*yuv[x].y + coe[1]*yuv[x].u + coe[2]*yuv[x].v)/10000);
            bgr[x].g = sat_cast<uchar>((coe[3]*yuv[x].y + coe[4]*yuv[x].u + coe[5]*yuv[x].v)/10000);
            bgr[x].b = sat_cast<uchar>((coe[6]*yuv[x].y + coe[7]*yuv[x].u + coe[8]*yuv[x].v)/10000);
        }
    }else{
        fprintf(stderr, " Unsupported DTYP in ycc2rgb func.\n");
        return Mat();
    }
    return A;
}


Mat rgb2gray(const Mat& rgbIm, const int32 HowToGray){
/*
 *  if HowToGray = 0 (BT 601)
 *  Y =  0.299 * r   +  0.587 * g +  0.114 * b
 *
 *  if HowToGray = 1 (BT 709)
 *  Y =  0.2126 * r +  0.7152 * g +  0.0722 * b
 *
 *  if HowToGray = 2 (3 equal-weight)
 *  Y = 0.3333 * r + 0.3334 * g + 0.3333 * b
 */
    uint32 row = rgbIm.getRow();
    uint32 col = rgbIm.getCol();
    uint32 chsize = rgbIm.getChannel();
    assert(chsize==3);

    Mat  A(rgbIm.getDatType(), row, col, 1);
    int32 coe1[]={    2990,    5870,    1140};
    int32 coe2[]={    2126,    7152,     722};
    int32 coe3[]={    3333,    3334,    3333};
    int32* coe ;
    switch (HowToGray) {
    case 0 : coe = coe1; break;
    case 1 : coe = coe2; break;
    default: coe = coe3;
    }
    uint32 x;
    if(rgbIm.getDatType()==DTYP::DOUBLE){
        bgr_d* bgr = (bgr_d*)rgbIm.getElptr().f64_ptr;
        double* Y  = A.getElptr().f64_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            Y[x] = (coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/10000; // Y
        }
    }else if(rgbIm.getDatType()==DTYP::FLOAT){
        bgr_f* bgr = (bgr_f*)rgbIm.getElptr().f32_ptr;
        float* Y   = A.getElptr().f32_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            Y[x] = (coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/10000  ; // Y
        }
    }else if(rgbIm.getDatType()==DTYP::INT){
        bgr_i* bgr = (bgr_i*)rgbIm.getElptr().int_ptr;
        int32*   Y = A.getElptr().int_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            Y[x] = (coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/10000 ; // Y
        }
    }else if(rgbIm.getDatType()==DTYP::UCHAR){
        bgr_uc* bgr = (bgr_uc*)rgbIm.getElptr().uch_ptr;
        uchar*  Y   = A.getElptr().uch_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            Y[x] = sat_cast<uchar>((coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/10000  ); // Y
        }
    }else{
        fprintf(stderr, " Unsupported DTYP in rgb2gray func.\n");
        return Mat();
    }
    return A;
}

Mat rgb2xyz(const Mat& rgbIm){
   /*
    * rgb -> XYZ
    * [0.4124564, 0.3575761, 0.1804375] [ r ]
    * [0.2126729, 0.7151522, 0.0721750]*[ g ]
    * [0.0193339, 0.1191920, 0.9503041] [ b ]
    */
    uint32 row = rgbIm.getRow();
    uint32 col = rgbIm.getCol();
    uint32 chsize = rgbIm.getChannel();
    assert(chsize==3);

    Mat  A(rgbIm.getDatType(), row, col, chsize);
    int32 coe[]={412456, 357576, 180438,
                 212673, 715152,  72175,
                  19334, 119192, 950304};
    uint32 x;
    if(rgbIm.getDatType()==DTYP::DOUBLE){
        bgr_d* bgr = (bgr_d*)rgbIm.getElptr().f64_ptr;
        xyz_d* xyz = (xyz_d*)A.getElptr().f64_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            xyz[x].x = (coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/1000000; // X
            xyz[x].y = (coe[3]*bgr[x].r + coe[4]*bgr[x].g + coe[5]*bgr[x].b)/1000000; // Y
            xyz[x].z = (coe[6]*bgr[x].r + coe[7]*bgr[x].g + coe[8]*bgr[x].b)/1000000; // Z
        }
    }else if(rgbIm.getDatType()==DTYP::FLOAT){
        bgr_f* bgr = (bgr_f*)rgbIm.getElptr().f32_ptr;
        xyz_f* xyz = (xyz_f*)A.getElptr().f32_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            xyz[x].x = (coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/1000000; // X
            xyz[x].y = (coe[3]*bgr[x].r + coe[4]*bgr[x].g + coe[5]*bgr[x].b)/1000000; // Y
            xyz[x].z = (coe[6]*bgr[x].r + coe[7]*bgr[x].g + coe[8]*bgr[x].b)/1000000; // Z
        }
    }else if(rgbIm.getDatType()==DTYP::INT){
        bgr_i* bgr = (bgr_i*)rgbIm.getElptr().int_ptr;
        xyz_i* xyz = (xyz_i*)A.getElptr().int_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            xyz[x].x = (coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/1000000; // X
            xyz[x].y = (coe[3]*bgr[x].r + coe[4]*bgr[x].g + coe[5]*bgr[x].b)/1000000; // Y
            xyz[x].z = (coe[6]*bgr[x].r + coe[7]*bgr[x].g + coe[8]*bgr[x].b)/1000000; // Z
        }
    }else if(rgbIm.getDatType()==DTYP::UCHAR){
        bgr_uc* bgr = (bgr_uc*)rgbIm.getElptr().uch_ptr;
        xyz_uc* xyz = (xyz_uc*)A.getElptr().uch_ptr;
        for(x=0 ; x < rgbIm.getRowColSize(); ++x ){
            xyz[x].x = sat_cast<uchar>((coe[0]*bgr[x].r + coe[1]*bgr[x].g + coe[2]*bgr[x].b)/1000000); // X
            xyz[x].y = sat_cast<uchar>((coe[3]*bgr[x].r + coe[4]*bgr[x].g + coe[5]*bgr[x].b)/1000000); // Y
            xyz[x].z = sat_cast<uchar>((coe[6]*bgr[x].r + coe[7]*bgr[x].g + coe[8]*bgr[x].b)/1000000); // Z
        }
    }else{
        fprintf(stderr, " Unsupported DTYP in rgb2xyz func.\n");
        return Mat();
    }
    return A;
}

Mat xyz2rgb(const Mat& xyzIm){
   /*
    * XYZ -> rgb
    * [ 3.2404542, -1.5371385, -0.4985314]  [ X ]
    * [-0.9692660,  1.8760108,  0.0415560] *[ Y ]
    * [ 0.0556434, -0.2040259,  1.0572252]  [ Z ]
    */
    uint32 row = xyzIm.getRow();
    uint32 col = xyzIm.getCol();
    uint32 chsize = xyzIm.getChannel();
    assert(chsize==3);

    Mat  A(xyzIm.getDatType(), row, col, chsize);
    int32 coe[]={ 3240454, -1537139, -498531,
                  -969266,  1876011,   41556,
                    55643,  -204026, 1057225};
    uint32 x;
    if(xyzIm.getDatType()==DTYP::DOUBLE){
        xyz_d* xyz = (xyz_d*)xyzIm.getElptr().f64_ptr;
        bgr_d* bgr = (bgr_d*)A.getElptr().f64_ptr;
        for(x=0 ; x < xyzIm.getRowColSize(); ++x ){
            bgr[x].r = (coe[0]*xyz[x].x + coe[1]*xyz[x].y + coe[2]*xyz[x].z)/1000000; // R
            bgr[x].g = (coe[3]*xyz[x].x + coe[4]*xyz[x].y + coe[5]*xyz[x].z)/1000000; // G
            bgr[x].b = (coe[6]*xyz[x].x + coe[7]*xyz[x].y + coe[8]*xyz[x].z)/1000000; // B
            bgr[x].r = bgr[x].r < 0 ? 0 : bgr[x].r;
            bgr[x].g = bgr[x].g < 0 ? 0 : bgr[x].g;
            bgr[x].b = bgr[x].b < 0 ? 0 : bgr[x].b;
        }
    }else if(xyzIm.getDatType()==DTYP::FLOAT){
        xyz_f* xyz = (xyz_f*)xyzIm.getElptr().f32_ptr;
        bgr_f* bgr = (bgr_f*)A.getElptr().f32_ptr;
        for(x=0 ; x < xyzIm.getRowColSize(); ++x ){
            bgr[x].r = (coe[0]*xyz[x].x + coe[1]*xyz[x].y + coe[2]*xyz[x].z)/1000000; // R
            bgr[x].g = (coe[3]*xyz[x].x + coe[4]*xyz[x].y + coe[5]*xyz[x].z)/1000000; // G
            bgr[x].b = (coe[6]*xyz[x].x + coe[7]*xyz[x].y + coe[8]*xyz[x].z)/1000000; // B
            bgr[x].r = bgr[x].r < 0 ? 0 : bgr[x].r;
            bgr[x].g = bgr[x].g < 0 ? 0 : bgr[x].g;
            bgr[x].b = bgr[x].b < 0 ? 0 : bgr[x].b;
        }
    }else if(xyzIm.getDatType()==DTYP::INT){
        xyz_i* xyz = (xyz_i*)xyzIm.getElptr().int_ptr;
        bgr_i* bgr = (bgr_i*)A.getElptr().int_ptr;
        for(x=0 ; x < xyzIm.getRowColSize(); ++x ){
            bgr[x].r = (coe[0]*xyz[x].x + coe[1]*xyz[x].y + coe[2]*xyz[x].z)/1000000; // R
            bgr[x].g = (coe[3]*xyz[x].x + coe[4]*xyz[x].y + coe[5]*xyz[x].z)/1000000; // G
            bgr[x].b = (coe[6]*xyz[x].x + coe[7]*xyz[x].y + coe[8]*xyz[x].z)/1000000; // B
            bgr[x].r = bgr[x].r < 0 ? 0 : bgr[x].r;
            bgr[x].g = bgr[x].g < 0 ? 0 : bgr[x].g;
            bgr[x].b = bgr[x].b < 0 ? 0 : bgr[x].b;
        }
    }else if(xyzIm.getDatType()==DTYP::UCHAR){
        xyz_uc* xyz = (xyz_uc*)xyzIm.getElptr().uch_ptr;
        bgr_uc* bgr = (bgr_uc*)A.getElptr().uch_ptr;
        for(x=0 ; x < xyzIm.getRowColSize(); ++x ){
            bgr[x].r = sat_cast<uchar>((coe[0]*xyz[x].x + coe[1]*xyz[x].y + coe[2]*xyz[x].z)/1000000); // R
            bgr[x].g = sat_cast<uchar>((coe[3]*xyz[x].x + coe[4]*xyz[x].y + coe[5]*xyz[x].z)/1000000); // G
            bgr[x].b = sat_cast<uchar>((coe[6]*xyz[x].x + coe[7]*xyz[x].y + coe[8]*xyz[x].z)/1000000); // B
        }
    }else{
        fprintf(stderr, " Unsupported DTYP in xyz2rgb func.\n");
        return Mat();
    }
    return A;
}

Mat rgb2Yxy(const Mat& rgbIm){
    uint32 chsize = rgbIm.getChannel();
    if( chsize != 3) return Mat();

    Mat A = rgb2xyz(rgbIm);

    double W, x, y;
    uint32 i;
    if(A.getDatType()==DTYP::DOUBLE){
        xyz_d* xyz = (xyz_d*)A.getElptr().f64_ptr;
        for( i=0 ; i < A.getRowColSize(); ++i ){
            W = xyz[i].x + xyz[i].y + xyz[i].z;
            x = xyz[i].x/W;
            y = xyz[i].y/W;
            xyz[i].x = xyz[i].y;
            xyz[i].y = x;
            xyz[i].z = y;
        }
    }else if(A.getDatType()==DTYP::FLOAT){
        xyz_f* xyz = (xyz_f*)A.getElptr().f32_ptr;
        for( i=0 ; i < A.getRowColSize(); ++i ){
            W = xyz[i].x + xyz[i].y + xyz[i].z;
            x = xyz[i].x/W;
            y = xyz[i].y/W;
            xyz[i].x = xyz[i].y;
            xyz[i].y = x;
            xyz[i].z = y;
        }
    }else if(A.getDatType()==DTYP::INT){
        xyz_i* xyz = (xyz_i*)A.getElptr().int_ptr;
        for( i=0 ; i < A.getRowColSize(); ++i ){
            W = xyz[i].x + xyz[i].y + xyz[i].z;
            x = (double)xyz[i].x/W;
            y = (double)xyz[i].y/W;
            xyz[i].x = xyz[i].y;
            xyz[i].y = x;
            xyz[i].z = y;
        }
    }else if(A.getDatType()==DTYP::FLOAT){
        xyz_f* xyz = (xyz_f*)A.getElptr().f32_ptr;
        for( i=0 ; i < A.getRowColSize(); ++i ){
            W = xyz[i].x + xyz[i].y + xyz[i].z;
            x = (double)xyz[i].x/W;
            y = (double)xyz[i].y/W;
            xyz[i].x = xyz[i].y;
            xyz[i].y = x;
            xyz[i].z = y;
        }
    }else{
        fprintf(stderr, " Unsupported DTYP in rgb2Yxy func.\n");
        return Mat();
    }
    return A;
}

Mat Yxy2rgb(const Mat& YxyIm){
    uint32 row    = YxyIm.getRow();
    uint32 col    = YxyIm.getCol();
    uint32 chsize = YxyIm.getChannel();
    if( chsize != 3) return Mat();

    Mat A(YxyIm.getDatType(), row, col, chsize);

    double W, X, Z;
    uint32 i;
    if(YxyIm.getDatType()==DTYP::DOUBLE){
        Yxy_d* Yxy = (Yxy_d*)YxyIm.getElptr().f64_ptr;
        xyz_d* xyz = (xyz_d*)A.getElptr().f64_ptr;
        for(i=0 ; i < YxyIm.getRowColSize(); ++i){
            W = Yxy[i].Y/Yxy[i].y;
            X = Yxy[i].x * W;
            Z = W-Yxy[i].Y-X;
            xyz[i].y = Yxy[i].Y;
            xyz[i].x = X;
            xyz[i].z = Z;
        }
    }else if(YxyIm.getDatType()==DTYP::FLOAT){
        Yxy_f* Yxy = (Yxy_f*)YxyIm.getElptr().f32_ptr;
        xyz_f* xyz = (xyz_f*)A.getElptr().f32_ptr;
        for(i=0 ; i < YxyIm.getRowColSize(); ++i){
            W = Yxy[i].Y/Yxy[i].y;
            X = Yxy[i].x * W;
            Z = W-Yxy[i].Y-X;
            xyz[i].y = Yxy[i].Y;
            xyz[i].x = X;
            xyz[i].z = Z;
        }
    }else if(YxyIm.getDatType()==DTYP::INT){
        Yxy_i* Yxy = (Yxy_i*)YxyIm.getElptr().int_ptr;
        xyz_i* xyz = (xyz_i*)A.getElptr().int_ptr;
        for(i=0 ; i < YxyIm.getRowColSize(); ++i){
            W = (double)Yxy[i].Y/Yxy[i].y;
            X = (double)Yxy[i].x * W;
            Z = (double)W-Yxy[i].Y-X;
            xyz[i].y = Yxy[i].Y;
            xyz[i].x = X;
            xyz[i].z = Z;
        }
    }else if(YxyIm.getDatType()==DTYP::UCHAR){
        Yxy_uc* Yxy = (Yxy_uc*)YxyIm.getElptr().uch_ptr;
        xyz_uc* xyz = (xyz_uc*)A.getElptr().uch_ptr;
        for(i=0 ; i < YxyIm.getRowColSize(); ++i){
            W = (double)Yxy[i].Y/Yxy[i].y;
            X = (double)Yxy[i].x * W;
            Z = (double)W-Yxy[i].Y-X;
            xyz[i].y = Yxy[i].Y;
            xyz[i].x = X;
            xyz[i].z = Z;
        }
    }else{
        fprintf(stderr, " Unsupported DTYP in Yxy2rgb func.\n");
        return Mat();
    }

    return xyz2rgb(A);
}


Mat histoPmf(const Mat& src, const uint32 bins, const double step, const double low_clipval){
    uint32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stderr,"histoPmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stderr,"histoPmf : src is not 1 channel matrix\n");
        return Mat();
    }else if( bins < 1) {
        fprintf(stderr,"histoPmf : 'bins' should be larger than or equal 1 \n");
        return Mat();
    }

    Mat A = Mat::zeros(1, bins, src.getChannel(), DTYP::DOUBLE);
    elemptr tarDat_pt = A.getElptr();
    elemptr srcDat_pt = src.getElptr();

    uint32 k, m, d ;

    if(src.getDatType()==DTYP::DOUBLE){
        for(m=0 ; m < ch ; ++m){
            for(k=m ; k < src.getLength() ; k+=ch ){
                d = uint32((srcDat_pt.f64_ptr[k]-low_clipval)/step);
                d = (d < 0) ? 0 : (d >= bins) ? bins-1: d;  // bin index range : 0 ~ bins-1
                tarDat_pt.f64_ptr[d*ch+m]++;
            }
        }
    }else if(src.getDatType()==DTYP::FLOAT){
        for(m=0 ; m < ch ; ++m){
            for(k=m ; k < src.getLength() ; k+=ch ){
                d = uint32((srcDat_pt.f32_ptr[k]-low_clipval)/step);
                d = (d < 0) ? 0 : (d >= bins) ? bins-1: d;  // bin index range : 0 ~ bins-1
                tarDat_pt.f64_ptr[d*ch+m]++;
            }
        }
    }else if(src.getDatType()==DTYP::INT){
        for(m=0 ; m < ch ; ++m){
            for(k=m ; k < src.getLength() ; k+=ch ){
                d = uint32((srcDat_pt.int_ptr[k]-low_clipval)/step);
                d = (d < 0) ? 0 : (d >= bins) ? bins-1: d;  // bin index range : 0 ~ bins-1
                tarDat_pt.f64_ptr[d*ch+m]++;
            }
        }
    }else if(src.getDatType()==DTYP::UCHAR){
        for(m=0 ; m < ch ; ++m){
            for(k=m ; k < src.getLength() ; k+=ch ){
                d = uint32((srcDat_pt.uch_ptr[k]-low_clipval)/step);
                d = (d < 0) ? 0 : (d >= bins) ? bins-1: d;  // bin index range : 0 ~ bins-1
                tarDat_pt.f64_ptr[d*ch+m]++;
            }
        }
    }else{
        fprintf(stderr, " Unsupported DTYP in histoPmf func.\n");
        return Mat();
    }
    return A;
}

Mat histoCmf(const Mat& src, const uint32 bins, const double step, const double low_clipval){
    uint32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stderr,"histoCmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stderr,"histoCmf : src is not 1 channel matrix\n");
        return Mat();
    }else if( bins < 1) {
        fprintf(stderr,"histoCmf : 'bins' should be larger than or equal 1 \n");
        return Mat();
    }else if( src.getDatType()==DTYP::CMPLX){
        fprintf(stderr, " Unsupported DTYP in histoCmf func.\n");
        return Mat();
    }

    Mat cmf = histoPmf(src, bins, step, low_clipval);

    double *srcDat_pt = cmf.getDataPtr<double>();
    uint32 ch_cmf = cmf.getChannel();
    uint32 m, k;
    for( m=0; m < ch_cmf; ++m){
        for( k=ch_cmf+m ; k < cmf.getLength() ; k+= ch_cmf)  // getLength() == bins*ch
            srcDat_pt[k] += srcDat_pt[k-ch];
    }

    return cmf;
}


Mat clip_HistoPmf(const Mat& src,const uint32 clipVal,const uint32 bins, const double step, const double low_clipval){
    uint32 ch  = src.getChannel();
    if( src.isEmpty() ){
        fprintf(stderr,"clip_histoPmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stderr,"clip_histoPmf : src is not 1 channel matrix\n");
        return Mat();
    }else if( bins < 1) {
        fprintf(stderr,"clip_histoPmf : 'bins' should be larger than or equal 1 \n");
        return Mat();
    }else if( step < 1) {
        fprintf(stderr,"clip_histoPmf : 'step' should be larger than or equal 1 \n");
        return Mat();
    }else if( src.getDatType()==DTYP::CMPLX){
        fprintf(stderr, " Unsupported DTYP in clip_HistoPmf func.\n");
        return Mat();
    }

    Mat pmf = histoPmf(src, bins, step, low_clipval );
    double *srcDat_pt = pmf.getDataPtr<double>();

    uint32 ch_pmf = pmf.getChannel();
    uint32 sum_clipped ;
    uint32 binval;
    uint32 k, m;

    for( m =0 ; m < ch_pmf ; ++m){
        // clipping on each channel
        sum_clipped = 0;
        for( k=m ; k < pmf.getLength() ; k+= ch_pmf){ // getLength() == bins*ch
            binval = int32(srcDat_pt[k]);
            if( binval > clipVal){
                sum_clipped += binval - clipVal;
                srcDat_pt[k] = clipVal;
            }
        }
        sum_clipped /= bins;
        // distributing the clipped sum
        for( k=m ; k < pmf.getLength() ; k+= ch_pmf){
            srcDat_pt[k] += sum_clipped;
        }
    }
    return pmf;
}

Mat clip_HistoCmf(const Mat& src,const uint32 clipVal,const uint32 bins, const double step, const double low_clipval){
    uint32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stderr,"clip_histoCmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stderr,"clip_histoCmf : src is not 1 channel matrix\n");
        return Mat();
    }else if( bins < 1) {
        fprintf(stderr,"clip_histoPmf : 'bins' should be larger than or equal 1 \n");
        return Mat();
    }else if( step < 1) {
        fprintf(stderr,"clip_histoPmf : 'step' should be larger than or equal 1 \n");
        return Mat();
    }else if( src.getDatType()==DTYP::CMPLX){
        fprintf(stderr, " Unsupported DTYP in clip_HistoCmf func.\n");
        return Mat();
    }

    Mat cmf = clip_HistoPmf(src, clipVal, bins, step, low_clipval);
    double *srcDat_pt = cmf.getDataPtr<double>();
    uint32 ch_cmf = cmf.getChannel();
    uint32 m, k;

    // making cumiltive data
    for( m = 0 ; m < ch_cmf ; ++m){
        for(k=ch_cmf+m ; k < cmf.getLength() ; k+= ch) // getLength() == bins*ch
            srcDat_pt[k] += srcDat_pt[k-ch];
    }
    return cmf;
}

Mat clip_HistoEq(const Mat& src,const Mat& histCmf, const double step){
    uint32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stdout,"histoEqual : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stdout,"histoEqual : src is not 1 channel matrix\n");
        return Mat();
    }else if( histCmf.isEmpty()){
        fprintf(stdout,"histCmf is empty \n");
        return Mat();
    }
    switch(src.getDatType()){
    case DTYP::DOUBLE : return _clip_HistoEq<double>(src, histCmf, step);
    case DTYP::FLOAT  : return _clip_HistoEq<float >(src, histCmf, step);
    case DTYP::INT    : return _clip_HistoEq<int32 >(src, histCmf, step);
    case DTYP::UCHAR  : return _clip_HistoEq<uchar >(src, histCmf, step);
    default           : {
        fprintf(stderr, " Unsupported DTYP in clip_HistoEqual func.\n");
        return Mat();
        }
    }
}


Mat gaussMaskGen (const double sigma, const double factor, const uint32 ch ){
    // mask size sigma*factor*2+1
    uint32 hp = static_cast<uint32>( sigma*factor);
    uint32 sz = static_cast<uint32>( (hp<<1) +1 );

    Mat mask = Mat::zeros(sz, sz, ch, DTYP::DOUBLE);
    uint32 y, x, cc;
    int32 yp, xp;
    for( cc =0 ; cc < ch ; ++cc ){
        for (y=0, yp=-static_cast<int32>(hp) ; y < sz ; ++y, ++yp){
            for(x=0, xp=-static_cast<int32>(hp) ; x < sz; ++x, ++xp)
                mask.at<double>(y, x, cc) = exp(-((xp*xp + yp*yp)/(2*sigma*sigma)));
        }
    }
    Mat sum = mask.sum();
    for( cc =0 ; cc < ch ; ++cc ){
        for (y=0, yp=-static_cast<int32>(hp) ; y < sz ; ++y, ++yp){
            for(x=0, xp=-static_cast<int32>(hp) ; x < sz; ++x, ++xp)
                mask.at<double>(y, x, cc) /= sum.at<double>(cc);
        }
    }

    return mask;
}

Mat boxMaskGen( const uint32 sz, const uint32 ch){
    Mat mask = Mat::ones(sz,sz,ch,DTYP::DOUBLE);
    mask /= (sz*sz);
    return mask;
}

Mat nakaSigTonemap( Mat& src, Mat& localmean, const double globalmean, const double Imax){
    DTYP srcDtype  = src.getDatType();
    DTYP lmeanDtype = localmean.getDatType();

    if(srcDtype != lmeanDtype){
        fprintf(stderr,"the data types between src and localmean aren't the same!\n ");
        return Mat();
    }

    Mat A;
    switch( srcDtype ){
    case DTYP::DOUBLE :
        A   = _nakaSigTm<double>( src, localmean, globalmean, Imax );
        break;
    case DTYP::FLOAT  :
        A   = _nakaSigTm<float >( src, localmean, globalmean, Imax);
        break;
    case DTYP::INT    :
        A   = _nakaSigTm<int32 >( src, localmean, globalmean, Imax);
        break;
    default:
        fprintf(stderr,"Unsupproted data type in nakaSigtonemap\n ");
    }

    return A;
}

Mat nakaSig3MeanTonemap( Mat& src, Mat& s_localmean, Mat& l_localmean, const double globalmean, const double Imax){
    DTYP srcDtype  = src.getDatType();
    DTYP smeanDtype = s_localmean.getDatType();
    DTYP lmeanDtype = l_localmean.getDatType();

    if(srcDtype != smeanDtype){
        fprintf(stderr,"the data types between src and s_localmean aren't the same!\n ");
        return Mat();
    }else if( srcDtype != lmeanDtype){
        fprintf(stderr,"the data types between src and l_localmean aren't the same!\n ");
        return Mat();
    }

    Mat A;
    switch( srcDtype ){
    case DTYP::DOUBLE :
        A   = _nakaSig3MeanTm<double>( src, s_localmean, l_localmean, globalmean, Imax );
        break;
    case DTYP::FLOAT  :
        A   = _nakaSig3MeanTm<float >( src, s_localmean, l_localmean, globalmean, Imax);
        break;
    case DTYP::INT    :
        A   = _nakaSig3MeanTm<int32 >( src, s_localmean, l_localmean, globalmean, Imax);
        break;
    default:
        fprintf(stderr,"Unsupproted data type in nakaSig3MeanTonemap\n ");
    }

    return A;
}

Mat logRetinexTonemap( Mat& src, Mat& surround){
    DTYP srcDtype  = src.getDatType();
    DTYP surrDtype = surround.getDatType();

    if(srcDtype != surrDtype){
        fprintf(stderr,"the data types between src and surround aren't the same!\n ");
        return Mat();
    }

    Mat A;
    switch( srcDtype ){
    case DTYP::DOUBLE :
        A   = _logRetinexTm<double>( src, surround );
        break;
    case DTYP::FLOAT  :
        A   = _logRetinexTm<float >( src, surround );
        break;
    case DTYP::INT    :
        A   = _logRetinexTm<int32 >( src, surround );
        break;
    default:
        fprintf(stderr,"Unsupproted data type in nakaSigtonemap\n ");
    }

    return A;
}

Mat gamma( const Mat& src, const double gmvalue){
    DTYP srcDtype = src.getDatType();

    switch( srcDtype){
    case DTYP::DOUBLE : return _gamma<double>(src, gmvalue);
    case DTYP::FLOAT  : return _gamma<float >(src, gmvalue);
    case DTYP::INT    : return _gamma<int32 >(src, gmvalue);
    case DTYP::UCHAR  : return _gamma<uchar >(src, gmvalue);
    default: fprintf(stderr, "Unsuppretd data type in gamma func.\n"); return Mat();
    }
}
inline void bitrev_permute(_complex* dat, int32 len){
    //-- data shuffle by bit reverse  -----------
    int32 i, k;
    int32 hlen = len >> 1;
    int32 t = 0;
    for (i=0; i< len-1; i++) {
        if (i < t) std::swap(dat[i], dat[t]);
         k = hlen;
         while (k <= t) {
             t -= k;
             k >>= 1;
         }
         t += k;
    }
    //-------- shuffle finished -----------------
}

void fft_dit2( _complex *dat, int32 len, bool backward){
// len has to be a number of power of 2
    assert( dat != nullptr);
    assert( len > 0 );

    int32 k, i;
    int32 n=0;
    k = len;
    while(k > 1){
        k >>= 1; ++n;
    }

    int32 t, pn;
    int32 mh = len >> 2;
#ifdef FFT_EXP_TABLE
    _complex *e_arr;
 #ifdef MALLOC_F
    //e_arr = (_complex *)calloc(len+1,sizeof(_complex));
    e_arr = (_complex *)calloc(mh+1,sizeof(_complex));
    if( e_arr == nullptr){
        fprintf(stderr, " memories for sin and con table are not allocated ! \n");
        return ;
    }
 #else
    try{
        //e_arr = new _complex[len+1];
        e_arr = new _complex[mh+1];
    }catch(std::bad_alloc& ex){
        fprintf(stderr, " memories for sin and con table are not allocated : %s\n", ex.what());
        return ;
    }
 #endif // MALLOC_F
    double PI_2 = M_PI*2;
    double cosq ;
    for(k=1; k <= mh ; ++k){
        cosq = cos(PI_2 * k/len); // cos is even function
        e_arr[k    ].re =  cosq;

        cosq = (backward) ? cosq : -cosq; // make sin table with cos
        e_arr[mh-k ].im =  cosq;
    }
    e_arr[0  ] = _complex( 1, 0);
    e_arr[mh ] = (backward)? _complex(0, 1) : _complex( 0,-1);
    //for(k=0;k < mh; ++k)
    //    fprintf(stdout, "%d : %8.8f + %8.8f j\n", k, e_arr[k].re, e_arr[k].im);
#else
    double theta ;
    double dirPI = (backward) ? 2*M_PI : -2*M_PI;
#endif
    //-- data shuffle by bit reverse  -----------
    bitrev_permute(dat,len);
    //-------- shuffle finished -----------------

    _complex tmp;
   // extract first loop from following loop
    for(i=0; i < len; i+=2){
        tmp      = dat[i] + dat[i+1];
        dat[i+1] = dat[i] - dat[i+1];
        dat[i  ] = tmp;
    }
    //for(pn = 2; pn <= len; pn <<=1 ){ // 'pn' is partial len at the step.
    for(pn = 4; pn <= len; pn <<=1 ){ // 'pn' is partial len at the step.
#ifdef FFT_EXP_TABLE
        t = len/pn;
        _complex ws = e_arr[t];
#else
        theta = dirPI /pn ;     // (direc * 2) * M_PI /pn;
        _complex ws(cos(theta), sin(theta));
#endif
        for(i=0; i < len; i += pn){   // group-wise at n-th step loop
            _complex w(1,0);
            int32 half_pn = pn >> 1;
            //-- Butterfly
            for(k=0; k < half_pn ; ++k){
                tmp = dat[i+k+ half_pn] * w;
                dat[i+k+half_pn] = dat[i+k] - tmp;
                dat[i+k] += tmp;
                w *= ws;
            }
        }
    }

    if( backward ){
        for( i=0 ; i < len; ++i)
            dat[i] /= len;
    }
#ifdef FFT_EXP_TABLE
 #ifdef MALLOC_F
    free(e_arr);
 #else
    delete [] e_arr;
 #endif
#endif
}

void fft_czt( _complex *dat, int32 len, bool inverse){
    if( dat == nullptr){
        fprintf(stderr, " argument (_complex *dat) of fft_p2 is NULL\n ");
        return ;
    }else if(len <= 0){
        fprintf(stderr, " argument len of fft_p2 is less than or equal to zero\n ");
        return ;
    }

    int32 N  = len;
    int32 cN = (N << 1)-1; // 2N - 1
    int32 N2 = 0;          // where N2 > 2N-1 and N2 is the power of 2.
    int32 k, i;

    k = cN;
    while(k > 1){
        k >>= 1; ++N2;
    }
    N2 = (1 << N2) < cN ? N2+1 : N2;
    N2 = 1 << N2;

    _complex*  chirp;
    _complex* ichirp;
    _complex* extdat;
#ifdef MALLOC_F
    chirp  = (_complex *)calloc(cN, sizeof(_complex));
    ichirp = (_complex *)calloc(N2, sizeof(_complex));
    extdat = (_complex *)calloc(N2, sizeof(_complex));
    if( chirp == nullptr || ichirp == nullptr || extdat == nullptr){
        fprintf(stderr,"memory allocation error!\n");
        return;
    }
#else
    try{
        chirp= new _complex[static_cast<uint32>(cN)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"Chirp memory bad allocation!: %s\n",ex.what());
        return;
    }
    try{
        ichirp= new _complex[static_cast<uint32>(N2)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"IChirp memory bad allocation!: %s\n",ex.what());
        return;
    }
    try{
        extdat= new _complex[static_cast<uint32>(N2)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"extdat memory bad allocation!: %s\n",ex.what());
        return;
    }
#endif


    // filling chirp values ; its indices: -N+1, -N+2, ..., 0 , ..., N-1
    // forward  : W**(k**2)/(-2) = exp(j2*PI/N*(k**2)/(-2)) = exp(-j*PI/N*(k**2))
    // backward : W**(k**2)/2 = exp(j2*PI/N*(k**2)/2) = exp(j*PI/N*(k**2))
    double unit_theta = inverse  ? -M_PI/N : M_PI/N;
    double theta;
    for(k=1-N , i=0 ; k <N ; ++k, ++i ){
        theta = unit_theta*(k*k);
        chirp[i] = _complex(cos(theta),-sin(theta));
    }

    // filling ichirp( inverse of chirp )
    // we need zero padding but 'new' operator calls the default constructor so each element of array is zero initilized.
    for(i=0; i < cN; ++i)
        ichirp[i] = 1.0/chirp[i];

    // filling extdat; its indices: 0, 1, 2,..., N-1
    // forward  : x(n)*exp(2j*PI/N*(n**2)/(-2)) = x(n)*exp(-j*PI/N*(n**2))
    // backward : X(n)*exp(2j*PI/N*(n**2)/2) = X(n)*exp(j*PI/N*(n**2))
    // we need zero padding but 'new' operator calls the default constructor so each element of array is zero initilized.
    int32 N_minus_1 = N-1;
    for(i=0; i < N ; ++i)
        extdat[i] = dat[i] * chirp[i+N_minus_1];

    // fft through fft_dit2
    fft_dit2(extdat, N2, false);
    fft_dit2(ichirp, N2, false);

    // convolution in time or sample domain
    for(i=0; i < N2; ++i)
        extdat[i] *= ichirp[i];

    // inverse fft
    fft_dit2(extdat, N2, true);

    for(i=0, k=N-1; i < N; ++i, ++k)
        dat[i] = extdat[k] * chirp[k];

    if( inverse ){
        for( i=0 ; i < N; ++i)
            dat[i] /= N;
    }

#ifdef MALLOC_F
    free( chirp);
    free(ichirp);
    free(extdat);
#else
    delete [] chirp;
    delete [] ichirp;
    delete [] extdat;
#endif

}

void  fft(_complex* dat, int32 len){
    int32 k, i;

    k = len;
    i = 0;
    while(k > 1){
        k >>= 1; ++i;
    }
    i = 1 << i;
    if( i == len ) // N is a power of 2
        fft_dit2(dat, len, false);
    else{
        std::vector<int32> fac;
        factorizeN(len,fac);
        if( !fac.empty() && fac.at(fac.size()-1) <= 61 )
            fft_compositN(dat,len,fac,false);
        else
            fft_czt(dat, len, false);
    }

}
void  ifft(_complex* dat, int32 len){
    int32 k, i;

    k = len;
    i = 0;
    while(k > 1){
        k >>= 1; ++i;
    }
    i = 1 << i;
    if( i == len ) // N is a power of 2
        fft_dit2(dat, len, true);
    else{
        std::vector<int32> fac;
        factorizeN(len,fac);
        if( !fac.empty() && fac.at(fac.size()-1) <= 61 )
            fft_compositN(dat,len,fac, true);
        else
            fft_czt(dat, len, true);
    }
}

void fft2d(_complex* dat, int32 r_len, int32 c_len){

    if( dat == nullptr ){
        fprintf(stderr, " dat argument into fft2d is NULL!\n"); return ;
    }else if( r_len <=0 || c_len <=0 ){
        fprintf(stderr, " both or one of r_len or c_len arguments of fft2d is zero or negative!\n"); return ;
    }


    _complex* rdat;
#ifdef MALLOC_F
    rdat = (_complex*)calloc(r_len, sizeof(_complex));
    if( rdat == nullptr){
        fprintf(stderr, " data memory for fft2d is not allocated \n"); return ;
    }
#else
    try{
        rdat = new _complex[static_cast<uint32>(r_len)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr, " data memory for fft2d is not allocated : %s\n", ex.what()); return ;
    }
#endif

    int32 i, k, r;
    _complex* cdat;
    // fft on column direction
    for(i=0, cdat = dat; i < r_len; ++i, cdat += c_len){
        fft(cdat, c_len);
    }
    // fft on row direction
    for(i=0; i < c_len; ++i ){
        for(k=0, r=i; k < r_len; ++k, r+= c_len){
            rdat[k] = dat[r];
        }
        fft(rdat, r_len);
        for(k=0, r=i; k < r_len; ++k, r+= c_len){
            dat[r] = rdat[k];
        }
    }
#ifdef MALLOC_F
    free( rdat );
#else
    delete [] rdat;
#endif
}

void ifft2d(_complex* dat, int32 r_len, int32 c_len){

    if( dat == nullptr ){
        fprintf(stderr, " dat argument into ifft2d is NULL!\n"); return ;
    }else if( r_len <=0 || c_len <=0 ){
        fprintf(stderr, " both or one of r_len or c_len arguments of ifft2d is zero or negative!\n"); return ;
    }

    _complex* rdat;
#ifdef MALLOC_F
    rdat = (_complex*)calloc(r_len, sizeof(_complex));
    if( rdat == nullptr){
        fprintf(stderr, " data memory for fft2d is not allocated \n"); return ;
    }
#else
    try{
        rdat = new _complex[static_cast<uint32>(r_len)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr, " data memory for ifft2d is not allocated : %s\n", ex.what()); return ;
    }
#endif

    if( rdat == nullptr){
        fprintf(stderr, " data memory for ifft2d is not allocated\n"); return ;
    }

    int32 i, k, r;
    _complex* cdat;
    // ifft on column direction
    for(i=0, cdat = dat; i < r_len; ++i, cdat += c_len){
        ifft(cdat, c_len);
    }
    // fft on row direction
    for(i=0; i < c_len; ++i ){
        for(k=0, r=i; k < r_len; ++k, r+= c_len){
            rdat[k] = dat[r];
        }
        ifft(rdat, r_len);
        for(k=0, r=i; k < r_len; ++k, r+= c_len){
            dat[r] = rdat[k];
        }
    }
#ifdef MALLOC_F
    free( rdat);
#else
    delete [] rdat;
#endif

}
inline int32 digit4_rev(int32 x, int32 ldn, int32 radbit, int32 andBit){
    int32 j = 0;
    while ( ldn > 0){
        j <<= radbit;
        j += (x & andBit);
        x >>= radbit;
        ldn -= radbit;
    }
    return j;
}
void permute_radix4(_complex *a, int32 len){
    int32 r;

    int32 k = len;
    int32 ldn = 0;
    while(k > 1){
        k >>= 1; ++ldn; // ldn -> log2(len)
    }

    int32 radbit = 2; // log2(radix)
    int32 andbit = 3; // radix-1
    for(int32 x=0; x < len ; x++){
        r = digit4_rev(x, ldn, radbit, andbit);
        if( r > x )
            std::swap(a[x], a[r]);
    }
}
void fft_dif4(_complex *dat, int32 len, bool backward){
/* referece :
 *  1) FFTs for programmers: alorithms and source code. Jorg Arndt
 *  2) Radix-4 DIF FFT Algorithm, https://hackmd.io/@akshayk07/ryn-yR7qr
 */

// len has to be a number of power of 4
    assert( dat != nullptr);
    assert( len > 0 );

    int32 k;
    int32 n=0;
    k = len;
    while(k > 1){
        k >>= 1; ++n;
    }
    int32 ldn = n >> 1; // log4(n)

    int32 mh, mh2, mh3;
    int32 m = len;
    _complex e, e2, e3;
    _complex u0, u1, u2, u3;
    _complex x, y, t0, t1, t2, t3;
#ifdef FFT_EXP_TABLE
    _complex *e_arr;
 #ifdef MALLOC_F
    e_arr = (_complex *)calloc(len+1, sizeof(_complex));
    if( e_arr == nullptr ){
        fprintf(stderr, "sin, cos table memory allocation error\n!");
        return ;
    }
 #else
    try {
        e_arr = new _complex[len+1];
    } catch ( std::bad_alloc& ex) {
        fprintf(stderr, "sin, cos table memory allocation error!: %s\n", ex.what());
        return;
    }
 #endif // end of MALLOC_F
    double PI_2 = M_PI*2;
    double cosq ;
    mh  = len >> 2;
    mh2 = len >> 1;
    mh3 = mh2 + mh;
    for(k=1; k <= mh ; ++k){
        cosq = cos(PI_2 * k/len); // cos is even function
        e_arr[mh2-k].re = -cosq;
        e_arr[mh2+k].re = -cosq;
        e_arr[len-k].re =  cosq;
        e_arr[k    ].re =  cosq;

        cosq = (backward) ? cosq : -cosq; // make sin table with cos
        e_arr[mh3-k].im = -cosq;
        e_arr[mh3+k].im = -cosq;
        e_arr[mh-k ].im =  cosq;
        e_arr[mh+k ].im =  cosq;
    }
    e_arr[0  ] = _complex( 1, 0);
    e_arr[mh2] = _complex(-1, 0);
    e_arr[mh ] = (backward)? _complex(0, 1) : _complex( 0,-1);
    e_arr[mh3] = (backward)? _complex(0,-1) : _complex( 0, 1);
    //for(k=0;k < len; ++k)
    //    fprintf(stdout, "%d : %8.8f + %8.8f j\n", k, e_arr[k].re, e_arr[k].im);
#else
    double dir2PI = (backward) ? 2*M_PI : -2*M_PI;
#endif

    for (int32 ldm = ldn; ldm >= 1; ldm--, m >>= 2){ // m >>=2 -> m /=4
        //m = pow(4,ldm);
        mh  = m >> 2; // mh = m /4;
        mh2 = mh+mh;
        mh3 = mh2+mh;
        // extract loop for k = 0 from following loop
        for ( int32 r = 0; r < len ; r += m){
            u0 = dat[r    ];
            u1 = dat[r+mh ];
            u2 = dat[r+mh2];
            u3 = dat[r+mh3];

            x = u0 + u2;
            y = u1 + u3;
            t0 = x + y;
            t2 = x - y;

            x = u0 - u2;
            y = u1 - u3;
            y = (backward) ? _complex(-y.im, y.re) : _complex(y.im,-y.re); // y * j * dir ;
            t1 = x + y;
            t3 = x - y;

            dat[r    ] = t0;
            dat[r+mh ] = t1;
            dat[r+mh2] = t2;
            dat[r+mh3] = t3;
        }
        //for ( k=0; k < mh; k++){ // k : 0, 1, 2, ... , m/4-1
        for ( k=1; k < mh; k++){ // k : 0, 1, 2, ... , m/4-1
#ifdef FFT_EXP_TABLE
            int32 idx_exp = k << (ldn-ldm)*2;
            e = e_arr[idx_exp];
#else
            e.re = cos(dir2PI*k/m);
            e.im = sin(dir2PI*k/m);
#endif
            e2 = e * e;
            e3 = e2* e;
            for ( int32 r = 0; r < len ; r += m){
                u0 = dat[r+k    ];
                u1 = dat[r+k+mh ];
                u2 = dat[r+k+mh2];
                u3 = dat[r+k+mh3];

                x = u0 + u2;
                y = u1 + u3;
                t0 = x + y;
                t2 = x - y;

                x = u0 - u2;
                y = u1 - u3;
                y = (backward) ? _complex(-y.im, y.re) : _complex(y.im,-y.re); // y * j * dir ;
                t1 = x + y;
                t3 = x - y;

                t1 = t1 * e ;
                t2 = t2 * e2;
                t3 = t3 * e3;

                dat[r+k    ] = t0;
                dat[r+k+mh ] = t1;
                dat[r+k+mh2] = t2;
                dat[r+k+mh3] = t3;
            }
        }
    }
    permute_radix4( dat, len);
    if( backward ){
        for( k=0 ; k < len; ++k)
            dat[k] /= len;
    }

#ifdef FFT_EXP_TABLE
 #ifdef MALLOC_F
    free( e_arr);
 #else
    delete [] e_arr;
 #endif
#endif
}

void fft_compositN(_complex *dat, int32 len, std::vector<int32>& fac, bool backward){
/* reference : Inside The Fft Black Box- Serial And Parallel Fast Fourier Transform Algorithms, chapter 15*/

    int32 B,F,Q,BF,QB;

    _complex *buf = new _complex[len];
    _complex *Y, *C , *X;
    //printf("factors: \n");
    //for(int32 t=0; t< fac.size(); ++t)
    //    printf("%d ", fac.at(t));
    //printf("\n");

    double pi_x2 = M_PI * 2;
    double theta;
    _complex WBF, z, zf, Yp;
    int32 l, r, v;
    int32 qB, ffB;
    int32 ff, bb, q, f;

    B = 1;
    Y = dat; C = buf;
    for( v=0; v < fac.size() ; ++v){
        F = fac[v];
        BF= B*F;
        Q = len/BF;
        QB= B*Q;
        theta = pi_x2/BF;
        // exchange C and Y unconditionally
        X = C;
        C = Y;
        Y = X;
        WBF = backward ? _complex(cos(theta), sin(theta)) : _complex(cos(theta), -sin(theta));
        z = 1;
        for( ff=0, ffB=0; ff < F; ++ff, ffB+=B){
            for( bb=0; bb < B; ++bb){
                for( q=0, qB=0, r=ffB+bb; q < Q; ++q, qB+=B, r+=BF){
                    zf = 1;
                    Yp = 0;
                    for( f=0, l=qB+bb; f < F; ++f, l+=QB ){
                        //l = (f*Q+q)*B+bb; l = f*Q*B+q*B+bb
                        Yp += C[l]*zf;
                        zf *= z ;
                    }
                    //r = (q*F+ff)*B+bb; r = q*F*B+ff*B+bb
                    Y[r] = Yp;
                }
                z*=WBF;
            }
        }
        B *= F;
    }

    if( Y == dat){
        for(v=0; v < len; ++v)
            dat[v] = Y[v];
    }

    if( backward){
        for(v=0; v < len; ++v)
            dat[v] /= len;
    }

    delete [] buf;
}

void factorizeN(int32 N, std::vector<int32>& fac){
/* reference : https://cp-algorithms.com/algebra/factorization.html  */

    int32 somePrm[169] = {
               3,     5,     7,    11,    13,    17,    19,    23,    29,
       31,    37,    41,    43,    47,    53,    59,    61,    67,    71,
       73,    79,    83,    89,    97,   101,   103,   107,   109,   113,
      127,   131,   137,   139,   149,   151,   157,   163,   167,   173,
      179,   181,   191,   193,   197,   199,   211,   223,   227,   229,
      233,   239,   241,   251,   257,   263,   269,   271,   277,   281,
      283,   293,   307,   311,   313,   317,   331,   337,   347,   349,
      353,   359,   367,   373,   379,   383,   389,   397,   401,   409,
      419,   421,   431,   433,   439,   443,   449,   457,   461,   463,
      467,   479,   487,   491,   499,   503,   509,   521,   523,   541,
      547,   557,   563,   569,   571,   577,   587,   593,   599,   601,
      607,   613,   617,   619,   631,   641,   643,   647,   653,   659,
      661,   673,   677,   683,   691,   701,   709,   719,   727,   733,
      739,   743,   751,   757,   761,   769,   773,   787,   797,   809,
      811,   821,   823,   827,   829,   839,   853,   857,   859,   863,
      877,   881,   883,   887,   907,   911,   919,   929,   937,   941,
      947,   953,   967,   971,   977,   983,   991,   997,  1009,  1013
    };
    while( (N & 0x00000001) == 0){ // find factor 2
            fac.push_back(2);
            N >>= 1;
    }
    for (int32 d : somePrm ) {
        if(d*d > N) break;
        while(N % d == 0){
            fac.push_back(d);
            N /= d;
        }
    }
    for (int32 d = somePrm[168]+2; ; d += 2) {
        if(d*d > N) break;
        while (N % d == 0) {
            fac.push_back(d);
            N /= d;
        }
    }
    if (N > 1)
        fac.push_back(N);
}

float cubic1d(float a_1, float a0, float a1, float a2, float t){
    float A= a2 -a_1 + 3*(a0- a1);
    float B= (2*a_1+4*a1)-(5*a0+a2); //2*a_1 -5*a0 + 4*a1 -a2;
    float C= a1 -a_1 ;
    float D= 2*a0;
    float r= A*t*t*t + B*t*t + C*t + D;
    return r/2.0f;
}
Mat copy_padding(const Mat& src, int32 pad_size){
    int32 r  = src.getRow();
    int32 c  = src.getCol();
    int32 xR = r+pad_size*2;
    int32 xC = c+pad_size*2;
    Mat des(src.getDatType(),xR,xC,src.getChannel());
    matRect srcRec, desRec;
    int32 pd = pad_size-1;
    int32 i,k,m;
    int32 rm1 = r-1;
    int32 cm1 = c-1;
    //copy main data
    srcRec.set(0,0,rm1,cm1);
    desRec.set(pad_size,pad_size,r+pd,c+pd);
    Mat::sliceCopyMat(src,srcRec,des,desRec);

    for(i=0,k=r+pad_size,m=c+pad_size; i < pad_size; ++i){
        //top
        srcRec.set(i   ,0       ,i   ,cm1 );
        desRec.set(pd-i,pad_size,pd-i,c+pd);
        Mat::sliceCopyMat(src,srcRec,des,desRec);
        //left
        srcRec.set(0       ,i   ,rm1 ,i   );
        desRec.set(pad_size,pd-i,r+pd,pd-i);
        Mat::sliceCopyMat(src,srcRec,des,desRec);
        //bottom
        srcRec.set(rm1-i,0       ,rm1-i,cm1 );
        desRec.set(k+i  ,pad_size,k+i  ,c+pd);
        Mat::sliceCopyMat(src,srcRec,des,desRec);
        //right
        srcRec.set(0       ,cm1-i,rm1 ,cm1-i);
        desRec.set(pad_size,m+i  ,r+pd,m+i  );
        Mat::sliceCopyMat(src,srcRec,des,desRec);
    }
    for(i=0; i< pad_size; ++i){
        //top left corner
        srcRec.set(0 ,pad_size+i,pd ,pad_size+i);
        desRec.set(0 ,pd-i      ,pd ,pd-i      );
        Mat::sliceCopyMat(des,srcRec,des,desRec);
        //top right corner
        srcRec.set(0 ,c+pd-i      ,pd ,c+pd-i      );
        desRec.set(0 ,c+pad_size+i,pd ,c+pad_size+i);
        Mat::sliceCopyMat(des,srcRec,des,desRec);
        //bottom left corner
        srcRec.set(r+pad_size ,pad_size+i,xR-1 ,pad_size+i);
        desRec.set(r+pad_size ,pd-i      ,xR-1 ,pd-i      );
        Mat::sliceCopyMat(des,srcRec,des,desRec);
        //bottom right corner
        srcRec.set(r+pad_size ,c+pd-i      ,xR-1 ,c+pd-i      );
        desRec.set(r+pad_size ,c+pad_size+i,xR-1 ,c+pad_size+i);
        Mat::sliceCopyMat(des,srcRec,des,desRec);
    }
    return des;
}

Mat bicubicIntp(const Mat& src,const int32 s){

    switch(src.getDatType()){
    case DTYP::DOUBLE : return _bicubicIntp<double>(src, s);
    case DTYP::FLOAT  : return _bicubicIntp<float >(src, s);
    case DTYP::INT    : return _bicubicIntp<int32 >(src, s);
    case DTYP::UCHAR  : return _bicubicIntp<uchar >(src, s);
    default           : {
        fprintf(stderr, " Unsupported DTYP in bicubicIntp func.\n");
        return Mat();
        }
    }
}


} // end of imgproc namespace
} // end of jmat namespace
