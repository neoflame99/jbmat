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


Mat histoPmf(const Mat& src, const uint32 bins, const double step, const double lowLmtVal){
    uint32 ch_src  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stderr,"histoPmf : src argument is empty matrix\n");
        return Mat();
    }else if( bins < 1) {
        fprintf(stderr,"histoPmf : 'bins' should be larger than or equal 1 \n");
        return Mat();
    }

    Mat A = Mat::zeros(1, bins, ch_src, DTYP::DOUBLE);
    uint32 ch = A.getChannel();
    elemptr tarDat_pt = A.getElptr();
    elemptr srcDat_pt = src.getElptr();

    uint32 k, m, d ;

    if(src.getDatType()==DTYP::DOUBLE){
        for(m=0 ; m < ch ; ++m){
            for(k=m ; k < src.getLength() ; k+=ch ){
                d = uint32((srcDat_pt.f64_ptr[k]-lowLmtVal)/step);
                d = (d < 0) ? 0 : (d >= bins) ? bins-1: d;  // bin index range : 0 ~ bins-1
                tarDat_pt.f64_ptr[d*ch+m]++;
            }
        }
    }else if(src.getDatType()==DTYP::FLOAT){
        for(m=0 ; m < ch ; ++m){
            for(k=m ; k < src.getLength() ; k+=ch ){
                d = uint32((srcDat_pt.f32_ptr[k]-lowLmtVal)/step);
                d = (d < 0) ? 0 : (d >= bins) ? bins-1: d;  // bin index range : 0 ~ bins-1
                tarDat_pt.f64_ptr[d*ch+m]++;
            }
        }
    }else if(src.getDatType()==DTYP::INT){
        for(m=0 ; m < ch ; ++m){
            for(k=m ; k < src.getLength() ; k+=ch ){
                d = uint32((srcDat_pt.int_ptr[k]-lowLmtVal)/step);
                d = (d < 0) ? 0 : (d >= bins) ? bins-1: d;  // bin index range : 0 ~ bins-1
                tarDat_pt.f64_ptr[d*ch+m]++;
            }
        }
    }else if(src.getDatType()==DTYP::UCHAR){
        for(m=0 ; m < ch ; ++m){
            for(k=m ; k < src.getLength() ; k+=ch ){
                d = uint32((srcDat_pt.uch_ptr[k]-lowLmtVal)/step);
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


Mat histoCmf(const Mat& src, const uint32 bins, const double step, const double lowLmtVal){
    if( src.isEmpty() ){
        fprintf(stderr,"histoCmf : src argument is empty matrix\n");
        return Mat();
    }else if( bins < 1) {
        fprintf(stderr,"histoCmf : 'bins' should be larger than or equal 1 \n");
        return Mat();
    }else if( src.getDatType()==DTYP::CMPLX){
        fprintf(stderr, " Unsupported DTYP in histoCmf func.\n");
        return Mat();
    }

    Mat cmf = histoPmf(src, bins, step, lowLmtVal);
    if( cmf.isEmpty()) return cmf;

    double *srcDat_pt = cmf.getDataPtr<double>();
    uint32 ch_cmf = cmf.getChannel();
    uint32 m, k;
    for( m=0; m < ch_cmf; ++m){
        for( k=ch_cmf+m ; k < cmf.getLength() ; k+= ch_cmf)  // getLength() == bins*ch
            srcDat_pt[k] += srcDat_pt[k-ch_cmf];
    }

    return cmf;
}

Mat clip_HistoPmf(const Mat& src,const uint32 clipVal,const uint32 bins, const double step, const double lowLmtVal){
    if( src.isEmpty() ){
        fprintf(stderr,"clip_histoPmf : src argument is empty matrix\n");
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

    Mat pmf = histoPmf(src, bins, step, lowLmtVal );
    if( pmf.isEmpty()) return pmf;
    double *srcDat_pt = pmf.getDataPtr<double>();

    uint32 ch_pmf = pmf.getChannel();
    uint32 sum_clipped, quot, rem ;
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
        quot = sum_clipped / bins;
        rem  = sum_clipped % bins;
        // distributing the clipped sum
        for( k=m ; k < pmf.getLength() ; k+= ch_pmf){
            srcDat_pt[k] += quot;
        }
        // distributing remainders to higher numbered bins
        for( k=(bins-rem)*ch_pmf+m ; k < pmf.getLength() ; k+= ch_pmf){
            srcDat_pt[k]++;
        }
    }
    return pmf;
}


Mat clip_HistoCmf(const Mat& src,const uint32 clipVal,const uint32 bins, const double step, const double lowLmtVal, const bool normalz){
    if( src.isEmpty() ){
        fprintf(stderr,"clip_histoCmf : src argument is empty matrix\n");
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

    Mat cmf = clip_HistoPmf(src, clipVal, bins, step, lowLmtVal);
    if(cmf.isEmpty()) return cmf;

    double *srcDat_pt = cmf.getDataPtr<double>();
    uint32 ch_cmf = cmf.getChannel();
    uint32 m, k;

    // making cumiltive data
    for( m = 0 ; m < ch_cmf ; ++m){
        for(k=ch_cmf+m ; k < cmf.getLength() ; k+= ch_cmf) // getLength() == bins*ch
            srcDat_pt[k] += srcDat_pt[k-ch_cmf];
    }
    // normalizing
    if(normalz){
        double denom;
        uint32 cmf_size = cmf.getLength();
        for( m = 0 ; m < ch_cmf ; ++m){
            denom = srcDat_pt[cmf_size-ch_cmf+m];
            for(k=m ; k < cmf_size-ch_cmf ; k+= ch_cmf)
                srcDat_pt[k] /= denom;
            srcDat_pt[cmf_size-ch_cmf+m]=1.0;
        }
    }
    return cmf;
}


Mat clip_HistoEq(const Mat& src,const Mat& histCmf, const double step, const double maxval, const double minval){
    uint32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stdout,"histoEqual : src argument is empty matrix\n");
        return Mat();
    }else if( histCmf.isEmpty()){
        fprintf(stdout,"histCmf is empty \n");
        return Mat();
    }else if( src.getDatType()==DTYP::CMPLX ){
        fprintf(stderr, " Unsupported DTYP in clip_HistoEqual func.\n");
        return Mat();
    }else if( histCmf.getDatType()!=DTYP::DOUBLE){
        fprintf(stderr, " HistoCmf's data type of clip_HistoEq is to be DTYP::DOUBLE. \n");
        return Mat();
    }

    Mat A = src.copy();
    elemptr srcDat_pt = src.getElptr();
    elemptr tarDat_pt = A.getElptr();

    Mat cmf     = histCmf.copy();
    uint32 ch_cmf = cmf.getChannel();
    uint32 cmf_sz = cmf.getLength();
    uint32 c, i;
    double* mapDat_pt = cmf.getDataPtr<double>();
    double leastbin_val;
    double scale_rat;
    // stratching cmf to fit full range
    for( c=0; c < ch_cmf; c++ ){
        leastbin_val = mapDat_pt[c];
        for(i=c; i < cmf_sz; i += ch_cmf ){
            mapDat_pt[i] -= leastbin_val;
        }
        scale_rat = (maxval-minval)/mapDat_pt[cmf_sz-ch+c];
        for(i=c; i < cmf_sz; i += ch_cmf ){
            mapDat_pt[i] *= scale_rat;
        }
    }

    uint32 bins = cmf.getRowColSize();
    double d0, d1, d2;
    int32  i0, i1;
    double mp1, mp2, mv;
    double lowlmt = 0;
    double upplmt = (bins-1)*step;
    uint32 ch_leastbin ;
    uint32 ch_mostbin ;
    if(src.getDatType()==DTYP::DOUBLE){
        for( c=0; c < ch; ++c){
            ch_leastbin = c;
            ch_mostbin  = (bins-1)*ch +c;
            for( i=c; i < src.getLength(); i+=ch){
                d0 = srcDat_pt.f64_ptr[i];
                if( d0 <= lowlmt ){
                    tarDat_pt.f64_ptr[i] = mapDat_pt[ch_leastbin];
                }else if( d0 >= upplmt){
                    tarDat_pt.f64_ptr[i] = mapDat_pt[ch_mostbin];
                }else {
                    d1 = floor(d0 / step);
                    d2 = d0 - d1*step;

                    i0 = d1*ch+c;
                    i1 = i0+ch;
                    mp1 = mapDat_pt[i0];
                    mp2 = mapDat_pt[i1];

                    mv  = (mp2-mp1)*d2/step;
                    tarDat_pt.f64_ptr[i] = (mp1 + mv);
                }
            }
        }
    }else if(src.getDatType()==DTYP::FLOAT){
        for( c=0; c < ch; ++c){
            ch_leastbin = c;
            ch_mostbin  = (bins-1)*ch +c;
            for( i=c; i < src.getLength(); i+=ch){
                d0 = srcDat_pt.f32_ptr[i];
                if( d0 <= lowlmt ){
                    tarDat_pt.f32_ptr[i] = mapDat_pt[ch_leastbin];
                }else if( d0 >= upplmt){
                    tarDat_pt.f32_ptr[i] = mapDat_pt[ch_mostbin];
                }else {
                    d1 = floor(d0 / step);
                    d2 = d0 - d1*step;

                    i0 = d1*ch+c;
                    i1 = i0+ch;
                    mp1 = mapDat_pt[i0];
                    mp2 = mapDat_pt[i1];

                    mv  = (mp2-mp1)*d2/step;
                    tarDat_pt.f32_ptr[i] = (mp1 + mv);
                }
            }
        }
    }else if(src.getDatType()==DTYP::INT){
        for( c=0; c < ch; ++c){
            ch_leastbin = c;
            ch_mostbin  = (bins-1)*ch +c;
            for( i=c; i < src.getLength(); i+=ch){
                d0 = srcDat_pt.int_ptr[i];
                if( d0 <= lowlmt ){
                    tarDat_pt.int_ptr[i] = mapDat_pt[ch_leastbin];
                }else if( d0 >= upplmt){
                    tarDat_pt.int_ptr[i] = mapDat_pt[ch_mostbin];
                }else {
                    d1 = floor(d0 / step);
                    d2 = d0 - d1*step;

                    i0 = d1*ch+c;
                    i1 = i0+ch;
                    mp1 = mapDat_pt[i0];
                    mp2 = mapDat_pt[i1];

                    mv  = (mp2-mp1)*d2/step;
                    tarDat_pt.int_ptr[i] = (mp1 + mv);
                }
            }
        }
    }else if(src.getDatType()==DTYP::UCHAR){
        for( c=0; c < ch; ++c){
            ch_leastbin = c;
            ch_mostbin  = (bins-1)*ch +c;
            for( i=c; i < src.getLength(); i+=ch){
                d0 = srcDat_pt.uch_ptr[i];
                if( d0 <= lowlmt ){
                    tarDat_pt.uch_ptr[i] = mapDat_pt[ch_leastbin];
                }else if( d0 >= upplmt){
                    tarDat_pt.uch_ptr[i] = mapDat_pt[ch_mostbin];
                }else {
                    d1 = floor(d0 / step);
                    d2 = d0 - d1*step;

                    i0 = d1*ch+c;
                    i1 = i0+ch;
                    mp1 = mapDat_pt[i0];
                    mp2 = mapDat_pt[i1];

                    mv  = (mp2-mp1)*d2/step;
                    tarDat_pt.uch_ptr[i] = (mp1 + mv);
                }
            }
        }
    }
    return A;
}

bool gamma( const Mat& src, const double gmvalue, Mat& des){
    if(src.getDatType() != des.getDatType()){
        fprintf(stderr,"src Mat and des Mat are not equal types!\n") ;
        return false;
    }else if(src.getDatType()==DTYP::CMPLX){
        fprintf(stderr, "Unsuppretd data type in gamma func.\n");
        return false;
    }else if(src.isEmpty()){
        fprintf(stderr,"src Mat is empty!\n") ;
        return false;
    }else if(des.isEmpty()){
        fprintf(stderr,"des Mat is empty!\n") ;
        return false;
    }

    DTYP srcDtype = src.getDatType();

    double maxA = src.max().max().at<double>(0);
    Mat mx = src.max();
    double dval;
    elemptr src_ptr = src.getElptr();
    elemptr des_ptr = des.getElptr();

    uint32 len = src.getLength();
    if( srcDtype== DTYP::DOUBLE){
        for(uint32 i=0; i < len; ++i){
            dval = src_ptr.f64_ptr[i] / maxA;
            dval = pow( dval, gmvalue) * maxA;
            des_ptr.f64_ptr[i] = dval;
        }
    }else if(srcDtype== DTYP::FLOAT){
        for(uint32 i=0; i < len; ++i){
            dval = src_ptr.f32_ptr[i] / maxA;
            dval = pow( dval, gmvalue) * maxA;
            des_ptr.f32_ptr[i] = (float) dval;
        }
    }else if(srcDtype== DTYP::INT){
        for(uint32 i=0; i < len; ++i){
            dval = src_ptr.int_ptr[i] / maxA;
            dval = pow( dval, gmvalue) * maxA;
            des_ptr.int_ptr[i] = (int32) dval;
        }
    }else if(srcDtype== DTYP::UCHAR){
        for(uint32 i=0; i < len; ++i){
            dval = src_ptr.uch_ptr[i] / maxA;
            dval = pow( dval, gmvalue) * maxA;
            des_ptr.uch_ptr[i] = (uchar) dval;
        }
    }

    return true;
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
