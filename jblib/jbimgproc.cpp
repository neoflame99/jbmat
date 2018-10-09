#include "jbimgproc.h"


const double bt601_r2y[3][3]={{ 0.299    ,  0.587    ,  0.114 },
                              {-0.168736 , -0.331264 ,  0.5   },
                              { 0.5      , -0.418688 , -0.081312} };
const double bt709_r2y[3][3]={{ 0.2126  ,  0.7152  , 0.0722  },
                              {-0.11457 , -0.38543 , 0.5     },
                              { 0.5     , -0.45415 ,-0.04585 }};

const double bt601_y2r[3][3]={{ 1.0000 ,  -0.0000 ,   1.4020 },
                              { 1.0000 ,  -0.3441 ,  -0.7141 },
                              { 1.0000 ,   1.7720 ,   0.0000 }};

const double bt709_y2r[3][3]={{ 1.0000 ,   0.0000 ,   1.5748 },
                              { 1.0000 ,  -0.1873 ,  -0.4681 },
                              { 1.0000 ,   1.8556 ,  -0.0000 }};

jbImgproc::jbImgproc()
{

}

jbMat jbImgproc::rgb2ycc(const jbMat& rgbIm, const int sel_eq){
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

    int row = rgbIm.getRow();
    int col = rgbIm.getCol();
    int imsize = row * col;
    int chsize = rgbIm.getChannel();
    if( chsize != 3 ) {
        fprintf(stdout,"rgbIm is not three channel image\n");
        return jbMat();
    }

    jbMat  A(row, col, chsize);
    double *rgbDat_pt = rgbIm.getMat().get();
    double *tarDat_pt = A.getMat().get();

    int ch_offset1 = imsize;
    int ch_offset2 = imsize << 1;

    int x;
    if( sel_eq == 0){
        for(x= 0 ; x < imsize; x++ ){
            tarDat_pt[x           ] = bt601_r2y[0][0] * rgbDat_pt[x  ] + bt601_r2y[0][1] * rgbDat_pt[x+ch_offset1] + bt601_r2y[0][2] * rgbDat_pt[x+ch_offset2];
            tarDat_pt[x+ch_offset1] = bt601_r2y[1][0] * rgbDat_pt[x  ] + bt601_r2y[1][1] * rgbDat_pt[x+ch_offset1] + bt601_r2y[1][2] * rgbDat_pt[x+ch_offset2];
            tarDat_pt[x+ch_offset2] = bt601_r2y[2][0] * rgbDat_pt[x  ] + bt601_r2y[2][1] * rgbDat_pt[x+ch_offset1] + bt601_r2y[2][2] * rgbDat_pt[x+ch_offset2];
        }
    }else if( sel_eq == 1){
        for(x= 0 ; x < imsize; x++ ){
            tarDat_pt[x           ] = bt709_r2y[0][0] * rgbDat_pt[x  ] + bt709_r2y[0][1] * rgbDat_pt[x+ch_offset1] + bt709_r2y[0][2] * rgbDat_pt[x+ch_offset2];
            tarDat_pt[x+ch_offset1] = bt709_r2y[1][0] * rgbDat_pt[x  ] + bt709_r2y[1][1] * rgbDat_pt[x+ch_offset1] + bt709_r2y[1][2] * rgbDat_pt[x+ch_offset2];
            tarDat_pt[x+ch_offset2] = bt709_r2y[2][0] * rgbDat_pt[x  ] + bt709_r2y[2][1] * rgbDat_pt[x+ch_offset1] + bt709_r2y[2][2] * rgbDat_pt[x+ch_offset2];
        }
    }
    return A;
}


jbMat jbImgproc::ycc2rgb(const jbMat& rgbIm, const int sel_eq){
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

    int row = rgbIm.getRow();
    int col = rgbIm.getCol();
    int imsize = row * col;
    int chsize = rgbIm.getChannel();
    if( chsize != 3 ) {
        fprintf(stdout,"rgbIm is not three channel image\n");
        return jbMat();
    }

    jbMat A(row, col, chsize);
    double *rgbDat_pt = rgbIm.getMat().get();
    double *tarDat_pt = A.getMat().get();

    int ch_offset1 = imsize ;
    int ch_offset2 = imsize << 1;
    int x;
    if( sel_eq == 0){
        for(x= 0 ; x < imsize; x++){
            tarDat_pt[x           ] = bt601_y2r[0][0] * rgbDat_pt[x  ] + bt601_y2r[0][1] * rgbDat_pt[x+ch_offset1] + bt601_y2r[0][2] * rgbDat_pt[x+ch_offset2];
            tarDat_pt[x+ch_offset1] = bt601_y2r[1][0] * rgbDat_pt[x  ] + bt601_y2r[1][1] * rgbDat_pt[x+ch_offset1] + bt601_y2r[1][2] * rgbDat_pt[x+ch_offset2];
            tarDat_pt[x+ch_offset2] = bt601_y2r[2][0] * rgbDat_pt[x  ] + bt601_y2r[2][1] * rgbDat_pt[x+ch_offset1] + bt601_y2r[2][2] * rgbDat_pt[x+ch_offset2];
        }
    }else if( sel_eq == 1){
        for(x= 0 ; x < imsize; x++){
            tarDat_pt[x           ] = bt709_y2r[0][0] * rgbDat_pt[x  ] + bt709_y2r[0][1] * rgbDat_pt[x+ch_offset1] + bt709_y2r[0][2] * rgbDat_pt[x+ch_offset2];
            tarDat_pt[x+ch_offset1] = bt709_y2r[1][0] * rgbDat_pt[x  ] + bt709_y2r[1][1] * rgbDat_pt[x+ch_offset1] + bt709_y2r[1][2] * rgbDat_pt[x+ch_offset2];
            tarDat_pt[x+ch_offset2] = bt709_y2r[2][0] * rgbDat_pt[x  ] + bt709_y2r[2][1] * rgbDat_pt[x+ch_offset1] + bt709_y2r[2][2] * rgbDat_pt[x+ch_offset2];
        }
    }
    return A;
}

jbMat jbImgproc::rgb2gray(const jbMat& rgbIm, const int HowToGray){
/*
 *  if HowToGray = 0 (BT 601)
 *  Y =  0.299 * r   +  0.587 * g +  0.114 * b
 *
 *  if HowToGray = 1 (BT 709)
 *  Y =  0.2126 * r +  0.7152 * g +  0.0722 * b
 *
 *  if HowToGray = 2 (3 equal-weight)
 *  Y = 0.333 * r + 0.334 * g + 0.333 * b
 */

    int row = rgbIm.getRow();
    int col = rgbIm.getCol();
    int imsize = row * col;
    int chsize = rgbIm.getChannel();
    if( chsize != 3 ) {
        fprintf(stdout,"rgbIm is not three channel image\n");
        return jbMat();
    }

    jbMat A(row, col, 1);
    double *rgbDat_pt = rgbIm.getMat().get();
    double *tarDat_pt = A.getMat().get();

    int ch_offset1 = imsize;
    int ch_offset2 = imsize << 1;
    int x,k;

    if( HowToGray==0){
        for(x= 0, k=0 ; x < imsize; x++, k++)
            tarDat_pt[k  ] =  bt601_r2y[0][0] * rgbDat_pt[x  ] + bt601_r2y[0][1] * rgbDat_pt[x+ch_offset1] + bt601_r2y[0][2] * rgbDat_pt[x+ch_offset2];
    }else if( HowToGray==1){
        for(x= 0, k=0 ; x < imsize; x++, k++)
            tarDat_pt[k  ] =  bt709_r2y[0][0] * rgbDat_pt[x  ] + bt709_r2y[0][1] * rgbDat_pt[x+ch_offset1] + bt709_r2y[0][2] * rgbDat_pt[x+ch_offset2];
    }else if( HowToGray==2){
        for(x= 0, k=0 ; x < imsize; x++, k++)
            tarDat_pt[k  ] =  0.333 * rgbDat_pt[x  ] + 0.334 * rgbDat_pt[x+ch_offset1] + 0.333 * rgbDat_pt[x+ch_offset2];
    }

    return A;
}

jbMat jbImgproc::histoPmf(const jbMat& src){
    int ch  = src.getChannel();
//    int row = src.getRow();
//    int col = src.getCol();
    if( src.isEmpty() ){
        fprintf(stdout,"histoPmf : src argument is empty matrix\n");
        return jbMat();
    }else if( ch != 1) {
        fprintf(stdout,"histoPmf : src is not 1 channel matrix\n");
        return jbMat();
    }

    jbMat A = jbMat::ones(1,255,1);
    double *tarDat_pt = A.getMat().get();
    double *srcDat_pt = src.getMat().get();

    int len = src.getLength();
    int k, d ;
    for(k=0; k < len ; k++){
        d = static_cast<int>(srcDat_pt[k]);
        d = (d > 255) ? 255 : (d < 0) ? 0 : d;
        (tarDat_pt[d])++;
    }
    return A;
}
