#include "jbimgproc.h"

namespace jbimgproc {

jbMat rgb2ycc(const jbMat& rgbIm, const int sel_eq){
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
    jbMat A;
    switch(rgbIm.getDatType()){
    case DTYP::DOUBLE : A = _rgb2ycc<double>(rgbIm, sel_eq); break;
    case DTYP::FLOAT  : A = _rgb2ycc<float >(rgbIm, sel_eq); break;
    case DTYP::INT    : A = _rgb2ycc<int   >(rgbIm, sel_eq); break;
    case DTYP::UCHAR  : A = _rgb2ycc<uchar >(rgbIm, sel_eq); break;
    }

    return A;
}


jbMat ycc2rgb(const jbMat& rgbIm, const int sel_eq){
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

    jbMat A = rgbIm.copy();
     *rgbDat_pt = rgbIm.getMat().get();
    double *tarDat_pt = A.getMat().get();

    int ch_offset1 = imsize ;
    int ch_offset2 = imsize << 1;
    int x;
    double tmp1, tmp2, tmp3;
    if( sel_eq == 0){
        for(x= 0 ; x < imsize; x++){
            tmp1 = bt601_y2r[0][0] * rgbDat_pt[x  ] + bt601_y2r[0][1] * rgbDat_pt[x+ch_offset1] + bt601_y2r[0][2] * rgbDat_pt[x+ch_offset2];
            tmp2 = bt601_y2r[1][0] * rgbDat_pt[x  ] + bt601_y2r[1][1] * rgbDat_pt[x+ch_offset1] + bt601_y2r[1][2] * rgbDat_pt[x+ch_offset2];
            tmp3 = bt601_y2r[2][0] * rgbDat_pt[x  ] + bt601_y2r[2][1] * rgbDat_pt[x+ch_offset1] + bt601_y2r[2][2] * rgbDat_pt[x+ch_offset2];

            tarDat_pt[x           ] = (tmp1 < 0 ) ? 0 : tmp1;
            tarDat_pt[x+ch_offset1] = (tmp2 < 0 ) ? 0 : tmp2;
            tarDat_pt[x+ch_offset2] = (tmp3 < 0 ) ? 0 : tmp3;
        }
    }else if( sel_eq == 1){
        for(x= 0 ; x < imsize; x++){
            tmp1 = bt709_y2r[0][0] * rgbDat_pt[x  ] + bt709_y2r[0][1] * rgbDat_pt[x+ch_offset1] + bt709_y2r[0][2] * rgbDat_pt[x+ch_offset2];
            tmp2 = bt709_y2r[1][0] * rgbDat_pt[x  ] + bt709_y2r[1][1] * rgbDat_pt[x+ch_offset1] + bt709_y2r[1][2] * rgbDat_pt[x+ch_offset2];
            tmp3 = bt709_y2r[2][0] * rgbDat_pt[x  ] + bt709_y2r[2][1] * rgbDat_pt[x+ch_offset1] + bt709_y2r[2][2] * rgbDat_pt[x+ch_offset2];

            tarDat_pt[x           ] = (tmp1 < 0 ) ? 0 : tmp1;
            tarDat_pt[x+ch_offset1] = (tmp2 < 0 ) ? 0 : tmp2;
            tarDat_pt[x+ch_offset2] = (tmp3 < 0 ) ? 0 : tmp3;

        }
    }
    return A;
}

jbMat jbimgproc::rgb2gray(const jbMat& rgbIm, const int HowToGray){
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

jbMat jbimgproc::histoPmf(const jbMat& src){
    int ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stdout,"histoPmf : src argument is empty matrix\n");
        return jbMat();
    }else if( ch != 1) {
        fprintf(stdout,"histoPmf : src is not 1 channel matrix\n");
        return jbMat();
    }

    jbMat A = jbMat::zeros(1,256,1);
    double *tarDat_pt = A.getMat().get();
    double *srcDat_pt = src.getMat().get();

    int len = src.getLength();
    int k, d ;
    for(k=0; k < len ; k++){
        d = static_cast<int>(srcDat_pt[k]);
        d = (d > 255) ? 255 : (d < 0) ? 0 : d;
        tarDat_pt[d]++;
    }
    return A;
}

jbMat jbimgproc::histoCmf(const jbMat& src){
    int ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stdout,"histoCmf : src argument is empty matrix\n");
        return jbMat();
    }else if( ch != 1) {
        fprintf(stdout,"histoCmf : src is not 1 channel matrix\n");
        return jbMat();
    }

    jbMat cmf = jbImgproc::histoPmf(src);
    if( cmf.isEmpty()){
        fprintf(stdout,"histoPmf is empty matrix\n");
        return jbMat();
    }

    double *srcDat_pt = cmf.getMat().get();

    for(int k=1; k < 256 ; k++)
        srcDat_pt[k] += srcDat_pt[k-1];


    return cmf;
}

jbMat jbimgproc::clip_HistoPmf(const jbMat& src,const unsigned int clipVal){
    int ch  = src.getChannel();
    if( src.isEmpty() ){
        fprintf(stdout,"clip_histoPmf : src argument is empty matrix\n");
        return jbMat();
    }else if( ch != 1) {
        fprintf(stdout,"clip_histoPmf : src is not 1 channel matrix\n");
        return jbMat();
    }

    jbMat pmf = jbImgproc::histoPmf(src);
    if( pmf.isEmpty()){
        fprintf(stdout,"histoPmf is empty matrix\n");
        return jbMat();
    }

    double *srcDat_pt = pmf.getMat().get();

    // clipping
    unsigned int sum_clipped =0;
    unsigned int binval;
    for(int k=0; k < 256 ; k++){
        binval = (unsigned int)srcDat_pt[k];
        if( binval > clipVal){
            sum_clipped += binval - clipVal;
            srcDat_pt[k] = clipVal;
        }
    }
    sum_clipped >>= 8; // divided by 256
    // distributing the clipped sum
    for(int k=0; k < 256 ; k++){
        srcDat_pt[k] += sum_clipped;
    }

    return pmf;
}

jbMat jbimgproc::clip_HistoCmf(const jbMat& src,const unsigned int clipVal){
    int ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stdout,"histoCmf : src argument is empty matrix\n");
        return jbMat();
    }else if( ch != 1) {
        fprintf(stdout,"histoCmf : src is not 1 channel matrix\n");
        return jbMat();
    }

    jbMat cmf = jbImgproc::histoPmf(src);
    if( cmf.isEmpty()){
        fprintf(stdout,"histoPmf is empty matrix\n");
        return jbMat();
    }

    double *srcDat_pt = cmf.getMat().get();

    // clipping
    unsigned int sum_clipped =0;
    unsigned int binval;
    for(int k=0; k < 256 ; k++){
        binval = (unsigned int)srcDat_pt[k];
        if( binval > clipVal){
            sum_clipped += binval - clipVal;
            srcDat_pt[k] = clipVal;
        }
    }
    sum_clipped >>= 8; // divided by 256
    // distributing the clipped sum
    for(int k=0; k < 256 ; k++){
        srcDat_pt[k] += sum_clipped;
    }

    // making cumiltive data
    for(int k=1; k < 256 ; k++)
        srcDat_pt[k] += srcDat_pt[k-1];

    return cmf;
}

jbMat jbimgproc::clip_HistoEqual(const jbMat& src,const jbMat& histCmf){
    int ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stdout,"histoEqual : src argument is empty matrix\n");
        return jbMat();
    }else if( ch != 1) {
        fprintf(stdout,"histoEqual : src is not 1 channel matrix\n");
        return jbMat();
    }else if( histCmf.isEmpty()){
        fprintf(stdout,"histCmf is empty \n");
        return jbMat();
    }

    jbMat A = src.copy();

    double *srcDat_pt = src.getMat().get();
    double *mapDat_pt = histCmf.getMat().get();
    double *tarDat_pt = A.getMat().get();
    unsigned int dat;
    unsigned int binCnt = histCmf.getLength();
    for(int i=0; i < src.getLength(); i++){
        dat = static_cast<unsigned int>(srcDat_pt[i]);
        if(binCnt < dat )
            dat = binCnt;
        tarDat_pt[i] = mapDat_pt[dat];
    }

    return A;
}


}
