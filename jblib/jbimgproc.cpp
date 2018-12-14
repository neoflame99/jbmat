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

    switch(rgbIm.getDatType()){
    case DTYP::DOUBLE : return _rgb2ycc<double>(rgbIm, sel_eq);
    case DTYP::FLOAT  : return _rgb2ycc<float >(rgbIm, sel_eq);
    case DTYP::INT    : return _rgb2ycc<int32 >(rgbIm, sel_eq);
    case DTYP::UCHAR  : return _rgb2ycc<uchar >(rgbIm, sel_eq);
    }
}


Mat ycc2rgb(const Mat& rgbIm, const int32 sel_eq){
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

    switch(rgbIm.getDatType()){
    case DTYP::DOUBLE : return _ycc2rgb<double>(rgbIm, sel_eq);
    case DTYP::FLOAT  : return _ycc2rgb<float >(rgbIm, sel_eq);
    case DTYP::INT    : return _ycc2rgb<int32 >(rgbIm, sel_eq);
    case DTYP::UCHAR  : return _ycc2rgb<uchar >(rgbIm, sel_eq);
    }
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
 *  Y = 0.333 * r + 0.334 * g + 0.333 * b
 */
    switch(rgbIm.getDatType()){
    case DTYP::DOUBLE : return _rgb2gray<double>(rgbIm, HowToGray);
    case DTYP::FLOAT  : return _rgb2gray<float >(rgbIm, HowToGray);
    case DTYP::INT    : return _rgb2gray<int32 >(rgbIm, HowToGray);
    case DTYP::UCHAR  : return _rgb2gray<uchar >(rgbIm, HowToGray);
    }
}
/*
Mat jbimgproc::histoPmf(const Mat& src){
    int32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stdout,"histoPmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stdout,"histoPmf : src is not 1 channel matrix\n");
        return Mat();
    }

    Mat A = Mat::zeros(1,256,1);
    double *tarDat_pt = A.getMat().get();
    double *srcDat_pt = src.getMat().get();

    int32 len = src.getLength();
    int32 k, d ;
    for(k=0; k < len ; k++){
        d = static_cast<int32>(srcDat_pt[k]);
        d = (d > 255) ? 255 : (d < 0) ? 0 : d;
        tarDat_pt[d]++;
    }
    return A;
}

Mat jbimgproc::histoCmf(const Mat& src){
    int32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stdout,"histoCmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stdout,"histoCmf : src is not 1 channel matrix\n");
        return Mat();
    }

    Mat cmf = jbImgproc::histoPmf(src);
    if( cmf.isEmpty()){
        fprintf(stdout,"histoPmf is empty matrix\n");
        return Mat();
    }

    double *srcDat_pt = cmf.getMat().get();

    for(int32 k=1; k < 256 ; k++)
        srcDat_pt[k] += srcDat_pt[k-1];


    return cmf;
}

Mat jbimgproc::clip_HistoPmf(const Mat& src,const unsigned int32 clipVal){
    int32 ch  = src.getChannel();
    if( src.isEmpty() ){
        fprintf(stdout,"clip_histoPmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stdout,"clip_histoPmf : src is not 1 channel matrix\n");
        return Mat();
    }

    Mat pmf = jbImgproc::histoPmf(src);
    if( pmf.isEmpty()){
        fprintf(stdout,"histoPmf is empty matrix\n");
        return Mat();
    }

    double *srcDat_pt = pmf.getMat().get();

    // clipping
    unsigned int32 sum_clipped =0;
    unsigned int32 binval;
    for(int32 k=0; k < 256 ; k++){
        binval = (unsigned int32)srcDat_pt[k];
        if( binval > clipVal){
            sum_clipped += binval - clipVal;
            srcDat_pt[k] = clipVal;
        }
    }
    sum_clipped >>= 8; // divided by 256
    // distributing the clipped sum
    for(int32 k=0; k < 256 ; k++){
        srcDat_pt[k] += sum_clipped;
    }

    return pmf;
}

Mat jbimgproc::clip_HistoCmf(const Mat& src,const unsigned int32 clipVal){
    int32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stdout,"histoCmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stdout,"histoCmf : src is not 1 channel matrix\n");
        return Mat();
    }

    Mat cmf = jbImgproc::histoPmf(src);
    if( cmf.isEmpty()){
        fprintf(stdout,"histoPmf is empty matrix\n");
        return Mat();
    }

    double *srcDat_pt = cmf.getMat().get();

    // clipping
    unsigned int32 sum_clipped =0;
    unsigned int32 binval;
    for(int32 k=0; k < 256 ; k++){
        binval = (unsigned int32)srcDat_pt[k];
        if( binval > clipVal){
            sum_clipped += binval - clipVal;
            srcDat_pt[k] = clipVal;
        }
    }
    sum_clipped >>= 8; // divided by 256
    // distributing the clipped sum
    for(int32 k=0; k < 256 ; k++){
        srcDat_pt[k] += sum_clipped;
    }

    // making cumiltive data
    for(int32 k=1; k < 256 ; k++)
        srcDat_pt[k] += srcDat_pt[k-1];

    return cmf;
}

Mat jbimgproc::clip_HistoEqual(const Mat& src,const Mat& histCmf){
    int32 ch  = src.getChannel();

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

    Mat A = src.copy();

    double *srcDat_pt = src.getMat().get();
    double *mapDat_pt = histCmf.getMat().get();
    double *tarDat_pt = A.getMat().get();
    unsigned int32 dat;
    unsigned int32 binCnt = histCmf.getLength();
    for(int32 i=0; i < src.getLength(); i++){
        dat = static_cast<unsigned int32>(srcDat_pt[i]);
        if(binCnt < dat )
            dat = binCnt;
        tarDat_pt[i] = mapDat_pt[dat];
    }

    return A;
}
*/

} // end of imgproc namespace
} // end of jmat namespace
