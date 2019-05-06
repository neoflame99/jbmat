#include "../../jblib/jbMat.h"
#include "../../jblib/jbmath.h"
#include "../../jblib/qimmat.h"
#include "../../jblib/jbimgproc.h"

#include <QImage>
#include <QString>
#include <QDir>
#include <QDebug>
#include <stdlib.h>
#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

#define _HDRIMG_

using namespace jmat;

Mat retinexSigmoidMethod1(Mat);
Mat retinexSigmoidMethod2(Mat);
Mat retinexLogMethod(Mat);


/* hdrfile reading and writing to bmp */
#ifdef _MACOS_
#ifdef _HDRIMG_
static QString hdrfile_path = QString("/Users/neoflame99/HDRI/memorial.hdr");
//static QString hdrfile_path = QString("/Users/neoflame99/HDRI/office.hdr");
#else
static QString hdrfile_path = QString("/Users/neoflame99/Workspace/Experiments/MSRetinex/example.bmp");
#endif
static QFileInfo f(hdrfile_path);
static QString Filename_only = f.fileName();
#else
//QString hdrfile_path = QString("/home/neoflame99/HDRI/memorial.hdr");
static QString hdrfile_path = QString("/home/neoflame99/HDRI/office.hdr");
static QFileInfo f(hdrfile_path);
static QString Filename_only = f.fileName();
#endif

int main(int argc, char *argv[])
{
#ifdef _HDRIMG_
    Mat hdrimg = qimmat::read_hdr(hdrfile_path);
#else
    QImage im(hdrfile_path);
    Mat hdrimg = qimmat::qim2mat(im);
#endif
    //QImage hdr_bmp = qimmat::mat2qim(hdrimg);
    //hdr_bmp.save("../output/hdr_bmp.bmp");


    fprintf(stdout," retinex method1\n");
    retinexSigmoidMethod1(hdrimg);
    fprintf(stdout," retinex method2\n");
    retinexSigmoidMethod2(hdrimg);

    fprintf(stdout," retinex log method \n");
    retinexLogMethod(hdrimg);
    return 0;
}

Mat retinexSigmoidMethod1(Mat hdrimg){

    /*--- gauss mask / box mask gen ---*/

    Mat gau = imgproc::gaussMaskGen(1,6);
    Mat box = imgproc::boxMaskGen(5);

    gau.printMat("gaussian Mat");
    box.printMat("box Mat");

    /*------------------*/

    /*----- tonemapping ------ */
    int i;
    double fac;
    Mat Yhdrimg, Yhdrimg_tm1, Yhdrimg_tm2;
    Mat Yhdrimg_tm3c, Yhdrimg_3c, hdrimg_tm;
    QImage hdrtm_bmp;


    double Imax;
    double globalMean;
    Yhdrimg   = imgproc::rgb2gray(hdrimg);
    Mat stdv = Mat::zeros(20,1,1, DTYP::DOUBLE);
    Mat stdvIdx = Mat::zeros(20,1,1,DTYP::INT);
    double maxv;
    uint32 idx=0;
    QString rfn;
    double gmv = 0.7; //0.45;

    Imax = Yhdrimg.max().at<double>(0);
    globalMean = Yhdrimg.mean().at<double>(0);
    Mat localMean = imgproc::localMeanMat(Yhdrimg, gau)*2;
    /* ----- 1st stage tonemapping ------ */
    Mat facV = Mat::zeros(20,1,1,DTYP::DOUBLE);
    for (i = 0, fac = 0.1 ; i < 20; ++i, fac += 0.1 ){
        facV.at<double>(static_cast<uint32>(i)) = fac;
        Yhdrimg_tm1  = imgproc::nakaSigTonemap(Yhdrimg, localMean, fac*globalMean, Imax);
        stdv.at<double>(static_cast<uint32>(i)) = Yhdrimg_tm1.std().at<double>(0);
    }
    maxv = stdv.max().at<double>(0);
    for (i=0; i < 20; ++i){
        if ( stdv.at<double>(static_cast<uint32>(i)) == maxv ){
            stdvIdx.at<int32>(static_cast<uint32>(i)) = 1;
            idx = static_cast<uint32>(i);
        }
    }
    stdv.printMat("stdv 1st:");
    stdvIdx.printMat("stdvIdx 1st: ");
    fac = facV.at<double>(idx);
    Yhdrimg_tm1 = imgproc::nakaSigTonemap(Yhdrimg, localMean, fac*globalMean, Imax);


    Yhdrimg_tm3c = Mat::repeat(Yhdrimg_tm1, 1, 1, 3);
    Yhdrimg_3c   = Mat::repeat(Yhdrimg, 1, 1, 3);

    hdrimg_tm = imgproc::gamma( hdrimg * Yhdrimg_tm3c / Yhdrimg_3c, gmv);
    hdrimg_tm = hdrimg_tm / hdrimg_tm.max().max().at<double>(0)*255.0;
    hdrtm_bmp = qimmat::mat2qim(hdrimg_tm);
    rfn = QString("../output/%1_method1_tm1_%2_gm_%3.bmp").arg(Filename_only).arg(fac).arg(gmv);
    hdrtm_bmp.save(rfn);

    /* ------- 2nd stage tonemapping ------- */
    /*
    facV = Mat::zeros(20,1,1,DTYP::DOUBLE);
    for (i = 0, fac = 0.1 ; i < 20; ++i, fac += 0.1 ){
        facV.at<double>(static_cast<uint32>(i)) = fac;
        Yhdrimg_tm2  = imgproc::nakaSigTonemap(Yhdrimg_tm1, gau, fac);
        stdv.at<double>(static_cast<uint32>(i)) = Yhdrimg_tm2.std().at<double>(0);
    }
    maxv = stdv.max().at<double>(0);
    for (i=0; i < 20; ++i){
        if ( stdv.at<double>(static_cast<uint32>(i)) == maxv ){
            stdvIdx.at<int32>(static_cast<uint32>(i)) = 1;
            idx = static_cast<uint32>(i);
        }
    }
    stdvIdx.printMat("stdvIdx 2nd: ");
    fac = facV.at<double>(1);
    Yhdrimg_tm2 = imgproc::nakaSigTonemap(Yhdrimg_tm1, gau, fac);

    Yhdrimg_tm3c = Mat::repeat(Yhdrimg_tm2, 1, 1, 3);

    hdrimg_tm = imgproc::gamma( hdrimg * Yhdrimg_tm3c / Yhdrimg_3c, 1.0);
    hdrimg_tm = hdrimg_tm / hdrimg_tm.max().max().at<double>(0)*255.0;
    hdrtm_bmp = qimmat::mat2qim(hdrimg_tm);
    rfn = QString("../output/%1_tm2_%2_gm_%3.bmp").arg(Filename_only).arg(fac).arg(1.0);
    hdrtm_bmp.save(rfn);
*/
    return hdrimg_tm;
}
Mat retinexSigmoidMethod2(Mat hdrimg){

    /*--- gauss mask / box mask gen ---*/
    int32 boxmasksize = 6;
    Mat gau = imgproc::gaussMaskGen(1,6);
    Mat box = imgproc::boxMaskGen(boxmasksize);

    gau.printMat("gaussian Mat");
    box.printMat("box Mat");

    /*------------------*/

    /*----- tonemapping ------ */
    int i;
    double fac;
    Mat Yhdrimg, Yhdrimg_tm1, Yhdrimg_tm2;
    Mat Yhdrimg_tm3c, Yhdrimg_3c, hdrimg_tm;
    QImage hdrtm_bmp;


    double Imax;
    double globalMean;
    Yhdrimg   = imgproc::rgb2gray(hdrimg);
    Mat stdv = Mat::zeros(20,1,1, DTYP::DOUBLE);
    Mat stdvIdx = Mat::zeros(20,1,1,DTYP::INT);
    double maxv;
    uint32 idx=0;
    QString rfn;
    double gmv = 0.7; //0.45;

    Imax = Yhdrimg.max().at<double>(0);
    globalMean = Yhdrimg.mean().at<double>(0);
    Mat s_localMean = imgproc::localMeanMat(Yhdrimg, gau);
    Mat l_localMean = imgproc::localMeanMat(Yhdrimg, box);
    /* ----- 1st stage tonemapping ------ */
    Mat facV = Mat::zeros(20,1,1,DTYP::DOUBLE);
    for (i = 0, fac = 0.1 ; i < 20; ++i, fac += 0.1 ){
        facV.at<double>(static_cast<uint32>(i)) = fac;
        Yhdrimg_tm1  = imgproc::nakaSig3MeanTonemap(Yhdrimg, s_localMean, l_localMean, fac*globalMean, Imax);
        stdv.at<double>(static_cast<uint32>(i)) = Yhdrimg_tm1.std().at<double>(0);
    }
    maxv = stdv.max().at<double>(0);
    for (i=0; i < 20; ++i){
        if ( stdv.at<double>(static_cast<uint32>(i)) == maxv ){
            stdvIdx.at<int32>(static_cast<uint32>(i)) = 1;
            idx = static_cast<uint32>(i);
        }
    }
    stdv.printMat("stdv 1st:");
    stdvIdx.printMat("stdvIdx 1st: ");
    fac = facV.at<double>(idx);
    Yhdrimg_tm1 = imgproc::nakaSig3MeanTonemap(Yhdrimg, s_localMean, l_localMean, fac*globalMean, Imax);


    Yhdrimg_tm3c = Mat::repeat(Yhdrimg_tm1, 1, 1, 3);
    Yhdrimg_3c   = Mat::repeat(Yhdrimg, 1, 1, 3);

    hdrimg_tm = imgproc::gamma( hdrimg * Yhdrimg_tm3c / Yhdrimg_3c, gmv);
    hdrimg_tm = hdrimg_tm / hdrimg_tm.max().max().at<double>(0)*255.0;
    hdrtm_bmp = qimmat::mat2qim(hdrimg_tm);
    rfn = QString("../output/%1_method2_gls_tm1_%2_gm_%3_llm_%4.bmp").arg(Filename_only).arg(fac).arg(gmv).arg(boxmasksize);
    hdrtm_bmp.save(rfn);

    return hdrimg_tm;
}

Mat retinexLogMethod(Mat hdrimg){

    /*--- gauss mask / box mask gen ---*/
    int32 boxmasksize = 300;
    Mat gau = imgproc::gaussMaskGen(2,6);
    Mat box = imgproc::boxMaskGen(boxmasksize);

    gau.printMat("gaussian Mat");
    //box.printMat("box Mat");

    /*------------------*/

    /*----- tonemapping ------ */
    Mat Yhdrimg, Yhdrimg_tm;
    Mat hdrimg_tm = hdrimg.copy();
    QImage hdrtm_bmp;
    Mat surround;

    /* extract gray */
    hdrimg += 1;
    Yhdrimg       = imgproc::rgb2gray(hdrimg,2);

    double maxv;
    QString rfn;
    double gmv = 0.45;


    // * ---- Intensity log tonemapping ---- * //
    surround   = imgproc::localMeanMat(Yhdrimg, gau);
    hdrtm_bmp  = qimmat::mat2qim(surround);
    hdrtm_bmp.save("../output/surround.bmp");
    Yhdrimg_tm = imgproc::logRetinexTonemap(Yhdrimg, surround );

    /* ---- simplest color balance by histogram ---- */
    uint32 len = Yhdrimg_tm.getLength(); // One channel Matrix
    double pl = 1; // low pixels of 1% among total pixels
    double ph = 1; // high piexls of 1% among total pixels
    double lo_pixels = pl*0.01*len;
    double hi_pixels = len - ph*0.01*len;
    double chistval;
    uint32 k;
    uint32 binIdx1=0, binIdx2=0;
    double scale;

    double *yhdrimgTm_datPt = Yhdrimg_tm.getDataPtr<double>();
    double *yhdrimg_datPt   = Yhdrimg.getDataPtr<double>();
    double *hdrimg_datPt    = hdrimg.getDataPtr<double>();
    double *hdrimgTm_datPt  = hdrimg_tm.getDataPtr<double>();
    double binstep = 0.01;
    Mat histCm = imgproc::histoCmf(Yhdrimg_tm,256,binstep);

    for( k=0; k < histCm.getRowColSize(); k++){
        chistval = histCm.at<double>(k);
        if( chistval <  lo_pixels) binIdx1 = k;
        if( chistval <= hi_pixels) binIdx2 = k;
    }

    scale = 255. /((double)(binIdx2 - binIdx1)*binstep);

    for( k = 0; k < len ; k++){
        if( yhdrimgTm_datPt[k] < binIdx1)      yhdrimgTm_datPt[k] = 0.;
        else if( yhdrimgTm_datPt[k] > binIdx2) yhdrimgTm_datPt[k] = 255.;
        else yhdrimgTm_datPt[k] = scale*(yhdrimgTm_datPt[k] - binIdx1);
    }

    hdrtm_bmp = qimmat::mat2qim(Yhdrimg_tm);
    hdrtm_bmp.save("../output/logTmGray.bmp");

    /* ----- color channel processing ------ */
    double ro,go,bo, r,g,b;
    uint32 len2 = len << 1;
    for(k =0; k < len ; k++){
        if(yhdrimg_datPt[k] <= 1.) yhdrimg_datPt[k]=1.;
        scale = yhdrimgTm_datPt[k]/ yhdrimg_datPt[k];
        if( scale > 3.) scale = 3.;

        r  = hdrimg_datPt[k     ];
        g  = hdrimg_datPt[k+len ];
        b  = hdrimg_datPt[k+len2];
        ro = scale*r;
        go = scale*g;
        bo = scale*b;

        if( ro > 255. || go > 255. || bo > 255.)
        {
            maxv = r;
            if( g > maxv) maxv = g;
            if( b > maxv) maxv = b;
            scale = 255. /maxv;
            ro = scale*r;
            go = scale*g;
            bo = scale*b;
        }
        hdrimgTm_datPt[k     ] = ro;
        hdrimgTm_datPt[k+len ] = go;
        hdrimgTm_datPt[k+len2] = bo;
    }


    hdrtm_bmp = qimmat::mat2qim(hdrimg_tm);
    rfn = QString("../output/%1_logmethod_gm_%2.bmp").arg(Filename_only).arg(gmv);
    hdrtm_bmp.save(rfn);

    return hdrimg_tm;
}
