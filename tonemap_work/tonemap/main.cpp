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

using namespace jmat;

Mat retinexSigmoidMethod1(Mat);
Mat retinexSigmoidMethod2(Mat);
Mat retinexLogMethod(Mat);

/* hdrfile reading and writing to bmp */
#ifdef _MACOS_
static QString hdrfile_path = QString("/Users/neoflame99/HDRI/memorial.hdr");
//static QString hdrfile_path = QString("/Users/neoflame99/HDRI/office.hdr");
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
    Mat hdrimg = qimmat::read_hdr(hdrfile_path);
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
    double gmv = 0.45;

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
    double gmv = 0.45;

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
    int32 boxmasksize = 6;
    Mat gau = imgproc::gaussMaskGen(1,6);
    Mat box = imgproc::boxMaskGen(boxmasksize);

    //gau.printMat("gaussian Mat");
    //box.printMat("box Mat");

    /*------------------*/

    /*----- tonemapping ------ */
    Mat Yhdrimg, Yhdrimg_tm;
    Mat Yhdrimg_plus1;
    Mat Yhdrimg_tm3c, Yhdrimg_3c, hdrimg_tm;
    QImage hdrtm_bmp;
    Mat surround;

    Yhdrimg       = imgproc::rgb2gray(hdrimg,2);
    Yhdrimg_plus1 = Yhdrimg+1;

    double maxv;
    QString rfn;
    double gmv = 0.45;

    maxv = hdrimg.max().max().at<double>(0);

    // * ---- Intensity log tonemapping ---- * //
    surround   = imgproc::localMeanMat(Yhdrimg, gau);
    Yhdrimg_tm = imgproc::logRetinexTonemap(Yhdrimg_plus1, surround );
    double tm_max = Yhdrimg_tm.max().at<double>(0);
    double tm_min = Yhdrimg_tm.min().at<double>(0);

    //Yhdrimg_tm = (Yhdrimg_tm - tm_min) * (255.0 / (tm_max - tm_min));

    /* ----- color channel processing ------ */
    double scaleB = 255.0 / maxv;

    uint32 r  = hdrimg.getRow();
    uint32 c  = hdrimg.getCol();
    uint32 ch = hdrimg.getChannel();
    uint32 rc = hdrimg.getRowColSize();
    Mat ratioMat  = Yhdrimg_tm / Yhdrimg;
    Mat ratioMat2(DTYP::DOUBLE, r,c,ch);

    double ratioVal ;
    uint32 idx=0;
    for( uint32 k =0 ; k < ch; ++k){
        for ( uint32 m=0 ; m < rc; ++m){
            ratioVal = ratioMat.at<double>(m);
            ratioMat2.at<double>(idx++) = ratioVal; // ( scaleB < ratioVal ) ? scaleB : ratioVal;
        }
    }
    fprintf(stdout, "%f %f", Yhdrimg_tm.max().at<double>(), Yhdrimg_tm.min().at<double>());
    hdrimg_tm = imgproc::gamma( hdrimg * ratioMat2, gmv);
    //hdrimg_tm = hdrimg_tm / hdrimg_tm.max().max().at<double>(0)*255.0;
    hdrtm_bmp = qimmat::mat2qim(hdrimg_tm);
    rfn = QString("../output/%1_logmethod_gm_%2.bmp").arg(Filename_only).arg(gmv);
    hdrtm_bmp.save(rfn);

    return hdrimg_tm;
}
