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
/* hdrfile reading and writing to bmp */
#ifdef _MACOS_
//QString hdrfile_path = QString("/Users/neoflame99/HDRI/memorial.hdr");
static QString hdrfile_path = QString("/Users/neoflame99/HDRI/office.hdr");
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

    retinexSigmoidMethod1(hdrimg);

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

    Yhdrimg   = imgproc::rgb2gray(hdrimg);
    Mat stdv = Mat::zeros(20,1,1, DTYP::DOUBLE);
    Mat stdvIdx = Mat::zeros(20,1,1,DTYP::INT);
    double maxv;
    uint32 idx;
    QString rfn;
    double gmv = 0.45;
    /* ----- 1st stage tonemapping ------ */
    Mat facV = Mat::zeros(20,1,1,DTYP::DOUBLE);
    for (i = 0, fac = 0.1 ; i < 20; ++i, fac += 0.1 ){
        facV.at<double>(static_cast<uint32>(i)) = fac;
        Yhdrimg_tm1  = imgproc::nakaSigTonemap(Yhdrimg, gau, fac);
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
    fac = 0;
    Yhdrimg_tm1 = imgproc::nakaSigTonemap(Yhdrimg, gau, fac);


    Yhdrimg_tm3c = Mat::repeat(Yhdrimg_tm1, 1, 1, 3);
    Yhdrimg_3c   = Mat::repeat(Yhdrimg, 1, 1, 3);

    hdrimg_tm = imgproc::gamma( hdrimg * Yhdrimg_tm3c / Yhdrimg_3c, gmv);
    hdrimg_tm = hdrimg_tm / hdrimg_tm.max().max().at<double>(0)*255.0;
    hdrtm_bmp = qimmat::mat2qim(hdrimg_tm);
    rfn = QString("../output/%1_tm1_%2_gm_%3.bmp").arg(Filename_only).arg(fac).arg(gmv);
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

}
