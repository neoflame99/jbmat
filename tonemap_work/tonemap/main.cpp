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

int main(int argc, char *argv[])
{
    /* hdrfile reading and writing to bmp */
#ifdef _MACOS_
    QString hdrfile_path = QString("/Users/neoflame99/Workspace/Qt5/readhdrfile/readhdrfile/memorial.hdr");
#else
    QString hdrfile_path = QString("/home/neoflame99/memorial.hdr");
#endif
    Mat hdrimg = qimmat::read_hdr(hdrfile_path);
    QImage hdr_bmp = qimmat::mat2qim(hdrimg);
    hdr_bmp.save("../output/hdr_bmp.bmp");
    /* ----- */

    /*--- gauss mask / box mask gen ---*/

    Mat gau = imgproc::gaussMaskGen(1,4);
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
    for (i =0, fac = 0.2 ; i < 16; ++i, fac += 0.1 ){
        Yhdrimg      = imgproc::rgb2gray(hdrimg);
        Yhdrimg_tm1  = imgproc::nakaSigTonemap(Yhdrimg, gau, 0.5);
        Yhdrimg_tm2  = imgproc::nakaSigTonemap(Yhdrimg_tm1, gau, fac);

        Yhdrimg_tm3c = Mat::repeat(Yhdrimg_tm2, 1, 1, 3);
        Yhdrimg_3c   = Mat::repeat(Yhdrimg, 1, 1, 3);

        hdrimg_tm = imgproc::gamma( hdrimg * Yhdrimg_tm3c / Yhdrimg_3c, 0.5);
        hdrimg_tm = hdrimg_tm / hdrimg_tm.max().max().at<double>(0)*255.0;
        hdrtm_bmp = qimmat::mat2qim(hdrimg_tm);
        QString rfn = QString("../output/hdr_tm_%1.bmp").arg(fac);
        hdrtm_bmp.save(rfn);
        fprintf(stdout,"%d, %f\n",i, fac);
    }

    return 0;
}
