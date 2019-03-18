
#include "../jblib/jbMat.h"
#include "../jblib/jbmath.h"
#include "../jblib/qimmat.h"
#include "../jblib/jbimgproc.h"

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

int32 main(int32 argc, char *argv[])
{

    uint32 row = 7;
    uint32 col = 7;
    uint32 ch = 2;
    uint32 rowcol = row*col;
    Mat ma(DTYP::DOUBLE, row,col,ch,"ma");

    double* a = ma.getDataPtr<double>();
    uint32 k = 1;
    for(uint32 c=0; c < ch; c++){
        for(uint32 i=0; i < row*col; i++){
            a[c*rowcol + i] = k++;
        }
    }

    ma.printMat();
    ma = ma + 10;
    ma.setName("ma");
    //ma += 10;
    ma.printMat();
    Mat maa = ma.copy();
    //maa.plusMat(ma);
    maa += ma;
    maa.printMat();
    ma.at<double>(0) = 1000;
    ma.printMat();
    ma.changeDType(DTYP::INT);
    ma.printMat();
    ma.at<int32>(2,2,1) = 2000;
    ma.printMat();
    Mat mab = ma*ma;
    mab.printMat();
    Mat mac = ma/ma;
    mac.printMat();
    Mat mad = ma -maa;
    mad.printMat();
    Mat mae(DTYP::DOUBLE, 5, 5, 2);

    for(uint32 c = 0; c< mae.getLength(); c++)
        mae.at<double>(c) = rand()%100;
    mae.printMat();

    Mat maf = inverse(mae);
    maf.printMat("inversed mat");

    Mat mag(DTYP::DOUBLE, 3,4,2,"mag");
    Mat mah(DTYP::DOUBLE, 4,3,2,"mah");
    Mat man(DTYP::DOUBLE, 3,4,2,"man");
    double *mag_m = mag.getDataPtr<double>();
    double *mah_m = mah.getDataPtr<double>();
    for(uint32 c=0; c < mag.getLength(); c++){
        mag_m[c] = c;
        mah_m[c] = c + c%4;
        man.at<double>(c) = rand()%1000;
    }
    Mat mai = mul(mag,mah);
    Mat maj = mul(mah,mag);
    mai.setName("mai");
    mag.printMat();
    mah.printMat();
    mai.printMat();
    maj.printMat();
    Mat mak = mai.mean<double>();
    Mat mal = mai.std<double>();
    mak.printMat("mak : mean of mai");
    mal.printMat("mal : std of mai");

    man.printMat();
    man.max<double>().printMat("man max");
    man.min<double>().printMat("man min");


    //-- conv
    k=0;
    Mat mb(DTYP::DOUBLE,5,5,2,"mb");
    Mat mc(DTYP::DOUBLE,5,5,2,"mc");
    double* b = mb.getDataPtr<double>();
    for(uint32 c=0; c < mb.getChannel(); c++){
        for(uint32 i=0; i < 5; i++){
            for (uint32 j=0; j < 5; j++){
                b[i*5+j+c*25] = 100+k++;
                mc.at<double>(i,j,c) = k;
            }
        }
    }
    ma.printMat("ma");
    mb.printMat("mb");
    fprintf(stdout,"1) full\n");
    Mat convO = conv2d(ma,mb,"zero","full");
    convO.printMat("conv0: full");
    fprintf(stdout,"2) same\n");
    convO = conv2d(ma,mb,"zero","same");
    convO.printMat("conv0: same");

    Mat me = ma.copy();
    Mat mf = ma;
    fprintf(stdout," me_mA_ptr: %p, md_mA_ptr: %p , mf_mA_ptr: %p\n", me.getMat().get(), ma.getMat().get(), mf.getMat().get());
    me.printMat();
    mf.printMat();
    ma.printMat();

    qDebug() << QDir::currentPath();
    QString fname1 = QString("../jbmat_bench/test.jpg");
    QString fname2 = QString("../jbmat_bench/test_filt.bmp");
    QImage img(fname1);
    Mat matIm = qimmat::qim2mat(img);
    Mat filt  =  Mat::ones(5,5,matIm.getChannel())/25;
    Mat FiltIm = conv2d(matIm, filt,"zero","same");
    // Mat FiltIm = jbMath::conv2d(matIm, filt,"zero","full");
    QImage cvim = qimmat::mat2qim(FiltIm);
    //QImage cvim = QimMat::jbmat2qim(matIm);
    cvim.save(fname2);

    /* hdrfile reading and writing to bmp */

    QString hdrfile_path = QString("/Users/neoflame99/Workspace/Qt5/readhdrfile/readhdrfile/memorial.hdr");
    Mat hdrimg = qimmat::read_hdr(hdrfile_path)*128;

    QImage hdr_bmp = qimmat::mat2qim(hdrimg);
    hdr_bmp.save("../jbmat_bench/hdr_bmp.bmp");
    /* ----- */

    Mat mg(DTYP::DOUBLE,1,7,3) ;
    mg.at<double>(0,0,0) = 1;
    mg.at<double>(0,0,1) = 1;
    mg.at<double>(0,0,2) = 1;
    mg.at<double>(0,1,0) = 1;
    mg.at<double>(0,1,1) = 0;
    mg.at<double>(0,1,2) = 0;
    mg.at<double>(0,2,0) = 0;
    mg.at<double>(0,2,1) = 1;
    mg.at<double>(0,2,2) = 0;
    mg.at<double>(0,3,0) = 0;
    mg.at<double>(0,3,1) = 0;
    mg.at<double>(0,3,2) = 1;
    mg.at<double>(0,4,0) = 1;
    mg.at<double>(0,4,1) = 0;
    mg.at<double>(0,4,2) = 1;
    mg.at<double>(0,5,0) = 1;
    mg.at<double>(0,5,1) = 1;
    mg.at<double>(0,5,2) = 0;
    mg.at<double>(0,6,0) = 0;
    mg.at<double>(0,6,1) = 1;
    mg.at<double>(0,6,2) = 1;

    mg.printMat();
    Mat mh = imgproc::rgb2ycc(mg,1);
    mh.printMat();
    Mat mi = imgproc::ycc2rgb(mh,1);
    mi.printMat();
    Mat mj = imgproc::rgb2gray(mg);
    mj.printMat();

    /*dot product */

    Mat mdps1({1,2,3});
    Mat mdps2({4,5,6});
    //mdps2.transpose();
    Mat mdp1 = dot(mdps1, mdps2);
    mdp1.printMat("vector dot product");

    Mat mdps3({1,2,3,4,5,6});
    Mat mdps4({7,8,9,10,11,12});
    mdps3.reshape(3,2,1);
    mdps4.reshape(3,2,1);
    Mat mdp2 = dot(mdps3, mdps4, 0);
    Mat mdp3 = dot(mdps3, mdps4, 1);
    mdps3.printMat("mdps3");
    mdps4.printMat("mdps4");
    mdp2.printMat("dim=0 dot product");
    mdp3.printMat("dim=1 dot product");


    Mat tt;
    {
        tt= Mat(DTYP::INT,1,1,1);
        tt.printMat();
    }

    /*--- gauss mask / box mask gen ---*/

    Mat gau = imgproc::gaussMaskGen(1,4);
    Mat box = imgproc::boxMaskGen(5);

    gau.printMat("gaussian Mat");
    box.printMat("box Mat");

    /*------------------*/
    /*
    Mat mk(DTYP::DOUBLE,3,4,2,"mk");
    Mat ml(DTYP::DOUBLE,4,3,2,"ml");
    Mat mn(DTYP::DOUBLE,4,4,1,"mn");
    for(int32 i=0; i < 24; i++){
        mk[i] = rand() % 50;
        ml[i] = rand() % 80;
        if(i < 16)
            mn[i] = rand() % 100;
    }
     Mat mm = mulMatrix(mk,ml);
    mk.printMat(); //std::string("mk"));
    ml.printMat(); //std::string("ml"));
    mm.printMat(); //std::string("mm"));
    mm.transpose();
    mm.printMat(std::string("mm_transposed"));
    mk.reshape(12,2,1);
    mk.printMat(std::string("mk reshape"));
    mn.printMat(std::string("mn"));
    //mn = jbMath::inverse(mn);
    //mn.printMat(std::string("mn inverse"));
     Mat Y  = imgproc::rgb2gray(matIm);
     Mat yccIm = imgproc::rgb2ycc(matIm);
     Mat ms = imgproc::clip_HistoCmf(Y, 1000);
    ms /= ms[255];
    ms *= 255;
     Mat mv = imgproc::clip_HistoEqual(Y,ms);
    yccIm.setChannelN(mv,0,1,0);
     Mat mvv;
    mvv.setName("mvv");
    mvv.setChannelN(Y,0,1,0);
    //for(int32 i=0; i < Y.getRow()*Y.getCol(); i++)
    //    yccIm[i] = mv[i] ;
     Mat histEqIm = imgproc::ycc2rgb(yccIm);
    QImage cvim2 = QimMat::mat2qim(histEqIm);
    //QImage cvim = QimMat::jbmat2qim(matIm);
    cvim2.save(QString("../jbmat_bench/test_histEq.bmp"));
    QImage cvim3 = QimMat::mat2qim(mvv);
    cvim3.save(QString("../jbmat_bench/test_y.bmp"));
     Mat mx = matIm.copySubMat(1,4,100,104);
    mx.setName("mx_subMat");
    mx.printMat();
     Mat my = matIm.copyChannelN(0);
     Mat mz = my.copySubMat(1,4,100,104);
    mz.printMat();
*/
    return 0;
}
