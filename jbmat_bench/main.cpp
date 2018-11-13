
#include "../jblib/jbMat.h"
#include "../jblib/jbmath.h"
#include "../jblib/qimmat.h"
//#include "../jblib/jbimgproc.h"

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

int main(int argc, char *argv[])
{

    uint row = 7;
    uint col = 7;
    uint ch = 2;
    uint rowcol = row*col;
    jbMat ma(DTYP::DOUBLE, row,col,ch);

    double* a = ma.getDataPtr<double>();
    uint k = 1;
    for(uint c=0; c < ch; c++){
        for(uint i=0; i < row*col; i++){
            a[c*rowcol + i] = k++;
        }
    }

    ma.printMat();
    ma += 10;
    ma.printMat();
    jbMat maa = ma.copy();
    //maa.plusMat(ma);
    maa += ma;
    maa.printMat();
    ma.at<double>(0) = 1000;
    ma.printMat();
    ma.changeDType(DTYP::INT);
    ma.printMat();
    ma.at<int>(2,2,1) = 2000;
    ma.printMat();
    jbMat mab = ma*ma;
    mab.printMat();
    jbMat mac = ma/ma;
    mac.printMat();
    jbMat mad = ma -maa;
    mad.printMat();
    jbMat mae(DTYP::DOUBLE, 5, 5, 2);

    for(uint c = 0; c< mae.getLength(); c++)
        mae.at<double>(c) = rand()%100;
    mae.printMat();

    jbMat maf = jbmath::inverse(mae);
    maf.printMat("inversed mat");

    jbMat mag(DTYP::DOUBLE, 3,4,2);
    jbMat mah(DTYP::DOUBLE, 4,3,2);
    double *mag_m = mag.getDataPtr<double>();
    double *mah_m = mah.getDataPtr<double>();
    for(uint c=0; c < mag.getLength(); c++){
        mag_m[c] = c;
        mah_m[c] = c + c%4;
    }
    jbMat mai = jbmath::dot(mag,mah);
    jbMat maj = jbmath::dot(mah,mag);
    mag.printMat();
    mah.printMat();
    mai.printMat();
    maj.printMat();

    //-- conv
    k=0;
    jbMat mb(DTYP::DOUBLE,5,5,2,"mb");
    jbMat mc(DTYP::DOUBLE,5,5,2,"mc");
    double* b = mb.getDataPtr<double>();
    for(uint c=0; c < mb.getChannel(); c++){
        for(uint i=0; i < 5; i++){
            for (uint j=0; j < 5; j++){
                b[i*5+j+c*25] = 100+k++;
                mc.at<double>(i,j,c) = k;
            }
        }
    }
    ma.printMat("ma");
    mb.printMat("mb");
    fprintf(stdout,"1) full\n");
    jbMat convO = jbmath::conv2d(ma,mb,"zero","full");
    convO.printMat("conv0: full");
    fprintf(stdout,"2) same\n");
    convO = jbmath::conv2d(ma,mb,"zero","same");
    convO.printMat("conv0: same");


    jbMat me = ma.copy();
    jbMat mf = ma;
    fprintf(stdout," me_mA_ptr: %p, md_mA_ptr: %p , mf_mA_ptr: %p\n", me.getMat().get(), ma.getMat().get(), mf.getMat().get());
    me.printMat();
    mf.printMat();
    ma.printMat();

    qDebug() << QDir::currentPath();
    QString fname1 = QString("../jbmat_bench/test.jpg");
    QString fname2 = QString("../jbmat_bench/test_filt.bmp");
    QImage img(fname1);
    jbMat matIm = qimmat::qim2jbmat(img);
    jbMat filt  = jbMat::ones(5,5,matIm.getChannel())/25;
    jbMat FiltIm = jbmath::conv2d(matIm, filt,"zero","same");
    //jbMat FiltIm = jbMath::conv2d(matIm, filt,"zero","full");
    QImage cvim = qimmat::jbmat2qim(FiltIm);
    //QImage cvim = QimMat::jbmat2qim(matIm);
    cvim.save(fname2);

    jbMat mg(DTYP::DOUBLE,1,7,3) ;
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
/*
    mg.printMat();
    jbMat mh = jbimgproc::rgb2ycc(mg,1);
    mh.printMat();
    jbMat mi = jbimgproc::ycc2rgb(mh,1);
    mi.printMat();
    jbMat mj = jbimgproc::rgb2gray(mg);
    mj.printMat();

    jbMat mk(3,4,2,"mk");
    jbMat ml(4,3,2,"ml");
    jbMat mn(4,4,1,"mn");
    for(int i=0; i < 24; i++){
        mk[i] = rand() % 50;
        ml[i] = rand() % 80;
        if(i < 16)
            mn[i] = rand() % 100;
    }
    jbMat mm = jbMath::mulMatrix(mk,ml);
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
    jbMat Y  = jbimgproc::rgb2gray(matIm);
    jbMat yccIm = jbimgproc::rgb2ycc(matIm);
    jbMat ms = jbimgproc::clip_HistoCmf(Y, 1000);
    ms /= ms[255];
    ms *= 255;
    jbMat mv = jbimgproc::clip_HistoEqual(Y,ms);
    yccIm.setChannelN(mv,0,1,0);
    jbMat mvv;
    mvv.setName("mvv");
    mvv.setChannelN(Y,0,1,0);
    //for(int i=0; i < Y.getRow()*Y.getCol(); i++)
    //    yccIm[i] = mv[i] ;
    jbMat histEqIm = jbimgproc::ycc2rgb(yccIm);
    QImage cvim2 = QimMat::jbmat2qim(histEqIm);
    //QImage cvim = QimMat::jbmat2qim(matIm);
    cvim2.save(QString("../jbmat_bench/test_histEq.bmp"));
    QImage cvim3 = QimMat::jbmat2qim(mvv);
    cvim3.save(QString("../jbmat_bench/test_y.bmp"));
    jbMat mx = matIm.copySubMat(1,4,100,104);
    mx.setName("mx_subMat");
    mx.printMat();
    jbMat my = matIm.copyChannelN(0);
    jbMat mz = my.copySubMat(1,4,100,104);
    mz.printMat();
*/
    return 0;
}
