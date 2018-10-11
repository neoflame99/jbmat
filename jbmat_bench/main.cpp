//#include <QCoreApplication>
#include "../jblib/jbMat.h"
#include "../jblib/jbmath.h"
#include "../jblib/qimmat.h"
#include "../jblib/jbimgproc.h"

#include <QImage>
#include <QString>
#include <QDir>
#include <QDebug>
#include <stdlib.h>

int main(int argc, char *argv[])
{

//    QCoreApplication a(argc, argv);
//    return a.exec();
    int row = 7;
    int col = 7;
    int ch = 2;
    int rowcol = row*col;
    jbMat ma(row,col,ch);
    //double* a = ma.getMat();
    double* a = ma.getMat().get();
    int k = 1;
    for(int c=0; c < ch; c++){
        for(int i=0; i < row*col; i++){
            a[c*rowcol + i] = k++;
        }
    }

    k=0;
    ma.printMat();
    jbMat mb(5,5,2);
    jbMat mc(5,5,2);
    //double *b = mb.getMat();
    double* b = mb.getMat().get();
    for(int c=0; c < mb.getChannel(); c++){
        for(int i=0; i < 5; i++){
            for (int j=0; j < 5; j++){
                b[i*5+j+c*25] = 100+k++;
                mc(i,j,c) = k;
            }
        }
    }
    mb.printMat();
    fprintf(stdout,"1) full\n");
    jbMat convO = jbMath::conv2d(ma,mb,"zero","full");
    convO.printMat();
    fprintf(stdout,"2) same\n");
    convO = jbMath::conv2d(ma,mb,"zero","same");
    convO.printMat();

    jbMat md = mb/mc;
    mb.printMat();
    mc.printMat();
    md.printMat();

    jbMat me = md.copy();
    jbMat mf = md;
    fprintf(stdout," me_mA_ptr: %p, md_mA_ptr: %p , mf_mA_ptr: %p\n", me.getMat().get(), md.getMat().get(), mf.getMat().get());
    me.printMat();
    mf.printMat();
    md.printMat();
/*
    for(int k=1; k<=3;k++ ){
        QString fname1 = QString("../%1.bmp").arg(k);
        QString fname2 = QString("../%1_1.bmp").arg(k);
        QImage img(fname1);
        jbMat matIm = QimMat::qim2jbmat(img);
        QImage cvim = QimMat::jbmat2qim(matIm);
        cvim.save(fname2);
    }
*/

    qDebug() << QDir::currentPath();

    QString fname1 = QString("../jbmat_bench/test.jpg");
    QString fname2 = QString("../jbmat_bench/test_filt.bmp");
    QImage img(fname1);
    jbMat matIm = QimMat::qim2jbmat(img);
    jbMat filt  = jbMat::ones(5,5,matIm.getChannel())/25;
    jbMat FiltIm = jbMath::conv2d(matIm, filt,"zero","same");
    //jbMat FiltIm = jbMath::conv2d(matIm, filt,"zero","full");
    QImage cvim = QimMat::jbmat2qim(FiltIm);
    //QImage cvim = QimMat::jbmat2qim(matIm);
    cvim.save(fname2);

    jbMat mg(1,7,3) ;
    mg(0,0,0) = 1;
    mg(0,0,1) = 1;
    mg(0,0,2) = 1;
    mg(0,1,0) = 1;
    mg(0,1,1) = 0;
    mg(0,1,2) = 0;
    mg(0,2,0) = 0;
    mg(0,2,1) = 1;
    mg(0,2,2) = 0;
    mg(0,3,0) = 0;
    mg(0,3,1) = 0;
    mg(0,3,2) = 1;
    mg(0,4,0) = 1;
    mg(0,4,1) = 0;
    mg(0,4,2) = 1;
    mg(0,5,0) = 1;
    mg(0,5,1) = 1;
    mg(0,5,2) = 0;
    mg(0,6,0) = 0;
    mg(0,6,1) = 1;
    mg(0,6,2) = 1;

    mg.printMat();
    jbMat mh = jbImgproc::rgb2ycc(mg,1);
    mh.printMat();
    jbMat mi = jbImgproc::ycc2rgb(mh,1);
    mi.printMat();
    jbMat mj = jbImgproc::rgb2gray(mg);
    mj.printMat();
/**/
    jbMat mk(3,4,2);
    jbMat ml(4,3,2);
    jbMat mn(4,4,1);
    for(int i=0; i < 24; i++){
        mk[i] = rand() % 50;
        ml[i] = rand() % 80;
        if(i < 16)
            mn[i] = rand() % 100;
    }
    jbMat mm = jbMath::mulMatrix(mk,ml);
    mk.printMat(std::string("mk"));
    ml.printMat(std::string("ml"));
    mm.printMat(std::string("mm"));
    mm.transpose();
    mm.printMat(std::string("mm_transposed"));
    mk.reshape(12,2,1);
    mk.printMat(std::string("mk reshape"));
    mn.printMat(std::string("mn"));
    //mn = jbMath::inverse(mn);
    //mn.printMat(std::string("mn inverse"));
    jbMat Y  = jbImgproc::rgb2gray(matIm);
    jbMat yccIm = jbImgproc::rgb2ycc(matIm);
    jbMat ms = jbImgproc::clip_HistoCmf(Y, 1000);

    ms /= ms[255];
    ms *= 255;
    jbMat mv = jbImgproc::clip_HistoEqual(Y,ms);
    for(int i=0; i < Y.getRow()*Y.getCol(); i++)
        yccIm[i] = mv[i] ;
    jbMat histEqIm = jbImgproc::ycc2rgb(yccIm);
    QImage cvim2 = QimMat::jbmat2qim(histEqIm);
    //QImage cvim = QimMat::jbmat2qim(matIm);
    cvim2.save(QString("../jbmat_bench/test_histEq.bmp"));
    /*
    jbMat mt = jbImgproc::histoPmf(Y);
    jbMat ms = jbImgproc::histoCmf(Y);
    if(!mt.isEmpty()){
        for(int i=0; i < 256; i++ ){
            fprintf(stdout," %3d:%6d,%7d  ", i, (int)mt[i],(int)ms[i]);
            if(i%3==2)
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
    }
    */

    return 0;
}
