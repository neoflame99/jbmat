//#include <QCoreApplication>
#include "../../jblib/jbMat.h"
#include "../../jblib/jbmath.h"
#include "../../jblib/qimmat.h"

#include <QImage>
#include <QString>
//#include <memory

int main(int argc, char *argv[])
{
//    QCoreApplication a(argc, argv);
//    return a.exec();
    int row = 7;
    int col = 7;
    int ch = 2;
    jbMat ma(row,col,ch);
    //double* a = ma.getMat();
    double* a = ma.getMat().get();
    int k = 1;
    for(int c=0; c < ch; c++){
        for(int i=0; i < row; i++){
            for( int j=0; j < col; j++)
                a[(i*ma.col+j)*ch+c] = k++;
        }
    }

    k=0;
    ma.printMat();
    jbMat mb(5,5,2);
    jbMat mc(5,5,2);
    //double *b = mb.getMat();
    double* b = mb.getMat().get();
    for(int c=0; c < mb.Nch; c++)
        for(int i=0; i < 5; i++)
            for(int j=0; j < 5; j++){
                b[(i*mb.col+j)*mb.Nch+c] = 100+k++;
                mc(i,j,c) = k;
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

    for(int k=1; k<=3;k++ ){
        QString fname1 = QString("../%1.bmp").arg(k);
        QString fname2 = QString("../%1_1.bmp").arg(k);
        QImage img(fname1);
        jbMat matIm = QimMat::qim2jbmat(img);
        QImage cvim = QimMat::jbmat2qim(matIm);
        cvim.save(fname2);
    }


    return 0;
}
