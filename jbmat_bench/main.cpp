
#include "../jblib/jbMat.h"
#include "../jblib/jbmath.h"
#include "../jblib/qimmat.h"
#include "../jblib/jbimgproc.h"

#include <QImage>
#include <QString>
#include <QDir>
#include <QDebug>
#include <stdlib.h>
#include <sys/time.h>
#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

using namespace jmat;

void fft_test();
void jbMat_test();
void tonemap_test();
void imagerw_test();
void imgproc_test();
void math_test();
void bicubic_test();
int32 main(int32 argc, char *argv[])
{
    //jbMat_test();
    //bicubic_test();
    uint32 row = 7;
    uint32 col = 7;
    uint32 ch = 2;
    uint32 rowcol = row*col;
    Mat ma(DTYP::DOUBLE, row,col,ch,"ma");

    double* a = ma.getDataPtr<double>();
    uint32 k = 1;
    for(uint32 c=0; c < ch; c++){
        for(uint32 i=0; i < row*col; i++){
            if(c == 1)
                a[c*rowcol + i] = k;
            else
                a[c*rowcol + i] = -(int)k;
            a[c*rowcol + i] = k;
            ++k;
        }
    }

    Mat ma_std = ma.std();
    ma_std.printMat("ma's std");

    ma.printMat();
    ma = ma + (double)10;
    ma.setName("ma");
    ma.printMat();

    ma.changeDType(DTYP::CMPLX);
    ma.printMat();

    ma.changeDType(DTYP::INT);
    ma.printMat();

    Mat mb = {cmplx(1.0,2.0),cmplx(3.0,4.0),cmplx(5.0,6.0),cmplx(7.0,8.0)};
    mb.printMat("mb");
    printf("mb data type: %d\n", mb.getDatType());

    Mat mc = Mat::ones(2,2,2, DTYP::CMPLX);
    mc.printMat("mc cmpx ones");
    mc = Mat::ones(2,2,2, DTYP::DOUBLE);
    mc.printMat("mc f64 ones");

    mc = Mat::zeros(2,2,2, DTYP::CMPLX);
    mc.printMat("mc cmpx zero");
    mc = Mat::zeros(2,2,2, DTYP::DOUBLE);
    mc.printMat("mc f64 zero");

    Mat md = { cmplx(1,2),cmplx(3,4),cmplx(5,6),cmplx(7,8),cmplx(9,10),cmplx(11,12),cmplx(13,14),cmplx(15,16),cmplx(17,18),cmplx(19,20),cmplx(21,22),cmplx(23,24)};
    md.reshape(3,2,2);
    md.printMat("md reshape");
    md.transpose();
    md.printMat("md transpose");
    md.reshape(12,1,1);
    Mat md_std = md._std<cmplx>();
    md_std.printMat("md's std1");
    Mat md_std2= md.std();
    md_std.printMat("md's std2");

    cmplx cA = cmplx(2,1);
    cmplx cB = -cA;
    cmplx cC = cA.sqrtc();
    cmplx cD = cB.sqrtc();
    cmplx cE = cmplx(2,-1).sqrtc();
    cmplx cF = cmplx(-2,1).sqrtc();
    printf("%f + %fi, %f +%fi\n", cA.re, cA.im, cB.re, cB.im );
    printf("%f + %fi, %f +%fi\n", cC.re, cC.im, cD.re, cD.im );
    printf("%f + %fi, %f +%fi\n", cE.re, cE.im, cF.re, cF.im );


    Mat me = {1,2,3,4,5,6};
    me.reshape(2,3,1);
    Mat me_ex1 = Mat::repeat(me, 2, 2, 1);
    me.changeDType(DTYP::DOUBLE);
    Mat me_ex2 = Mat::repeat(me, 2, 2, 2);
    Mat me2 = {cmplx(1,2),cmplx(3,4),cmplx(5,6),cmplx(7,8)};
    Mat me2_ex1 = Mat::repeat(me2, 2, 1, 3);
    Mat me2_ex2 = Mat::repeat(me2, 1, 2, 3);
    me.printMat("me") ;
    me_ex1.printMat("me_ex1");
    me_ex2.printMat("me_ex2");
    me2.printMat("me2") ;
    me2_ex1.printMat("me2_ex1");
    me2_ex2.printMat("me2_ex2");

    return 0;
}
void bicubic_test(){
    int rc=7;
    Mat ma(DTYP::INT,rc,rc,1);
    for(int32 i=0; i < rc*rc; ++i)
       ma.at<int>(i) = rand()%100;
    ma.printMat();
    Mat mb = imgproc::bicubicIntp(ma,2);
    mb.printMat();

    QString fname1 = QString("../jbmat_bench/test.jpg");
    QString fname2 = QString("../jbmat_bench/test_scale2.bmp");
    QImage img(fname1);
    Mat matIm = qimmat::qim2mat(img);
    Mat imScale = imgproc::bicubicIntp(matIm,2);
    QImage cvim = qimmat::mat2qim(imScale);
    cvim.save(fname2);
}

void jbMat_test(){

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
    maa.printMat("maa");
    //maa.plusMat(ma);
    Mat ma_ch1 = ma.copyChannelN(1);
    ma_ch1.printMat("ma_ch1");
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

    Mat maf_submat = maf.copySubMat(0,3,1,3);
    maf_submat.printMat("maf submat [0:3, 1:3]");
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
    Mat mak = mai._mean<double>();
    Mat mal = mai._std<double>();
    mak.printMat("mak : mean of mai");
    mal.printMat("mal : std of mai");

    man.printMat();
    man._max<double>().printMat("man max");
    man._min<double>().printMat("man min");



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


    Mat mdps3({1,2,3,4,5,6});
    Mat mdps4({7,8,9,10,11,12});
    mdps3 *= 1000;
    mdps3.printMat("mdps3 transpose input");
    mdps3.transpose();
    mdps3.printMat("mdps3 transpose output");
    mdps3.transpose();
    mdps3.printMat("mdps3 reverse transpose output");

    Mat tt;
    {
        tt= Mat(DTYP::INT,1,1,1);
        tt.printMat();
    }
    /* ------- Mat repeat ------ */
    Mat mdps3rp = Mat::repeat(mdps3,2,3,2);
    mdps3rp.printMat("mdps3 repeat (2,3,2)");
    Mat mdps3rp2 = Mat::repeat(mdps3,1,1,3);
    mdps3rp2.printMat("mdps3 repeat 2nd (1,1,3)");

    _complex ca(100,5);
    _complex cb(4,-2);
    _complex da = ca * cb;
    _complex db = ca / cb;
    _complex dc = ca + cb;
    _complex dd = ca - cb;
    printf("ca*cb=%f+j%f, ca/cb=%f+j%f, ca+cb=%f+j%f, ca-cb=%f+j%f\n",da.re, da.im, db.re, db.im, dc.re, dc.im, dd.re, dd.im);

    /* ---- sliceCopyMat  ----*/
    printf("slice copy\n");
    Mat mB = Mat::zeros(8,8,2,ma.getDatType());
    matRect Ra(0,0,3,4);
    matRect Rb(2,2,5,6);

    mB.printMat();
    Mat::sliceCopyMat(ma,Ra,mB,Rb);
    mB.printMat();
    ma.printMat();

    /* --- copy_padding -----*/
    Mat mC = imgproc::copy_padding(ma,2);
    mC.printMat();
}
void fft_test(){

    int32 len = 16;
    _complex *dat1 = new _complex[len];
    _complex *dat2 = new _complex[len];

    for(int32 i=0; i < len; ++i){
        dat1[i] = _complex(i+1, len-i);
        dat2[i] = dat1[i];
        printf("%3d : %f %+fj\n",i, dat1[i].re, dat1[i].im);
    }
    imgproc::fft(dat1, len);  // fft
    imgproc::ifft(dat2, len);  // ifft

    FILE *fft_fd;
    fft_fd = fopen("../fft_result.txt","w");

    fprintf(fft_fd, "FFT: DIT \n");
    for(int32 i=0; i < len; ++i){
        fprintf(fft_fd, "%3d : % 8.6f %+8.6fj \n",i, dat1[i].re, dat1[i].im);
    }
    fprintf(fft_fd, "IFFT: DIT \n");
    for(int32 i=0; i < len; ++i){
        fprintf(fft_fd, "%3d : % 8.6f %+8.6fj \n",i, dat2[i].re, dat2[i].im);
    }
    fclose(fft_fd);

    double max_d1, max_d2, min_d1, min_d2;
    max_d1 = abs(dat1[0].re) > abs(dat1[0].im) ? abs(dat1[0].re) : abs(dat1[0].im);
    min_d1 = abs(dat1[0].re) < abs(dat1[0].im) ? abs(dat1[0].im) : abs(dat1[0].re);
    max_d2 = abs(dat2[0].re) > abs(dat2[0].im) ? abs(dat2[0].re) : abs(dat2[0].im);
    min_d2 = abs(dat2[0].re) < abs(dat2[0].im) ? abs(dat2[0].im) : abs(dat2[0].re);
    for(int32 i=1; i < len; ++i){
        if( max_d1 < abs(dat1[i].re)) max_d1 = abs(dat1[i].re);
        if( max_d1 < abs(dat1[i].im)) max_d1 = abs(dat1[i].im);
        if( max_d2 < abs(dat2[i].re)) max_d2 = abs(dat2[i].re);
        if( max_d2 < abs(dat2[i].im)) max_d2 = abs(dat2[i].im);

        if( min_d1 > abs(dat1[i].re)) min_d1 = abs(dat1[i].re);
        if( min_d1 > abs(dat1[i].im)) min_d1 = abs(dat1[i].im);
        if( min_d2 > abs(dat2[i].re)) min_d2 = abs(dat2[i].re);
        if( min_d2 > abs(dat2[i].im)) min_d2 = abs(dat2[i].im);
    }
    printf("FFT:\n");

    double rat;
    rat = 1.0;
    if (max_d1 > 10 ){
        rat = pow(10,floor(log10(max_d1)));
        printf("\t%.3E \n",rat);
    }
    for(int32 i=0; i < len; ++i){
        printf("%3d : % 8.4f %+8.4fj \n",i, dat1[i].re/rat, dat1[i].im/rat);
    }

    printf("IFFT:\n");
    rat = 1.0;
    if(max_d2 > 10){
        rat = pow(10,floor(log10(max_d2)));
        printf("\t%.3E \n",rat);
    }
    for(int32 i=0; i < len; ++i)
        printf("%3d : % 8.4f %+8.4fj \n",i, dat2[i].re/rat, dat2[i].im/rat);


    printf("data\n");
    for(int32 i=0; i < len; ++i){
        dat1[i] = _complex(i+1, len-i);
        dat2[i] = dat1[i];
    }
    for(int32 i=0; i < len; ++i)
        printf("%3d : % 8.4f %+8.4fj \n",i, dat2[i].re, dat2[i].im);

    imgproc::fft2d(dat1, 4, 3);

    max_d1 = abs(dat1[0].re) > abs(dat1[0].im) ? abs(dat1[0].re) : abs(dat1[0].im);
    min_d1 = abs(dat1[0].re) < abs(dat1[0].im) ? abs(dat1[0].im) : abs(dat1[0].re);
    printf("2D FFT:\n");
    rat = 1.0;
    if (max_d1 > 10 ){
        rat = pow(10,floor(log10(max_d1)));
        printf("\t%.3E \n",rat);
    }
    for(int32 i=0; i < len; ++i){
        printf("%3d : % 8.4f %+8.4fj \n",i, dat1[i].re/rat, dat1[i].im/rat);
    }

    imgproc::ifft2d(dat1, 4, 3);

    for(int32 i=1; i < len; ++i){
        if( max_d1 < abs(dat1[i].re)) max_d1 = abs(dat1[i].re);
        if( max_d1 < abs(dat1[i].im)) max_d1 = abs(dat1[i].im);
        if( max_d2 < abs(dat2[i].re)) max_d2 = abs(dat2[i].re);
        if( max_d2 < abs(dat2[i].im)) max_d2 = abs(dat2[i].im);
    }
    max_d2 = abs(dat2[0].re) > abs(dat2[0].im) ? abs(dat2[0].re) : abs(dat2[0].im);
    min_d2 = abs(dat2[0].re) < abs(dat2[0].im) ? abs(dat2[0].im) : abs(dat2[0].re);
    for(int32 i=1; i < len; ++i){
        if( min_d1 > abs(dat1[i].re)) min_d1 = abs(dat1[i].re);
        if( min_d1 > abs(dat1[i].im)) min_d1 = abs(dat1[i].im);
        if( min_d2 > abs(dat2[i].re)) min_d2 = abs(dat2[i].re);
        if( min_d2 > abs(dat2[i].im)) min_d2 = abs(dat2[i].im);
    }

    printf("2D IFFT:\n");
    rat = 1.0;
    if(max_d2 > 10){
        rat = pow(10,floor(log10(max_d2)));
        printf("\t%.3E \n",rat);
    }
    for(int32 i=0; i < len; ++i)
        printf("%3d : % 8.4f %+8.4fj \n",i, dat2[i].re/rat, dat2[i].im/rat);

    delete [] dat1;
    delete [] dat2;

    printf("FFT DIF4\n");
    len = 16;
    _complex *dat3 = new _complex[len];

    printf(" DATA:\n");
    for(int32 i=0; i < len; ++i){
        dat3[i] = _complex(i, len-i);
        printf("%3d : %f %+fj\n",i, dat3[i].re, dat3[i].im);
    }
    imgproc::fft_dif4(dat3, len, false); // fft

    printf(" FFT:\n");
    for(int32 i=0; i < len; ++i)
        printf("%3d : % 8.4f %+8.4fj \n",i, dat3[i].re, dat3[i].im);

    imgproc::fft_dif4(dat3, len, true); // ifft
    printf(" IFFT:\n");
    for(int32 i=0; i < len; ++i)
        printf("%3d : % 8.4f %+8.4fj \n",i, dat3[i].re, dat3[i].im);

    delete [] dat3;

    std::vector<int32> ff;
    imgproc::factorizeN(187, ff);
    printf(" factorizing:\n");
    for(int32 i=0; i < ff.size(); ++i)
        printf("%3d \n",ff.at(i));


    len = 64;
    struct timeval  tv1, tv2;

    _complex *dat4 = new _complex[len];
    _complex *dat5 = new _complex[len];

    printf(" DATA:\n");
    for(int32 i=0; i < len; ++i){
        dat4[i] = _complex(i+1, 0); //len-i);
        dat5[i] = _complex(i+1, 0); //len-i);
        printf("%3d : %f %+fj\n",i, dat4[i].re, dat4[i].im);
    }
    ff.clear();

    gettimeofday(&tv1, NULL);
    imgproc::factorizeN(len, ff);
    imgproc::fft_compositN(dat4, len, ff, false); // fft
    gettimeofday(&tv2, NULL);

    printf ("Total time = %f seconds\n",
             (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
             (double) (tv2.tv_sec - tv1.tv_sec));

    gettimeofday(&tv1, NULL);
    imgproc::fft_czt(dat5,len,false);
    gettimeofday(&tv2, NULL);

    printf ("Total time = %f seconds\n",
             (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
             (double) (tv2.tv_sec - tv1.tv_sec));

    gettimeofday(&tv1, NULL);
    imgproc::fft_dif4(dat5,len,false);
    gettimeofday(&tv2, NULL);

    printf ("Total time = %f seconds\n",
             (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
             (double) (tv2.tv_sec - tv1.tv_sec));
    gettimeofday(&tv1, NULL);
    imgproc::fft_dit2(dat5,len,false);
    gettimeofday(&tv2, NULL);

    printf ("Total time = %f seconds\n",
             (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
             (double) (tv2.tv_sec - tv1.tv_sec));


    imgproc::fft_compositN(dat4, len, ff, true); // fft
    for(int32 i=0; i < len; ++i)
        printf("%3d : % 8.4f %+8.4fj \n",i, dat4[i].re, dat4[i].im);

    delete [] dat4;
    delete [] dat5;
}

void tonemap_test(){
    /* ----- */

    /*--- gauss mask / box mask gen ---*/

    Mat gau = imgproc::gaussMaskGen(1,4);
    Mat box = imgproc::boxMaskGen(5);

    gau.printMat("gaussian Mat");
    box.printMat("box Mat");

    /*------------------*/

    /*----- tonemapping ------ */
    /*
    imgproc::gamma(hdrimg, 1.0);
    Mat Yhdrimg    = imgproc::rgb2gray(hdrimg);
    Mat Yhdrimg_tm = imgproc::nakaSigTonemap(Yhdrimg,gau,0.2);
    Mat Yhdrimg_tm3c= Mat::repeat(Yhdrimg_tm, 1, 1,3);
    Mat Yhdrimg_3c  = Mat::repeat(Yhdrimg, 1, 1, 3);

    Mat hdrimg_tm = imgproc::gamma( hdrimg * Yhdrimg_tm3c / Yhdrimg_3c, 0.45);
    hdrimg_tm = hdrimg_tm / hdrimg_tm.max().max().at<double>(0)*255.0;
    QImage hdrtm_bmp = qimmat::mat2qim(hdrimg_tm);
    hdrtm_bmp.save("../jbmat_bench/hdr_tm_0.2.bmp");
    */
    /*-------------------------*/
}

void imagerw_test(){
    /* hdrfile reading and writing to bmp */
    QString hdrfile_path = QString("/Users/neoflame99/Workspace/Qt5/readhdrfile/readhdrfile/memorial.hdr");
    Mat hdrimg = qimmat::read_hdr(hdrfile_path);
    Mat hdrimg_sub = hdrimg.copySubMat(0,1,0,1);
    hdrimg_sub.printMat("hdrimg_sub");
    QImage hdr_bmp = qimmat::mat2qim(hdrimg);
    hdr_bmp.save("../jbmat_bench/hdr_bmp.bmp");

}
void imgproc_test(){
   /*
    qDebug() << QDir::currentPath();
    QString fname1 = QString("../jbmat_bench/test.jpg");
    QString fname2 = QString("../jbmat_bench/test_filt.bmp");
    QImage img(fname1);
    Mat matIm = qimmat::qim2mat(img);
    Mat filt  =  Mat::ones(5,5,matIm.getChannel())/25;
    Mat FiltIm = conv2d(matIm, filt,"zero","same");
    Mat Y  = imgproc::rgb2gray(matIm);
    Mat yccIm = imgproc::rgb2ycc(matIm);
    Mat ms = imgproc::clip_HistoCmf(Y, 1000, 256, 1);
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
}
void math_test(){


    Mat mk(DTYP::DOUBLE,3,4,2,"mk");
    Mat ml(DTYP::DOUBLE,4,3,2,"ml");
    Mat mn(DTYP::DOUBLE,4,4,1,"mn");
    for(int32 i=0; i < 24; i++){
        mk.at<double>(i) = rand() % 50;
        ml.at<double>(i) = rand() % 80;
        if(i < 16)
            mn.at<double>(i) = rand() % 100;
    }
    Mat mm = mul(mk,ml);
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

    //-- conv
    int32 k=0;
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

    uint32 row = 7;
    uint32 col = 7;
    uint32 ch = 2;
    uint32 rowcol = row*col;
    Mat ma(DTYP::DOUBLE, row,col,ch,"ma");

    ma.printMat("ma");
    mb.printMat("mb");
    fprintf(stdout,"1) full\n");
    Mat convO = conv2d(ma,mb,"zero","full");
    convO.printMat("conv0: full");
    fprintf(stdout,"2) same\n");
    convO = conv2d(ma,mb,"zero","same");
    convO.printMat("conv0: same");
}
