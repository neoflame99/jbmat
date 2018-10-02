#include "tmo.hpp"

const double etha = 2.3e-5;

 int tmLog(Jreadhdr& hdrim){
    size_t ysize, xsize, len;
    ysize = hdrim.get_ysize();
    xsize = hdrim.get_xsize();
    len   = ysize * xsize * 3;

    hdrim.conv_rgb2Yxy();
    double* dat64f = hdrim.get_dat64f();
    if( dat64f == 0) return -1;

    double Lwmax = -1;
    double Lw;
    double Ylogsum, Ylogavg;
    Ylogsum = 0;
    for(size_t i=0 ; i < len; i +=3)
        Ylogsum += log(etha + dat64f[i]);

    Ylogavg = Ylogsum / (ysize*xsize);

    double Lwa = exp(Ylogavg);

    // applying world adaptation luminace
    for(size_t i=0; i < len; i +=3){
        dat64f[i] = dat64f[i]/Lwa;
        if(dat64f[i] > Lwmax)
            Lwmax = dat64f[i];
    }

    for(size_t i=0; i < len; i+=3){
        Lw = dat64f[i];
        dat64f[i] = log(Lw+1)/log(Lwmax+1);
    }
    hdrim.conv_Yxy2rgb();

    // gamma bt709
    for(size_t i=0; i < len; i++)
        dat64f[i] = (dat64f[i] <= 0.018) ? 4.5*dat64f[i] : pow(dat64f[i],0.45)*1.099-0.099;

    return 0;
}


 int tmDrago(Jreadhdr& hdrim, double bias_b){
    size_t ysize, xsize;
    ysize = hdrim.get_ysize();
    xsize = hdrim.get_xsize();
    double* dat64f = hdrim.get_dat64f();
    size_t len = ysize * xsize *3;
    hdrim.conv_rgb2Yxy();

    // find log average
    double Ylogsum=0, Ylogavg=0;
    double Lwmax = 0;
    double Ld, Lw ;
    double Ldmax = 1;
    double bias_exponent = log(bias_b) / log(0.5);

    for(size_t i=0 ; i < len; i +=3)
        Ylogsum += log(etha + dat64f[i]);

    Ylogavg = Ylogsum / (ysize*xsize);

    double Lwa = exp(Ylogavg);

    // applying world adaptation luminace
    for(size_t i=0; i < len; i +=3){
        dat64f[i] = dat64f[i]/Lwa;
        if(dat64f[i] > Lwmax)
            Lwmax = dat64f[i];
    }

    // Adaptive Tonemapping
    double denom   = log10(Lwmax+1);
    double Lw_norm = 0.0;
    double bias;
    for(size_t i=0; i < len; i+=3){
        Lw        = dat64f[i];
        Lw_norm   = Lw / Lwmax;
        bias      = pow( Lw_norm, bias_exponent);
        if(bias==0.0)
            Lw_norm = Lw_norm;
        Ld        = Ldmax * log(Lw+1)/(denom * log(2+bias*8));
        dat64f[i] = Ld;
    }

    hdrim.conv_Yxy2rgb();
    // gamma bt709
    for(size_t i=0; i < len; i++)
        dat64f[i] = (dat64f[i] <= 0.018) ? 4.5*dat64f[i] : pow(dat64f[i],0.45)*1.099-0.099;

    return 0;
}
