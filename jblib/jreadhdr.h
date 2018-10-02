#ifndef JREADHDR_H
#define JREADHDR_H

#include <QFile>

class Jreadhdr
{
    size_t ysize, xsize;
    int    rgbe_format;
    bool   isHdrfile;
    unsigned char *rgbe;
    double *dat64f;
    double exposure;
    bool   hasValidData;
    double maxval;
    int    colorformat;
    void   init();
public:
    Jreadhdr();
    Jreadhdr(QFile& file);
    ~Jreadhdr();

    int read_hdr(QFile& file); //return value  -2 : file open error, -1 : memory allocation error, 0 : reading normally
    int conv_rgbe2rgb64f(bool normalizing = true);    //return value  -2 : It doesn't have valid data, -1 : memory allocation error, 0 : conversion ok
    int conv_rgb2xyz();
    int conv_xyz2rgb();
    int conv_rgb2Yxy();
    int conv_Yxy2rgb();
    unsigned char* get_rgbe(){ return rgbe; }
    size_t get_ysize() const { return ysize; }
    size_t get_xsize() const { return xsize; }
    int get_rgbe_format() const { return rgbe_format; } // 0 : rgb, 1: xyz
    bool get_isHdrfile() const {return isHdrfile; }
    double get_exposure() const {return exposure; }
    double* get_dat64f() { return dat64f;}
    double get_maxval() const {return maxval; }
};

#endif // JREADHDR_H
