#ifndef JZNP_H
#define JZNP_H

#include <stddef.h>

class Jznp
{
    int xsize;
    int ysize;
    int cformat;  // color format
    int channels; // No. of channels
    int bpc;      // bytes per channel
    unsigned char* data;
    size_t data_len; // length of data in byte
    bool empty;
    int unitstep;
public:
    Jznp();
    Jznp(int xsz, int ysz, int cfmt, int ch, int bpc = 8);
    Jznp(Jznp& znp);
    ~Jznp();
    int get_xsize() const { return xsize;}
    int get_ysize() const { return ysize;}
    int get_cformat() const { return cformat;}
    int get_channels() const { return channels;}
    int get_bpc() const { return bpc;}
    unsigned char* get_data_ptr() { return data;}
    int get_unitstep() const { return unitstep; }
    size_t get_data_len() const { return data_len;}
    bool isEmpty() const { return empty; }
    void set_xsize(int xsz) { this->xsize = xsz;}
    void set_ysize(int ysz) { this->ysize = ysz;}
    void set_size (int xsz, int ysz){ this->xsize = xsz; this->ysize = ysz; set_data_len( xsz*ysz);}
    void set_cforamt(int cfmt){ this->cformat = cfmt;}
    void set_channels(int ch) { this->channels = ch;}
    void set_bpc(int bpc) { this->bpc = bpc;}
    void set_data_len(size_t len) { this->data_len = len;}
    void set_unitstep(int step) { this->unitstep = step;}    
};

#endif // JZNP_H
