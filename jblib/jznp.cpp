#include "jznp.h"

Jznp::Jznp()
{
    this->xsize     = 0;
    this->ysize     = 0;
    this->cformat   = 0;
    this->channels  = 0;
    this->bpc       = 0;
    this->data_len  = 0;
    this->data      = 0;
    this->empty     = 0;
    this->unitstep  = 0;
}

Jznp::Jznp(int xsz, int ysz, int cfmt, int ch, int bpc){
    this->xsize     = xsz;
    this->ysize     = ysz;
    this->cformat   = cfmt;
    this->channels  = ch;
    this->bpc       = bpc;

    this->unitstep  = ( bpc <= 1) ? sizeof(char) : sizeof(int);
    this->data_len  = xsz * ysz * unitstep * ch;

    this->data      = new unsigned char [data_len];
    this->empty     = false;
    if(this->data == 0){
        this->empty = true;
        this->data_len = 0;
    }
    else {
        for(size_t i=0; i< data_len; i++)
            this->data[i] = 0;
    }
}

Jznp::Jznp(Jznp& znp){
    this->xsize     = znp.get_xsize();
    this->ysize     = znp.get_ysize();
    this->cformat   = znp.get_cformat();
    this->channels  = znp.get_channels();
    this->bpc       = znp.get_bpc();
    this->data_len  = znp.get_data_len();
    this->data      = new unsigned char [data_len];
    this->unitstep  = znp.get_unitstep();
    this->empty     = false;

    if(this->data == 0){
        this->empty = true;
        this->data_len = 0;
    }else{
        // deep copying
        unsigned char *znp_dat_ptr = znp.get_data_ptr();
        for(size_t i=0; i < data_len; i++)
            this->data[i] = znp_dat_ptr[i];
    }

}

Jznp::~Jznp(){
    if( data != 0)
        delete data;
}
