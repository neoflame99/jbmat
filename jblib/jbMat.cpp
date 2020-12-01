#include <stdio.h>
#include "jbMat.h"
#include <iostream>
#include <utility>

namespace jmat {


int32 Mat::instant_count = 0;

//-- shallow copy version & using shared_ptr
inline void Mat::alloc(const uint32 len){

    mA = len==0 ? nullptr : shr_ptr (new uchar[len], std::default_delete<uchar[]>());
    /*
    try{
        mA = (len==0)? nullptr : new double[len];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"memory allocation error: %s\n",ex.what());
        mA = nullptr;
    }
    */
}
void Mat::init(uint32 r, uint32 c, uint32 ch, DTYP dt, bool do_alloc){
    row = r;
    col = c;
    Nch = ch;
    stepCol   = Nch;
    stepRow   = col*stepCol;
    length    = row*stepRow;
    lenRowCol = row*col;

    datT = dt;
    switch(datT){
    case DTYP::UCHAR  : byteStep = 1; break;
    case DTYP::INT    : byteStep = 4; break;
    case DTYP::FLOAT  : byteStep = 4; break;
    case DTYP::DOUBLE : byteStep = 8; break;
    case DTYP::CMPLX  : byteStep = sizeof(cmplx); break;
    }
    byteLen = length*byteStep;
    if(do_alloc)
        alloc(byteLen);
    sync_data_ptr();
}
void Mat::initName(){
    obj_name =std::string( "Mat_") + std::to_string (instant_count);
    instant_count++;
}

Mat::Mat(DTYP dt):mA(nullptr){
    init(0,0,0, dt, true);
    initName();
}
Mat::Mat(DTYP dt, uint32 rc ):mA(nullptr){
    init(rc, rc, 1, dt, true);
    initName();
}
Mat::Mat(DTYP dt, uint32 r, uint32 c, uint32 ch):mA(nullptr){
    init(r, c, ch, dt, true);
    initName();
}
Mat::Mat(shr_ptr ma, DTYP dt, uint32 r, uint32 c, uint32 ch):mA(ma){
    init(r, c, ch, dt, false);
    initName();
}
Mat::Mat(DTYP dt, uint32 r, uint32 c, uint32 ch, std::string name):mA(nullptr){
    init(r, c, ch, dt, true);
    obj_name = name;
}
Mat::Mat(const Mat& mat){
    row = mat.row;
    col = mat.col;
    Nch = mat.Nch;
    length    = mat.length;
    stepCol   = mat.stepCol;
    stepRow   = mat.stepRow;
    lenRowCol = mat.lenRowCol;
    mA        = mat.mA;
    datT      = mat.datT;
    byteStep  = mat.byteStep;
    byteLen   = mat.byteLen;
    sync_data_ptr();
    initName();

    /* -- deep copying
    alloc(length);
    double *pt_matdat = mat.getMat().get();
    double *pt_thisdat = mA.get();
    std::copy(pt_matdat,pt_matdat+mat.length,pt_thisdat);
    */
#ifdef _DEBUG_MODE_
    fprintf(stdout,"copy constructor\n");
#endif
}

Mat::Mat( std::initializer_list<double> list ){
    //-- Making a column vector
    init(list.size(),1,1,DTYP::DOUBLE, true);

    if(elptr.f64_ptr !=nullptr){
        memcpy(elptr.f64_ptr, list.begin(), list.size()*sizeof(double));
    }
}
Mat::Mat( std::initializer_list<float> list ){
    //-- Making a column vector
    init(list.size(),1,1,DTYP::FLOAT, true);

    if(elptr.f32_ptr!=nullptr){
        memcpy(elptr.f32_ptr, list.begin(), list.size()*sizeof(float));
    }
}
Mat::Mat( std::initializer_list<int32> list ){
    //-- Making a column vector
    init(list.size(),1,1,DTYP::INT, true);

    if(elptr.int_ptr!=nullptr){
        memcpy(elptr.int_ptr, list.begin(), list.size()*sizeof(int32));
    }
}
Mat::Mat( std::initializer_list<uchar> list ){
    //-- Making a column vector
    init(list.size(),1,1,DTYP::UCHAR, true);

    if(elptr.uch_ptr!=nullptr){
        memcpy(elptr.uch_ptr, list.begin(), list.size()*sizeof(uchar));
    }
}
Mat::Mat( std::initializer_list<cmplx> list ){
    //-- Making a column vector
    init(list.size(),1,1,DTYP::CMPLX, true);

    if(elptr.cmx_ptr!=nullptr){
        //memcpy(elptr.cmx_ptr, list.begin(), list.size()*sizeof(cmplx));
        std::vector<cmplx> v(list);
        uint32 k = 0;
        for(std::vector<cmplx>::iterator itr = v.begin(); itr != v.end(); itr++){
            elptr.cmx_ptr[k++] = *itr;
        }
    }
}

Mat::~Mat(){}

void Mat::setRowCol(uint32 r, uint32 c, uint32 ch){
    uint32 lenrc    = r*c;
    uint32 step_col = ch;
    uint32 step_row = c*step_col;
    uint32 len      = r*step_row;

    if(length != len){
        mA.reset();
        alloc(len);
    }
    row = r;
    col = c;
    Nch = ch;
    lenRowCol = lenrc;
    stepCol   = step_col;
    stepRow   = step_row;
    length    = len;
    sync_data_ptr();
}

/*
//-- overloading operators : it calls copy constructor
Mat& Mat::operator=( Mat other){
//The parameter to the ‘operator=()’ is passed by value which calls copy constructor
//to create an object local to the ‘operator=()’.
//Than the value of the temp object is swapped with ‘*this’ object
//
    std::swap(row,other.row);
    std::swap(col,other.col);
    std::swap(Nch,other.Nch);

    std::swap(length,other.length);
    std::swap(lenRowCol,other.lenRowCol);
    std::swap(mA,other.mA); // swapping mA pointer

    return *this;
}
*/

// rvalue reference
Mat& Mat::operator=(Mat&& other){

#ifdef _DEBUG_MODE_
    fprintf(stdout,"rvalue referance Assign operator\n");
#endif
    std::swap(row,other.row);
    std::swap(col,other.col);
    std::swap(Nch,other.Nch);
    std::swap(length,other.length);
    std::swap(stepCol,other.stepCol);
    std::swap(stepRow,other.stepRow);
    std::swap(lenRowCol,other.lenRowCol);
    std::swap(byteStep,other.byteStep);
    std::swap(byteLen,other.byteLen);
    std::swap(datT,other.datT);
    std::swap(mA, other.mA);
    sync_data_ptr();

    return *this;
}

Mat& Mat::operator+=(const Mat& other){
    return plusMat(other);
}
Mat& Mat::operator+=(const double scalar){
    return plusScalar(scalar);
}

Mat& Mat::operator-=(const Mat& other){
    return minusMat(other);
}
Mat& Mat::operator-=(const double scalar){
    return minusScalar(scalar);
}

Mat& Mat::operator*=(const Mat& other){
    return mulMat(other);
}
Mat& Mat::operator*=(const double scalar){
    return mulScalar(scalar);
}

Mat& Mat::operator/=(const Mat& other){
    return divMat(other);
}
Mat& Mat::operator/=(const double scalar){
    return divByScalar(scalar);
}

// unary minus or minus sign op.
Mat Mat::operator-() const{
    if(this->isEmpty()){
        fprintf(stdout,"Mat::operator- (unary): Stop adding two operands because the both operands are empty\n");
        return Mat();
    }

    Mat A = this->copy();
    uint32 k = 0;
    switch(datT){
    case DTYP::CMPLX : for(; k < length; ++k) { A.elptr.cmx_ptr[k] = -A.elptr.cmx_ptr[k];} break;
    case DTYP::DOUBLE: for(; k < length; ++k) { A.elptr.f64_ptr[k] = -A.elptr.f64_ptr[k];} break;
    case DTYP::FLOAT : for(; k < length; ++k) { A.elptr.f32_ptr[k] = -A.elptr.f32_ptr[k];} break;
    case DTYP::INT   : for(; k < length; ++k) { A.elptr.int_ptr[k] = -A.elptr.int_ptr[k];} break;
    case DTYP::UCHAR : for(; k < length; ++k) { A.elptr.uch_ptr[k] = -A.elptr.uch_ptr[k];} break;
    }
    return A;
}

Mat Mat::operator+(const Mat& other) const{
    if( row != other.getRow() || col != other.getCol()|| Nch != other.getChannel() ){
        fprintf(stdout,"Mat::operator+ : Stop adding two operands because the both operands are not same size\n");
        return Mat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Mat::operator+ : Stop adding two operands because the both operands are empty\n");
        return Mat();
    }

    Mat sum = this->copy();
    sum.plusMat(other);

    return sum;
}
Mat Mat::operator-(const Mat& other) const{
    if( row != other.getRow() || col != other.getCol() || Nch != other.getChannel() ){
        fprintf(stdout,"Mat::operator- : Stop adding two operands because the both operands are not same size\n");
        return Mat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Mat::operator- : Stop adding two operands because the both operands are empty\n");
        return Mat();
    }

    Mat subtract = this->copy();
    subtract.minusMat(other);

    return subtract;
}
Mat Mat::operator*(const Mat& other) const{
    if( row != other.getRow() || col != other.getCol() || Nch != other.getChannel() ){
        fprintf(stdout,"Mat::operator* : Stop adding two operands because the both operands are not same size\n");
        return Mat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Mat::operator* : Stop adding two operands because the both operands are empty\n");
        return Mat();
    }
    Mat product = this->copy();
    product.mulMat(other);

    return product;
}
Mat Mat::operator/(const Mat& other) const{
    if( row != other.getRow() || col != other.getCol() || Nch != other.getChannel() ){
        fprintf(stdout,"Mat::operator/ : Stop adding two operands because the both operands are not same size\n");
        return Mat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Mat::operator/ : Stop adding two operands because the both operands are empty\n");
        return Mat();
    }
    Mat div = this->copy();
    div.divMat(other);

    return div;
}

Mat operator+(const Mat& A, const double scalar){
    if(A.isEmpty()){
        fprintf(stderr,"Mat::operator+ with scalar : Mat A is empty \n");
        return Mat();
    }
    Mat sum = A.copy();
    sum.plusScalar(scalar);

    return sum;
}
Mat operator+(const double scalar , const Mat& A){
    if(A.isEmpty()){
        fprintf(stderr,"Mat::operator+ with scalar : Mat A is empty \n");
        return Mat();
    }
    Mat sum = A.copy();
    sum.plusScalar(scalar);

    return sum;
}

Mat operator-(const Mat& A, const double scalar) {
    if(A.isEmpty()){
        fprintf(stderr,"Mat::operator- with scalar : this mat is empty \n");
        return Mat();
    }
    Mat subtract = A.copy();
    subtract.minusScalar(scalar);

    return subtract;
}
Mat operator-(const double scalar,const Mat& A) {
    if(A.isEmpty()){
        fprintf(stderr,"Mat::operator- with scalar : this mat is empty \n");
        return Mat();
    }
    Mat subtract = A.copy();
    subtract.minusScalar(scalar);

    return subtract;
}

Mat operator*(const Mat& A, const double scalar){
    if(A.isEmpty()){
        fprintf(stderr,"Mat::operator* with scalar : this mat is empty \n");
        return Mat();
    }
    Mat product = A.copy();
    product.mulScalar(scalar);

    return product;
}
Mat operator*(const double scalar, const Mat& A){
    if(A.isEmpty()){
        fprintf(stderr,"Mat::operator* with scalar : this mat is empty \n");
        return Mat();
    }
    Mat product = A.copy();
    product.mulScalar(scalar);

    return product;
}

Mat operator/(const Mat& lhs, const double scalar){
    if(lhs.isEmpty()){
        fprintf(stderr,"Mat::operator/ with scalar : this mat is empty \n");
        return Mat();
    }
    Mat div = lhs.copy();
    div.divByScalar(scalar);

    return div;
}
Mat operator/(const double scalar, const Mat& rhs){
    if(rhs.isEmpty()){
        fprintf(stderr,"Mat::operator/ with scalar : this mat is empty \n");
        return Mat();
    }
    Mat div = rhs.copy();
    div.divScalar(scalar);

    return div;
}

int32 Mat::reshape(uint32 r, uint32 c, uint32 ch){
    uint32 rc   = r*c;
    uint32 tlen = rc*ch;
    if( tlen != length){
        fprintf(stderr,"Mat::reshape : reshape argument is not correct!\n");
        return -1;
    }
    row = r;
    col = c;
    Nch = ch;
    stepCol = Nch;
    stepRow = col * stepCol;
    lenRowCol = rc;

    return 0;
}

void Mat::changeDType(const DTYP dt){

    assert(static_cast<int>(dt) < 5 && static_cast<int>(datT) < 5 );

    if( dt == datT ) return ;

    elemptr t_elptr;
    t_elptr.uch_ptr = nullptr;
    try{
        switch(dt){
        case DTYP::UCHAR  : t_elptr.uch_ptr = new uchar[(U64)length * sizeof(uchar )]; break;
        case DTYP::INT    : t_elptr.uch_ptr = new uchar[(U64)length * sizeof(int32 )]; break;
        case DTYP::FLOAT  : t_elptr.uch_ptr = new uchar[(U64)length * sizeof(float )]; break;
        case DTYP::DOUBLE : t_elptr.uch_ptr = new uchar[(U64)length * sizeof(double)]; break;
        case DTYP::CMPLX  : t_elptr.uch_ptr = new uchar[(U64)length * sizeof(cmplx )];
        }
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"Mat::changeDType : memory allocation error: %s\n",ex.what());
        return;
    }

    uint32  k = 0;
    if(dt == DTYP::DOUBLE){
        switch (datT) {
        case DTYP::FLOAT : for( ; k < length ; ++k ){ t_elptr.f64_ptr[k] = (double) elptr.f32_ptr[k]; } break;
        case DTYP::INT   : for( ; k < length ; ++k ){ t_elptr.f64_ptr[k] = (double) elptr.int_ptr[k]; } break;
        case DTYP::UCHAR : for( ; k < length ; ++k ){ t_elptr.f64_ptr[k] = (double) elptr.uch_ptr[k]; } break;
        case DTYP::CMPLX : for( ; k < length ; ++k ){ t_elptr.f64_ptr[k] =          elptr.cmx_ptr[k].re; } break;
        case DTYP::DOUBLE: ;
        }
        byteStep = sizeof(double);
    }else if( dt== DTYP::FLOAT){
        switch (datT) {
        case DTYP::DOUBLE: for( ; k < length ; ++k ){ t_elptr.f32_ptr[k] = (float) elptr.f64_ptr[k]; } break;
        case DTYP::INT   : for( ; k < length ; ++k ){ t_elptr.f32_ptr[k] = (float) elptr.int_ptr[k]; } break;
        case DTYP::UCHAR : for( ; k < length ; ++k ){ t_elptr.f32_ptr[k] = (float) elptr.uch_ptr[k]; } break;
        case DTYP::CMPLX : for( ; k < length ; ++k ){ t_elptr.f32_ptr[k] = (float) elptr.cmx_ptr[k].re; } break;
        case DTYP::FLOAT : ;
        }
        byteStep = sizeof(float);
    }else if( dt== DTYP::INT){
        switch (datT) {
        case DTYP::DOUBLE: for( ; k < length ; ++k ){ t_elptr.int_ptr[k] = (int32) elptr.f64_ptr[k]; } break;
        case DTYP::FLOAT : for( ; k < length ; ++k ){ t_elptr.int_ptr[k] = (int32) elptr.f32_ptr[k]; } break;
        case DTYP::UCHAR : for( ; k < length ; ++k ){ t_elptr.int_ptr[k] = (int32) elptr.uch_ptr[k]; } break;
        case DTYP::CMPLX : for( ; k < length ; ++k ){ t_elptr.int_ptr[k] = (int32) elptr.cmx_ptr[k].re; } break;
        case DTYP::INT   : ;
        }
        byteStep = sizeof(int);
    }else if( dt== DTYP::UCHAR){
        switch (datT) {
        case DTYP::DOUBLE: for( ; k < length ; ++k ){ t_elptr.uch_ptr[k] = (uchar) elptr.f64_ptr[k]; } break;
        case DTYP::FLOAT : for( ; k < length ; ++k ){ t_elptr.uch_ptr[k] = (uchar) elptr.f32_ptr[k]; } break;
        case DTYP::INT   : for( ; k < length ; ++k ){ t_elptr.uch_ptr[k] = (uchar) elptr.int_ptr[k]; } break;
        case DTYP::CMPLX : for( ; k < length ; ++k ){ t_elptr.uch_ptr[k] = (uchar) elptr.cmx_ptr[k].re; } break;
        case DTYP::UCHAR : ;
        }
        byteStep = sizeof(uchar);
    }else if( dt== DTYP::CMPLX){
        switch (datT) {
        case DTYP::DOUBLE: for( ; k < length ; ++k ){ t_elptr.cmx_ptr[k] = cmplx( elptr.f64_ptr[k], 0.0); } break;
        case DTYP::FLOAT : for( ; k < length ; ++k ){ t_elptr.cmx_ptr[k] = cmplx( elptr.f32_ptr[k], 0.0); } break;
        case DTYP::INT   : for( ; k < length ; ++k ){ t_elptr.cmx_ptr[k] = cmplx( elptr.int_ptr[k], 0.0); } break;
        case DTYP::UCHAR : for( ; k < length ; ++k ){ t_elptr.cmx_ptr[k] = cmplx( elptr.uch_ptr[k], 0.0); } break;
        case DTYP::CMPLX : ;
        }
        byteStep = sizeof(cmplx);
    }
    mA.reset( t_elptr.uch_ptr, std::default_delete<uchar[]>());
    sync_data_ptr();

    datT    = dt;
    byteLen = length * byteStep;

#ifdef _DEBUG
    switch(datT){
    case DTYP::DOUBLE: fprintf(stdout,"datT is changed to double\n"); break;
    case DTYP::FLOAT : fprintf(stdout,"datT is changed to float \n"); break;
    case DTYP::INT   : fprintf(stdout,"datT is changed to int32 \n"); break;
    case DTYP::UCHAR : fprintf(stdout,"datT is changed to uchar \n"); break;
    case DTYP::CMPLX : fprintf(stdout,"datT is changed to uchar \n"); break;
    }
#endif
}

void Mat::transpose(){
    if(isEmpty()){
        fprintf(stderr,"Mat::transpose: This Mat is empty\n");
        return ;
    }

    uchar *tmA;
    try{
        tmA = new uchar[static_cast<unsigned long>(byteLen)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"Mat::transpose Error: %s\n",ex.what());
        return;
    }

    uint32 i, j, lhs_idx, rhs_idx;
    uint32 rr, rf, lf, lr;
    uint32 byte_step_col= stepCol * byteStep;
    uint32 rr_step      = col*byte_step_col;  // rough step of right-hand index
    uint32 lr_step      = row*byte_step_col;  // rough step of left-hand index
    for(i=0, rr=0, lf=0; i < row; i++, rr+=rr_step, lf+=byte_step_col){     // row step
        for(j=0, lr=0, rf=0; j < col; j++, lr+=lr_step, rf+=byte_step_col){ // column step
            lhs_idx = lr + lf;  // (j*row+i)*byteStep*Nch(=stepCol)
            rhs_idx = rr + rf;  // (i*col+j)*byteStep*Nch(=stepCol)

            memcpy(&tmA[lhs_idx], &dat_ptr[rhs_idx], byte_step_col);
        }
    }
    mA.reset(tmA, std::default_delete<uchar[]>());
    sync_data_ptr();

    i   = col;
    col = row;
    row = i;
    stepRow = col*stepCol;
}

void Mat::printMat(const std::string objname) {
    if(!objname.empty())
        fprintf(stdout,"object : %s \n", objname.c_str());
    if(isEmpty()){
        fprintf(stdout,"Empty matrix\n");
        return ;
    }
    const int32 bufsz = 2049;
    const int32 bsz   = 32;
    const double neg_max_double = -DBL_EPSILON ;
    const double pos_min_double =  DBL_EPSILON ;
    char buf[bufsz]="\0";
    char tmp[bsz];
    uint32 i,j,k, m;

    m = 0;
    if(datT==DTYP::CMPLX){
        cmplx val;
        for( i = 0; i < lenRowCol; i += col){ // rows
            snprintf(buf,bufsz,"[");
            for( j=0; j < col; j++){          // columns
                strncat(buf," (", 3);
                for( k = 0 ; k < Nch; k++){
                    val = elptr.cmx_ptr[m++];
                    if( val.re >= neg_max_double && val.re <= pos_min_double)
                        val.re = 0.0;
                    if( val.im >= neg_max_double && val.im <= pos_min_double)
                        val.im = 0.0;

                    snprintf(tmp,bsz,"% 6.2g %-+6.2gi",val.re, val.im);
                    strncat(buf,tmp,bufsz);
                    if( k < Nch-1)
                        strncat(buf, ",", 2);
                }
                strncat(buf,") ", 3);
            }
            strncat(buf,"]",2);
            fprintf(stdout,"%s\n",buf);
        }
    }else if(datT==DTYP::DOUBLE || datT==DTYP::FLOAT){
        double val;
        for( i = 0; i < length; i += stepRow){ // rows
            snprintf(buf,bufsz,"[");
            for( j=0 ; j < stepRow; j+= stepCol ){ // columns
                strncat(buf, " (", 3 );
                for( k =0; k < Nch; ++k){
                    if(datT == DTYP::FLOAT)
                        val = elptr.f32_ptr[m++];
                    else
                        val = elptr.f64_ptr[m++];

                    if( val >= neg_max_double && val <= pos_min_double)
                        val = 0.0;

                    snprintf(tmp, bsz,"% 6.3g",val);
                    strncat(buf, tmp, bsz);
                    if( k < Nch-1)
                        strncat(buf, ",", 2);
                }
                strncat(buf, ") ", 3);
            }
            strncat(buf,"]",2);
            fprintf(stdout,"%s\n",buf);
        }
    }else{
        int32 val;
        for( i = 0; i < length; i += stepRow){ // rows
            snprintf(buf,bufsz,"[");
            for( j=0; j < stepRow; j+= stepCol){ // columns
                strncat(buf," (", 3);
                for( k = 0 ; k < Nch; k++){
                    if(datT == DTYP::UCHAR)
                        val = elptr.uch_ptr[m++];
                    else
                        val = elptr.int_ptr[m++];

                    if( val > -10000 && val < 10000 )
                        snprintf(tmp,bsz,"% 5d",val);
                    else
                        snprintf(tmp,bsz,"% 5g",(double)val);

                    strncat(buf,tmp,bufsz);
                    if( k < Nch-1)
                        strncat(buf, ",", 2);
                }
                strncat(buf,") ",3);
            }
            strncat(buf,"]",2);
            fprintf(stdout,"%s\n",buf);
        }
    }
}
void Mat::printMat()  {
    printMat(obj_name);
}

Mat Mat::copy() const{
    Mat A(this->datT, row, col, Nch);
    uchar *pt_matdat  = A.getDataPtr<uchar>();

    std::copy(dat_ptr, dat_ptr + length*byteStep, pt_matdat);
    return A;
}

Mat Mat::ones(uint32 r, uint32 c, uint32 ch, DTYP dt){
    if( r == 0 || c == 0 || ch == 0){
        fprintf(stdout,"In ones method: arguments r , c or ch are to be larger than 0 ");
        return Mat();
    }

    Mat A(dt, r, c, ch);
    elemptr t_ptrs;
    t_ptrs.uch_ptr = A.getDataPtr();
    uint32 len     = A.getLength();
    uint32 k       = 0;
    switch(dt){
    case DTYP::DOUBLE : while( k < len ) { t_ptrs.f64_ptr[k++] = 1.0; } break;
    case DTYP::FLOAT  : while( k < len ) { t_ptrs.f32_ptr[k++] = 1.f; } break;
    case DTYP::INT    : while( k < len ) { t_ptrs.int_ptr[k++] = 1  ; } break;
    case DTYP::UCHAR  : while( k < len ) { t_ptrs.uch_ptr[k++] = 1  ; } break;
    case DTYP::CMPLX  : while( k < len ) { t_ptrs.cmx_ptr[k++] = cmplx(1.0); }
    }
    return A;
}

Mat Mat::zeros(uint32 r, uint32 c, uint32 ch, DTYP dt){
    if( r == 0 || c == 0 || ch == 0){
        fprintf(stderr,"In zeros method: arguments r , c or ch are to be larger than 0 ");
        return Mat();
    }

    Mat A(dt, r, c, ch);
    uchar* pt_dat = A.getDataPtr();

    memset(pt_dat, 0, A.getByteLen());
    return A;
}

Mat Mat::copyChannelN(const uint32 NoCh) const{
    uint32 numCh = NoCh;
    if( isEmpty() ){
        fprintf(stderr,"Current Mat has no data! So It cannot operate copyChannelN \n");

    }else if(numCh >= Nch ){
        fprintf(stdout,"channel number of copyChannelN() is out of bound\n The last channel is selected\n");
        numCh = Nch -1;
    }

    Mat A(this->datT, row,col,1);

    uint32 cplen  = lenRowCol*byteStep;
    uint32 offset = numCh*cplen;
    uchar *srcDat_pt = mA.get()+offset;
    uchar *tarDat_pt = A.mA.get();

    std::copy(srcDat_pt, srcDat_pt+cplen, tarDat_pt);
    /*
    for(uint32 i=offset; i < lenRowCol*byteStep; i++ )
        tarDat_pt[i] = srcDat_pt[i];
    */
    return A;
}

void Mat::setChannelN(const Mat& src, const uint32 srcFromCh,const uint32 Channels, const uint32 tarToCh){
    if(src.getChannel() < srcFromCh+Channels ){
        fprintf(stderr,"setChannelN(): srcFromCh and Channels are not correct! \n");
        return ;
    }

    uint32 tar_ch;
    if(isEmpty()){
        init(src.getRow(),src.getCol(), Channels, src.getDatType() );
        tar_ch = 0;
    }else if( Nch < tarToCh+Channels){
        fprintf(stderr,"setChannelN(): tarToCh and Channels are not correct! \n");
        return ;
    }else if( src.getRow() != row && src.getCol() != col){
        fprintf(stderr,"*this and src are not equal size. *this:(%d, %d, %d) != src(%d, %d, %d)\n",row,col,Nch,src.getRow(),src.getCol(),src.getChannel());
        return ;
    }else
        tar_ch = tarToCh;

    sync_data_ptr();
    uint32 cplen = src.getRowColSize() * src.getByteStep();
    uint32 start = cplen * srcFromCh;
    uint32 last  = cplen * (srcFromCh+Channels);
    uchar *srcdat_ptr_start = src.getMat().get() + start;
    uchar *srcdat_ptr_last  = src.getMat().get() + last ;
    uchar *tardat_ptr_start = mA.get() + cplen*tarToCh;
    std::copy( srcdat_ptr_start, srcdat_ptr_last, tardat_ptr_start );
/*
    uint32 src_idx = srcFromCh*lenRowCol*byteStep;
    uint32 tar_idx = tar_ch*lenRowCol*byteStep;    
    for(uint32 j=0; j < Channels ; j++){
        for(uint32 i=0; i < lenRowCol*byteStep; i++ ){
            dat_ptr[tar_idx++] = srcdat_ptr[src_idx++];
        }
    }
*/
}

void Mat::setName(std::string name){
    obj_name = name;
}

Mat Mat::copySubMat(const uint32 startRow, const uint32 endRow, const uint32 startCol, const uint32 endCol) const {
    if( startRow > row || endRow > row || startCol > col || endCol > col){
        fprintf(stderr,"copySubMat() : one or more arguments are out of bound from *this mat \n");
        return Mat();
    }else if( endRow < startRow){
        fprintf(stderr,"copySubMat() : endRow argument is less than startRow argument into copySubMat \n");
        return Mat();
    }else if( endCol < startCol){
        fprintf(stderr,"copySubMat() : endCol argument is less than startCol argument into copySubMat \n");
        return Mat();
    }

    uint32 new_row = endRow-startRow+1;
    uint32 new_col = endCol-startCol+1;
    Mat A( datT, new_row, new_col, Nch);
    uchar *tardat_ptr = A.getMat().get();

    uint32 ch_offset = 0 ;
    uint32 r, ch;
    uint32 offset;
    uint32 lenRCByteStep = lenRowCol*byteStep;
    uint32 rowByteStep   = col*byteStep;
    uint32 startColByte  = startCol*byteStep;
    uint32 endColByte    = (endCol+1)*byteStep;
    uint32 colBytes      = endColByte - startColByte;
    uint32 rowstart      = startRow*rowByteStep;
    uint32 rowend        = (endRow+1)*rowByteStep;
    uchar* k             = tardat_ptr;
    uchar* colstart;
    uchar* colend;
    for( ch=0, ch_offset=0; ch < Nch; ++ch, ch_offset += lenRCByteStep){
        for( r = rowstart; r < rowend; r+=rowByteStep ){
            offset   = ch_offset + r;
            colstart = dat_ptr + offset + startColByte;
            colend   = colstart+ endColByte;
            /*
            for( c = colstart; c < colend ; c++){
                tardat_ptr[k++] = dat_ptr[ c ];
            } */
            // =>
            std::copy(colstart,colend, k);
            k += colBytes;
        }
    }
    return A;
}

int32 Mat::sliceCopyMat(const Mat& src, const matRect& srcSlice,const Mat& des, const matRect& desSlice ){
    int32 srcR = src.getRow();
    int32 srcC = src.getCol();
    int32 desR = des.getRow();
    int32 desC = des.getCol();

    if( srcSlice.sR > srcR || srcSlice.eR > srcR || srcSlice.sC > srcC || srcSlice.eC > srcC){
        fprintf(stderr,"SliceCopyMat() : one or more arguments are out of bound from \'src mat\' \n");
        return -1;
    }else if( desSlice.sR > desR || desSlice.eR > desR || desSlice.sC > desC || desSlice.eC > desC){
        fprintf(stderr,"SliceCopyMat() : one or more arguments are out of bound from \'des mat\' \n");
        return -1;
    }else if( srcSlice.eR < srcSlice.sR || srcSlice.eC < srcSlice.sC){
        fprintf(stderr,"SliceCopyMat() : eR(or eC)of srcSlice has to be larger than sR(or sC) \n");
        return -1;
    }else if( desSlice.eR < desSlice.sR || desSlice.eC < desSlice.sC){
        fprintf(stderr,"SliceCopyMat() : eR(or eC)of desSlice has to be larger than sR(or sC) \n");
        return -1;
    }else if( (desSlice.eR-desSlice.sR != srcSlice.eR-srcSlice.sR) || (desSlice.eC-desSlice.sC != srcSlice.eC-srcSlice.sC)) {
        fprintf(stderr,"SliceCopyMat() : The sizes between srcSlice and desSlice are not matched! \n");
        return -1;
    }else if( des.getDatType() != src.getDatType()){
        fprintf(stderr,"SliceCopyMat() : The data types of srcSlice and desSlice are not matched! \n");
        return -1;
    }else if( des.getChannel() != src.getChannel()){
        fprintf(stderr,"SliceCopyMat() : The channels of srcSlice and desSlice are not matched! \n");
        return -1;
    }else if( src.isEmpty() || des.isEmpty()){
        fprintf(stderr,"SliceCopyMat() : At least one of src and des mat is empty! \n");
        return -1;
    }

    uint32 src_lenRowCol = srcR * srcC;
    uint32 des_lenRowCol = desR * desC;
    uint32 byteStep  = src.getByteStep();

    uchar *des_dat_ptr = des.getMat().get();
    uchar *src_dat_ptr = src.getMat().get();

    uint32 src_ch_offset, des_ch_offset ;
    uint32 src_r, des_r;
    uint32 ch;
    uint32 src_offset, des_offset;
    uint32 src_lenRCByteStep = src_lenRowCol*byteStep;
    uint32 src_rowByteStep   = srcC*byteStep;
    uint32 src_startColByte  = srcSlice.sC*byteStep;
    uint32 src_endColByte    = (srcSlice.eC+1)*byteStep;
    uint32 src_rowstart      = srcSlice.sR*src_rowByteStep;
    uint32 src_rowend        = (srcSlice.eR+1)*src_rowByteStep;
    uchar *src_colstart_p;
    uchar *src_colend_p;

    uint32 des_lenRCByteStep = des_lenRowCol*byteStep;
    uint32 des_rowByteStep   = desC*byteStep;
    uint32 des_startColByte  = desSlice.sC*byteStep;
    uint32 des_rowstart      = desSlice.sR*des_rowByteStep;
    uchar *des_colstart_p;
    src_ch_offset = 0;
    des_ch_offset = 0;
    for( ch=0; ch < src.getChannel(); ++ch){
        for( src_r = src_rowstart, des_r = des_rowstart; src_r < src_rowend; src_r+=src_rowByteStep, des_r+=des_rowByteStep ){
            src_offset     = src_ch_offset + src_r;
            src_colstart_p = src_dat_ptr + (src_offset + src_startColByte);
            src_colend_p   = src_dat_ptr + (src_offset + src_endColByte);

            des_offset     = des_ch_offset + des_r;
            des_colstart_p = des_dat_ptr + des_offset + des_startColByte;
            /*
            for( c = colstart; c < colend ; c++){
                tardat_ptr[k++] = dat_ptr[ c ];
            } */
            // =>
            std::copy(src_colstart_p, src_colend_p, des_colstart_p);
            //printf("%p, %p, %p\n",src_colstart_p, src_colend_p, des_colstart_p );
        }
        src_ch_offset += src_lenRCByteStep;
        des_ch_offset += des_lenRCByteStep;
    }
    return 1;
}

Mat& Mat::plusMat(const Mat& other){
    if(isEmpty() || other.isEmpty()) {
        fprintf(stderr, "Mat::plusMat method : Either of this or other is empty\n");
        return *this;
    }else if(row != other.getRow() || col != other.getCol() || Nch != other.getChannel()){
        fprintf(stderr, "Mat::plusMat method : The sizes of this and other are not the same!\n");
        return *this;
    }

    DTYP othrdatT = other.getDatType();
    if(othrdatT != datT){
        fprintf(stderr, "Waring, Mat::plusMat method : data types between self and other are different\n");
    }
    uint32 k = 0;
    elemptr slf_ptr, oth_ptr;
    slf_ptr.uch_ptr = dat_ptr;
    oth_ptr.uch_ptr = other.getDataPtr();

    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] += oth_ptr.f64_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] += oth_ptr.f64_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] += oth_ptr.f64_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] += oth_ptr.f64_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] += oth_ptr.f64_ptr[k]; }
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] += oth_ptr.f32_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] += oth_ptr.f32_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] += oth_ptr.f32_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] += oth_ptr.f32_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] += oth_ptr.f32_ptr[k]; }
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] += oth_ptr.int_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] += oth_ptr.int_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] += oth_ptr.int_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] += oth_ptr.int_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] += oth_ptr.int_ptr[k]; }
        }
    }else if (othrdatT== DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] += oth_ptr.uch_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] += oth_ptr.uch_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] += oth_ptr.uch_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] += oth_ptr.uch_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] += oth_ptr.uch_ptr[k]; }
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] += oth_ptr.cmx_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] += oth_ptr.cmx_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] += oth_ptr.cmx_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] += oth_ptr.cmx_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] += oth_ptr.cmx_ptr[k]; }
        }
    }
    return *this;
}

Mat& Mat::minusMat(const Mat& other){
    if(isEmpty() || other.isEmpty()) {
        fprintf(stdout, "Mat::minusMat method : either of this or other is empty\n");
        return *this;
    }else if(row != other.getRow() || col != other.getCol() || Nch != other.getChannel()){
        fprintf(stderr, "Mat::minusMat method : The sizes of this and other are not the same!\n");
        return *this;
    }

    DTYP othrdatT = other.getDatType();
    if(othrdatT != datT){
        fprintf(stderr, "Waring, Mat::minusMat method : data types between self and other are different\n");
    }
    uint32 k = 0;
    elemptr slf_ptr, oth_ptr;
    slf_ptr.uch_ptr = dat_ptr;
    oth_ptr.uch_ptr = other.getDataPtr();

    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] -= oth_ptr.f64_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] -= oth_ptr.f64_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] -= oth_ptr.f64_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] -= oth_ptr.f64_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] -= oth_ptr.f64_ptr[k]; }
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] -= oth_ptr.f32_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] -= oth_ptr.f32_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] -= oth_ptr.f32_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] -= oth_ptr.f32_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] -= oth_ptr.f32_ptr[k]; }
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] -= oth_ptr.int_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] -= oth_ptr.int_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] -= oth_ptr.int_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] -= oth_ptr.int_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] -= oth_ptr.int_ptr[k]; }
        }
    }else if (othrdatT== DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] -= oth_ptr.uch_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] -= oth_ptr.uch_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] -= oth_ptr.uch_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] -= oth_ptr.uch_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] -= oth_ptr.uch_ptr[k]; }
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] -= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] -= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] -= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] -= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] -= oth_ptr.cmx_ptr[k]; }
        }
    }
    return *this;
}

Mat& Mat::mulMat(const Mat& other){
    if(isEmpty() || other.isEmpty()) {
        fprintf(stdout, "Mat::mulMat method : either of this or other is empty\n");
        return *this;
    }else if(row != other.getRow() || col != other.getCol() || Nch != other.getChannel()){
        fprintf(stderr, "Mat::mulMat method : The sizes of this and other are not the same!\n");
        return *this;
    }

    DTYP othrdatT = other.getDatType();
    if(othrdatT != datT){
        fprintf(stderr, "Waring, Mat::mulMat method : data types between self and other are different\n");
    }
    uint32 k = 0;
    elemptr slf_ptr, oth_ptr;
    slf_ptr.uch_ptr = dat_ptr;
    oth_ptr.uch_ptr = other.getDataPtr();

    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] *= oth_ptr.f64_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] *= oth_ptr.f64_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] *= oth_ptr.f64_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] *= oth_ptr.f64_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] *= oth_ptr.f64_ptr[k]; }
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] *= oth_ptr.f32_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] *= oth_ptr.f32_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] *= oth_ptr.f32_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] *= oth_ptr.f32_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] *= oth_ptr.f32_ptr[k]; }
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] *= oth_ptr.int_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] *= oth_ptr.int_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] *= oth_ptr.int_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] *= oth_ptr.int_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] *= oth_ptr.int_ptr[k]; }
        }
    }else if (othrdatT== DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] *= oth_ptr.uch_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] *= oth_ptr.uch_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] *= oth_ptr.uch_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] *= oth_ptr.uch_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] *= oth_ptr.uch_ptr[k]; }
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] *= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] *= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] *= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] *= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] *= oth_ptr.cmx_ptr[k]; }
        }
    }
    return *this;
}

Mat& Mat::divMat(const Mat& other){
    if(isEmpty() || other.isEmpty()) {
        fprintf(stdout, "Mat::divMat method : either of this or other is empty\n");
        return *this;
    }else if(row != other.getRow() || col != other.getCol() || Nch != other.getChannel()){
        fprintf(stderr, "Mat::divMat method : The sizes of this and other are not the same!\n");
        return *this;
    }

    DTYP othrdatT = other.getDatType();
    if(othrdatT != datT){
        fprintf(stderr, "Waring, Mat::divMat method : data types between self and other are different\n");
    }
    uint32 k = 0;
    elemptr slf_ptr, oth_ptr;
    slf_ptr.uch_ptr = dat_ptr;
    oth_ptr.uch_ptr = other.getDataPtr();

    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] /= oth_ptr.f64_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] /= oth_ptr.f64_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] /= oth_ptr.f64_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] /= oth_ptr.f64_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] /= oth_ptr.f64_ptr[k]; }
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] /= oth_ptr.f32_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] /= oth_ptr.f32_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] /= oth_ptr.f32_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] /= oth_ptr.f32_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] /= oth_ptr.f32_ptr[k]; }
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] /= oth_ptr.int_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] /= oth_ptr.int_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] /= oth_ptr.int_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] /= oth_ptr.int_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] /= oth_ptr.int_ptr[k]; }
        }
    }else if (othrdatT== DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] /= oth_ptr.uch_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] /= oth_ptr.uch_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] /= oth_ptr.uch_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] /= oth_ptr.uch_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] /= oth_ptr.uch_ptr[k]; }
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { slf_ptr.f64_ptr[k] /= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { slf_ptr.f32_ptr[k] /= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { slf_ptr.int_ptr[k] /= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { slf_ptr.uch_ptr[k] /= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { slf_ptr.cmx_ptr[k] /= oth_ptr.cmx_ptr[k]; }
        }
    }
    return *this;
}

Mat& Mat::plusScalar(const double scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] += scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] += scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] += scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] += scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] += scalar; }
    }
    return *this;
}
Mat& Mat::plusScalar(const float scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] += scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] += scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] += scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] += scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] += scalar; }
    }
    return *this;
}
Mat& Mat::plusScalar(const int32 scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] += scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] += scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] += scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] += scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] += scalar; }
    }
    return *this;
}
Mat& Mat::plusScalar(const uchar scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] += scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] += scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] += scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] += scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] += scalar; }
    }
    return *this;
}

Mat& Mat::minusScalar(const double scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] -= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] -= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] -= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] -= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] -= scalar; }
    }
    return *this;
}

Mat& Mat::minusScalar(const float scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] -= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] -= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] -= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] -= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] -= scalar; }
    }
    return *this;
}

Mat& Mat::minusScalar(const int32 scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] -= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] -= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] -= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] -= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] -= scalar; }
    }
    return *this;
}
Mat& Mat::minusScalar(const uchar scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] -= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] -= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] -= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] -= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] -= scalar; }
    }
    return *this;
}

Mat& Mat::mulScalar(const double scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] *= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] *= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] *= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] *= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] *= scalar; }
    }
    return *this;
}

Mat& Mat::mulScalar(const float scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] *= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] *= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] *= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] *= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] *= scalar; }
    }
    return *this;
}

Mat& Mat::mulScalar(const int32 scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] *= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] *= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] *= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] *= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] *= scalar; }
    }
    return *this;
}
Mat& Mat::mulScalar(const uchar scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] *= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] *= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] *= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] *= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] *= scalar; }
    }
    return *this;
}

Mat& Mat::divByScalar(const double scalar){
    if(isEmpty()) return *this;
    assert(scalar != 0.0);

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] /= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] /= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] /= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] /= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] /= scalar; }
    }
    return *this;
}

Mat& Mat::divByScalar(const float scalar){
    if(isEmpty()) return *this;
    assert(scalar == 0.0);

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] /= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] /= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] /= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] /= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] /= scalar; }
    }
    return *this;
}

Mat& Mat::divByScalar(const int32 scalar){
    if(isEmpty()) return *this;
    assert(scalar == 0.0);

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] /= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] /= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] /= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] /= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] /= scalar; }
    }
    return *this;
}
Mat& Mat::divByScalar(const uchar scalar){
    if(isEmpty()) return *this;
    assert(scalar == 0.0);

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : while( k < length ) { ptr.f64_ptr[k++] /= scalar; } break;
    case DTYP::FLOAT  : while( k < length ) { ptr.f32_ptr[k++] /= scalar; } break;
    case DTYP::INT    : while( k < length ) { ptr.int_ptr[k++] /= scalar; } break;
    case DTYP::UCHAR  : while( k < length ) { ptr.uch_ptr[k++] /= scalar; } break;
    case DTYP::CMPLX  : while( k < length ) { ptr.cmx_ptr[k++] /= scalar; }
    }
    return *this;
}

Mat& Mat::divScalar(const double scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : for(; k < length; ++k ) { ptr.f64_ptr[k] = scalar / ptr.f64_ptr[k]; } break;
    case DTYP::FLOAT  : for(; k < length; ++k ) { ptr.f32_ptr[k] = scalar / ptr.f32_ptr[k]; } break;
    case DTYP::INT    : for(; k < length; ++k ) { ptr.int_ptr[k] = scalar / ptr.int_ptr[k]; } break;
    case DTYP::UCHAR  : for(; k < length; ++k ) { ptr.uch_ptr[k] = scalar / ptr.uch_ptr[k]; } break;
    case DTYP::CMPLX  : for(; k < length; ++k ) { ptr.cmx_ptr[k] = scalar / ptr.cmx_ptr[k]; }
    }
    return *this;
}

Mat& Mat::divScalar(const float scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : for(; k < length; ++k ) { ptr.f64_ptr[k] = scalar / ptr.f64_ptr[k]; } break;
    case DTYP::FLOAT  : for(; k < length; ++k ) { ptr.f32_ptr[k] = scalar / ptr.f32_ptr[k]; } break;
    case DTYP::INT    : for(; k < length; ++k ) { ptr.int_ptr[k] = scalar / ptr.int_ptr[k]; } break;
    case DTYP::UCHAR  : for(; k < length; ++k ) { ptr.uch_ptr[k] = scalar / ptr.uch_ptr[k]; } break;
    case DTYP::CMPLX  : for(; k < length; ++k ) { ptr.cmx_ptr[k] = scalar / ptr.cmx_ptr[k]; }
    }
    return *this;
}

Mat& Mat::divScalar(const int32 scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : for(; k < length; ++k ) { ptr.f64_ptr[k] = scalar / ptr.f64_ptr[k]; } break;
    case DTYP::FLOAT  : for(; k < length; ++k ) { ptr.f32_ptr[k] = scalar / ptr.f32_ptr[k]; } break;
    case DTYP::INT    : for(; k < length; ++k ) { ptr.int_ptr[k] = scalar / ptr.int_ptr[k]; } break;
    case DTYP::UCHAR  : for(; k < length; ++k ) { ptr.uch_ptr[k] = scalar / ptr.uch_ptr[k]; } break;
    case DTYP::CMPLX  : for(; k < length; ++k ) { ptr.cmx_ptr[k] = scalar / ptr.cmx_ptr[k]; }
    }
    return *this;
}
Mat& Mat::divScalar(const uchar scalar){
    if(isEmpty()) return *this;

    uint32 k = 0;
    elemptr ptr;
    ptr.uch_ptr = dat_ptr;
    switch(datT){
    case DTYP::DOUBLE : for(; k < length; ++k ) { ptr.f64_ptr[k] = scalar / ptr.f64_ptr[k]; } break;
    case DTYP::FLOAT  : for(; k < length; ++k ) { ptr.f32_ptr[k] = scalar / ptr.f32_ptr[k]; } break;
    case DTYP::INT    : for(; k < length; ++k ) { ptr.int_ptr[k] = scalar / ptr.int_ptr[k]; } break;
    case DTYP::UCHAR  : for(; k < length; ++k ) { ptr.uch_ptr[k] = scalar / ptr.uch_ptr[k]; } break;
    case DTYP::CMPLX  : for(; k < length; ++k ) { ptr.cmx_ptr[k] = scalar / ptr.cmx_ptr[k]; }
    }
    return *this;
}

Mat Mat::max(){
    if(isEmpty()) return Mat();

    uint32 ch = getChannel();
    uint32 k, m, n;
    elemptr Aptrs;
    Mat  A(getDatType(),1,1,ch);
    Aptrs.uch_ptr = A.getDataPtr();
    k=0;
    if(datT == DTYP::DOUBLE){
        for(; k<ch; ++k) { Aptrs.f64_ptr[k]= elptr.f64_ptr[k]; }
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                if( Aptrs.f64_ptr[k] < elptr.f64_ptr[n])
                    Aptrs.f64_ptr[k] = elptr.f64_ptr[n];
            }
        }
    }else if(datT == DTYP::FLOAT){
        for(; k<ch; ++k) { Aptrs.f32_ptr[k]= elptr.f32_ptr[k]; }
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                if( Aptrs.f32_ptr[k] < elptr.f32_ptr[n])
                    Aptrs.f32_ptr[k] = elptr.f32_ptr[n];
            }
        }
    }else if(datT == DTYP::INT){
        for(; k<ch; ++k) { Aptrs.int_ptr[k]= elptr.int_ptr[k]; }
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                if( Aptrs.int_ptr[k] < elptr.int_ptr[n])
                    Aptrs.int_ptr[k] = elptr.int_ptr[n];
            }
        }
    }else if(datT == DTYP::UCHAR){
        for(; k<ch; ++k) { Aptrs.uch_ptr[k]= elptr.uch_ptr[k]; }
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                if( Aptrs.uch_ptr[k] < elptr.uch_ptr[n])
                    Aptrs.uch_ptr[k] = elptr.uch_ptr[n];
            }
        }
    }else if(datT == DTYP::CMPLX){
        cmplx   cl, tmp;
        double *clmag_ch, large_mag, tmp_mag;
        clmag_ch = new double[ch];
        for(; k<ch; ++k) { cl = elptr.cmx_ptr[k]; Aptrs.cmx_ptr[k]= cl; clmag_ch[k] = cl.square();}
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                tmp       = elptr.cmx_ptr[n];
                tmp_mag   = tmp.square();
                large_mag = clmag_ch[k];
                if( clmag_ch[k] < tmp_mag){
                    Aptrs.cmx_ptr[k] = tmp;
                    clmag_ch[k]      = tmp_mag;
                }
            }
        }
        delete [] clmag_ch;
    }
    return A;
}

Mat Mat::min(){
    if(isEmpty()) return Mat();

    uint32 ch = getChannel();
    uint32 k, m, n;
    elemptr Aptrs;
    Mat  A(getDatType(),1,1,ch);
    Aptrs.uch_ptr = A.getDataPtr();
    k=0;
    if(datT == DTYP::DOUBLE){
        for(; k<ch; ++k) { Aptrs.f64_ptr[k]= elptr.f64_ptr[k]; }
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                if( Aptrs.f64_ptr[k] > elptr.f64_ptr[n])
                    Aptrs.f64_ptr[k] = elptr.f64_ptr[n];
            }
        }
    }else if(datT == DTYP::FLOAT){
        for(; k<ch; ++k) { Aptrs.f32_ptr[k]= elptr.f32_ptr[k]; }
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                if( Aptrs.f32_ptr[k] > elptr.f32_ptr[n])
                    Aptrs.f32_ptr[k] = elptr.f32_ptr[n];
            }
        }
    }else if(datT == DTYP::INT){
        for(; k<ch; ++k) { Aptrs.int_ptr[k]= elptr.int_ptr[k]; }
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                if( Aptrs.int_ptr[k] > elptr.int_ptr[n])
                    Aptrs.int_ptr[k] = elptr.int_ptr[n];
            }
        }
    }else if(datT == DTYP::UCHAR){
        for(; k<ch; ++k) { Aptrs.uch_ptr[k]= elptr.uch_ptr[k]; }
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                if( Aptrs.uch_ptr[k] > elptr.uch_ptr[n])
                    Aptrs.uch_ptr[k] = elptr.uch_ptr[n];
            }
        }
    }else if(datT == DTYP::CMPLX){
        cmplx   cl, tmp;
        double *clmag_ch, large_mag, tmp_mag;
        clmag_ch = new double[ch];
        for(; k<ch; ++k) { cl = elptr.cmx_ptr[k]; Aptrs.cmx_ptr[k]= cl; clmag_ch[k] = cl.square(); }
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                tmp       = elptr.cmx_ptr[n];
                tmp_mag   = tmp.square();
                large_mag = clmag_ch[k];
                if( clmag_ch[k] > tmp_mag){
                    Aptrs.cmx_ptr[k] = tmp;
                    clmag_ch[k]      = tmp_mag;
                }
            }
        }
        delete [] clmag_ch;
    }
    return A;
}

Mat Mat::mean(){
    if(isEmpty()) return Mat();

    return (sum()/getRowColSize());
}

Mat Mat::sum(){
    if(isEmpty()) return Mat();
    uint32 ch = getChannel();
    uint32 k, m, n;
    n=0;
    if(datT == DTYP::CMPLX){
        Mat  B = Mat::zeros(1,1,ch,DTYP::CMPLX);
        cmplx* Bptr = B.getDataPtr<cmplx>();

        for(m = 0 ; m < length; m+=ch ){
            for(k=0 ; k < ch; ++k, ++n)
                Bptr[k] += elptr.cmx_ptr[n];
        }
        return B;
    }else{
        Mat  A = Mat::zeros(1,1,ch,DTYP::DOUBLE);
        double* Aptr = A.getDataPtr<double>();

        if(datT == DTYP::DOUBLE){
            for(m = 0 ; m < length; m+=ch ){
                for( k=0 ; k < ch; ++k, ++n)
                    Aptr[k] += elptr.f64_ptr[n];
            }
        }else if(datT == DTYP::FLOAT){
            for(m = 0 ; m < length; m+=ch ){
                for( k=0 ; k < ch; ++k, ++n)
                    Aptr[k] += double(elptr.f32_ptr[n]);
            }
        }else if(datT == DTYP::INT){
            for(m = 0 ; m < length; m+=ch ){
                for( k=0 ; k < ch; ++k, ++n)
                    Aptr[k] += double(elptr.int_ptr[n]);
            }
        }else if(datT == DTYP::UCHAR){
            for(m = 0 ; m < length; m+=ch ){
                for( k=0 ; k < ch; ++k, ++n)
                    Aptr[k] += double(elptr.int_ptr[n]);
            }
        }
        return A;
    }
}

Mat Mat::var(){
    if(isEmpty()) return Mat();

    uint32 k, m, n;
    uint32 ch = getChannel();
    Mat avg   = mean();
    Mat A     = Mat::zeros(1,1,ch,DTYP::DOUBLE);
    uint32 Div= getRowColSize() -1;
    n = 0;
    if(datT==DTYP::CMPLX){
        // STD of complex numbers
        // E[ |Z - E[Z]|^2 ] = E[|Z|^2] - |E[Z]|^2 --> real valued Mat
        // https://en.wikipedia.org/wiki/Complex_random_variable#Expectation
        cmplx diff;
        for(m = 0 ; m < length; m+=ch){
            for(k=0; k < ch; ++k, ++n){
                diff = elptr.cmx_ptr[n] - avg.elptr.cmx_ptr[k];
                A.at<double>(k) += diff.square();
            }
        }
    }else{
        double diff;
        if(datT == DTYP::DOUBLE){
            for(m = 0 ; m < length; m+=ch){
                for(k=0; k < ch; ++k, ++n){
                    diff = double(elptr.f64_ptr[n]) - avg.at<double>(k);
                    A.at<double>(k) += diff*diff;
                }
            }
        }else if(datT == DTYP::FLOAT){
            for(m = 0 ; m < length; m+=ch){
                for(k=0; k < ch; ++k, ++n){
                    diff = double(elptr.f32_ptr[n]) - avg.at<double>(k);
                    A.at<double>(k) += diff*diff;
                }
            }
        }else if(datT == DTYP::INT){
            for(m = 0 ; m < length; m+=ch){
                for(k=0; k < ch; ++k, ++n){
                    diff = double(elptr.int_ptr[n]) - avg.at<double>(k);
                    A.at<double>(k) += diff*diff;
                }
            }
        }else if(datT == DTYP::UCHAR){
            for(m = 0 ; m < length; m+=ch){
                for(k=0; k < ch; ++k, ++n){
                    diff = double(elptr.uch_ptr[n]) - avg.at<double>(k);
                    A.at<double>(k) += diff*diff;
                }
            }
        }
    }
    A /= Div; // Varaince

    return A;
}

Mat Mat::std(){
    if(isEmpty()) return Mat();

    Mat V = var();
    for(uint32 k=0; k < V.getChannel(); ++k)  // standard deviation
        V.at<double>(k) = sqrt(V.at<double>(k));
    return V;
}

Mat Mat::sqrtm(){
    if(isEmpty()) return Mat();

    uint32 k=0;
    if( datT == DTYP::CMPLX ){
        Mat V = Mat::zeros(row, col, Nch, DTYP::CMPLX);
        for(; k < V.length; ++k)
            V.at<cmplx>(k) = elptr.cmx_ptr[k].sqrtc();
        return V;
    }else{
        Mat V = Mat::zeros(row, col, Nch, DTYP::DOUBLE);
        if( datT == DTYP::DOUBLE ){
            for(; k < V.length; ++k)  V.at<double>(k) = sqrt(elptr.f64_ptr[k]);
        }else if( datT == DTYP::FLOAT){
            for(; k < V.length; ++k)  V.at<double>(k) = sqrt(elptr.f32_ptr[k]);
        }else if( datT == DTYP::INT  ){
            for(; k < V.length; ++k)  V.at<double>(k) = sqrt(elptr.int_ptr[k]);
        }else if( datT == DTYP::UCHAR){
            for(; k < V.length; ++k)  V.at<double>(k) = sqrt(elptr.uch_ptr[k]);
        }
        return V;
    }
}

Mat Mat::repeat(const Mat& src, const uint32 rr, const uint32 rc, const uint32 rch){
    if(src.isEmpty()) return Mat();

    if(src.getChannel() > 1){
        fprintf(stderr,"src Mat argument of repeat method should have only 1 channel\n");
        return Mat();
    }

    return Mat();
    //switch(src.getDatType()){
    //case DTYP::DOUBLE : return _repeat<double>(src, rr, rc, rch);
    //case DTYP::FLOAT  : return _repeat<float >(src, rr, rc, rch);
    //case DTYP::INT    : return _repeat<int32 >(src, rr, rc, rch);
    //case DTYP::UCHAR  : return _repeat<uchar >(src, rr, rc, rch);
    //default           : return _repeat<cmplx >(src, rr, rc, rch); // case DTYP::CMPLX
    //}
}

} // end of jmat namespace
