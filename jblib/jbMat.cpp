/*
 * Copyright (C) 2020. Jong B. Choi
 * License : MIT License
 * contact : neoflame99@naver.com
 */

#include <stdio.h>
#include "jbMat.h"
#include <iostream>
#include <utility>

namespace jmat {



void Mat::init(const uint32 r,const uint32 c,const uint32 ch,const DTYP dt,const bool do_alloc){
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
Mat::Mat(){
    init(0,0,0,DTYP::UCHAR,true);
}
Mat::Mat(const DTYP dt, const uint32 rc ):mA(nullptr){
    init(rc, rc, 1, dt, true);
}
Mat::Mat(const DTYP dt, const uint32 r, const uint32 c, const uint32 ch):mA(nullptr){
    init(r, c, ch, dt, true);
}
Mat::Mat(const shr_ptr ma, const DTYP dt, const uint32 r, const uint32 c, const uint32 ch):mA(ma){
    init(r, c, ch, dt, false);
}
Mat::Mat(const DTYP dt, const uint32 r, const uint32 c, const uint32 ch, const std::string name):mA(nullptr){
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

Mat::Mat( const std::initializer_list<double> list ){
    //-- Making a column vector
    init(list.size(),1,1,DTYP::DOUBLE, true);

    if(elptr.f64_ptr !=nullptr){
        memcpy(elptr.f64_ptr, list.begin(), list.size()*sizeof(double));
    }
}
Mat::Mat( const std::initializer_list<float> list ){
    //-- Making a column vector
    init(list.size(),1,1,DTYP::FLOAT, true);

    if(elptr.f32_ptr!=nullptr){
        memcpy(elptr.f32_ptr, list.begin(), list.size()*sizeof(float));
    }
}
Mat::Mat( const std::initializer_list<int32> list ){
    //-- Making a column vector
    init(list.size(),1,1,DTYP::INT, true);

    if(elptr.int_ptr!=nullptr){
        memcpy(elptr.int_ptr, list.begin(), list.size()*sizeof(int32));
    }
}
Mat::Mat( const std::initializer_list<uchar> list ){
    //-- Making a column vector
    init(list.size(),1,1,DTYP::UCHAR, true);

    if(elptr.uch_ptr!=nullptr){
        memcpy(elptr.uch_ptr, list.begin(), list.size()*sizeof(uchar));
    }
}
Mat::Mat( const std::initializer_list<cmplx> list ){
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

void Mat::setRowCol(const uint32 r,const uint32 c,const uint32 ch){
    uint32 len = r*c*ch;
    if(this->length != len){
        mA.reset();
        init(r,c,ch,datT,true);
    }else{
        init(r,c,ch,datT,false);
    }
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
        fprintf(stdout,"Mat::ones() : In ones method: arguments r , c or ch are to be larger than 0 ");
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
        fprintf(stderr,"Mat::zeros() : In zeros method: arguments r , c or ch are to be larger than 0 ");
        return Mat();
    }

    Mat A(dt, r, c, ch);
    uchar* pt_dat = A.getDataPtr();

    memset(pt_dat, 0, A.getByteLen());
    return A;
}

Mat Mat::extractChannel(const Mat& src, const uint32 ch){
    if( src.isEmpty() ){
        fprintf(stderr,"Mat::extractChannel() : Source Mat is empty! \n");
        return Mat();
    }else if(ch >= src.getChannel() ){
        fprintf(stderr,"Mat::extractChannel() : The argument ch is out of channel numbers of src Mat\n");
        return Mat();
    }

    Mat A(src.datT, src.getRow(), src.getCol(), 1);
    uint32 k,n ;
    switch(src.getDatType()){
    case DTYP::DOUBLE: for(k=0, n=ch; k< A.getLength(); ++k, n+=src.getChannel()){ A.elptr.f64_ptr[k]= src.elptr.f64_ptr[n]; } break;
    case DTYP::FLOAT : for(k=0, n=ch; k< A.getLength(); ++k, n+=src.getChannel()){ A.elptr.f32_ptr[k]= src.elptr.f32_ptr[n]; } break;
    case DTYP::INT   : for(k=0, n=ch; k< A.getLength(); ++k, n+=src.getChannel()){ A.elptr.int_ptr[k]= src.elptr.int_ptr[n]; } break;
    case DTYP::UCHAR : for(k=0, n=ch; k< A.getLength(); ++k, n+=src.getChannel()){ A.elptr.uch_ptr[k]= src.elptr.uch_ptr[n]; } break;
    case DTYP::CMPLX : for(k=0, n=ch; k< A.getLength(); ++k, n+=src.getChannel()){ A.elptr.cmx_ptr[k]= src.elptr.cmx_ptr[n]; }
    }
    return A;
}

void Mat::setChannel(const Mat& src, const uint32 srcCh, const uint32 tarCh, const uint32 Channels){
    if(isEmpty()){
        fprintf(stderr,"Mat::setChannel(): this Mat is empty! \n");
        return ;
    }else if(src.getChannel() < srcCh+Channels ){
        fprintf(stderr,"Mat::setChannel() : srcCh and Channels are not correct! \n");
        return ;
    }else if( Nch < tarCh+Channels){
        fprintf(stderr,"Mat::setChannel(): tarToCh and Channels are not correct! \n");
        return ;
    }else if( src.getRow() != row && src.getCol() != col){
        fprintf(stderr,"Mat::setChannel(): this and src are not equal size. *this:(%d, %d, %d) != src(%d, %d, %d)\n",row,col,Nch,src.getRow(),src.getCol(),src.getChannel());
        return ;
    }else if( src.getDatType() != datT){
        fprintf(stderr,"Mat::setChannel(): data types between src and this are not equal! \n");
        return ;
    }

    uint32 k, n, i;
    for(i=0 ; i < Channels; ++i){
        k = srcCh +i;
        n = tarCh +i;
        switch(datT){
        case DTYP::DOUBLE: 	for(; k < src.getLength(); k+=src.getChannel(), n+=this->Nch){ elptr.f64_ptr[n] = src.elptr.f64_ptr[k]; } break;
        case DTYP::FLOAT : 	for(; k < src.getLength(); k+=src.getChannel(), n+=this->Nch){ elptr.f32_ptr[n] = src.elptr.f32_ptr[k]; } break;
        case DTYP::INT   : 	for(; k < src.getLength(); k+=src.getChannel(), n+=this->Nch){ elptr.int_ptr[n] = src.elptr.int_ptr[k]; } break;
        case DTYP::UCHAR : 	for(; k < src.getLength(); k+=src.getChannel(), n+=this->Nch){ elptr.uch_ptr[n] = src.elptr.uch_ptr[k]; } break;
        case DTYP::CMPLX : 	for(; k < src.getLength(); k+=src.getChannel(), n+=this->Nch){ elptr.cmx_ptr[n] = src.elptr.cmx_ptr[k]; }
        }
    }
}

void Mat::setName(const std::string name){
    obj_name = name;
}

Mat Mat::copySubMat(const uint32 startRow, const uint32 endRow, const uint32 startCol, const uint32 endCol) const {
    if( startRow > row || endRow > row || startCol > col || endCol > col){
        fprintf(stderr,"Mat::copySubMat() : one or more arguments are out of bound from *this mat \n");
        return Mat();
    }else if( endRow < startRow){
        fprintf(stderr,"Mat::copySubMat() : endRow argument is less than startRow argument into copySubMat \n");
        return Mat();
    }else if( endCol < startCol){
        fprintf(stderr,"Mat::copySubMat() : endCol argument is less than startCol argument into copySubMat \n");
        return Mat();
    }

    uint32 new_row = endRow-startRow+1;
    uint32 new_col = endCol-startCol+1;
    Mat A( datT, new_row, new_col, Nch);

    elemptr sRowPtr, tPtr;
    uint32 m, n;
    uint32 start_p = startCol*stepCol*byteStep;
    uint32 cp_byte_len = byteStep*new_col*stepCol;

    tPtr.uch_ptr = A.getDataPtr();
    for(m=startRow, n=0; m <= endRow; ++m, ++n){
        sRowPtr = getRowElptr(m);
        tPtr    = A.getRowElptr(n);
        memcpy(tPtr.uch_ptr, sRowPtr.uch_ptr+start_p, cp_byte_len);
    }

    return A;
}

bool Mat::sliceCopyMat(const Mat& src, const matRect& srcSlice,const Mat& des, const matRect& desSlice ){
    int32 srcR = src.getRow();
    int32 srcC = src.getCol();
    int32 desR = des.getRow();
    int32 desC = des.getCol();

    if( src.isEmpty() || des.isEmpty()){
        fprintf(stderr,"Mat::sliceCopyMat() : At least one of src and des mat is empty! \n");
        return -1;
    }else if( des.getDatType() != src.getDatType()){
        fprintf(stderr,"Mat::sliceCopyMat() : The data types of srcSlice and desSlice are not matched! \n");
        return -1;
    }else if( des.getChannel() != src.getChannel()){
        fprintf(stderr,"Mat::sliceCopyMat() : The channels of srcSlice and desSlice are not matched! \n");
        return -1;
    }else if( srcSlice.sR > srcR || srcSlice.eR > srcR || srcSlice.sC > srcC || srcSlice.eC > srcC){
        fprintf(stderr,"Mat::sliceCopyMat() : one or more arguments are out of bound from \'src mat\' \n");
        return -1;
    }else if( desSlice.sR > desR || desSlice.eR > desR || desSlice.sC > desC || desSlice.eC > desC){
        fprintf(stderr,"Mat::sliceCopyMat() : one or more arguments are out of bound from \'des mat\' \n");
        return -1;
    }else if( srcSlice.eR < srcSlice.sR || srcSlice.eC < srcSlice.sC){
        fprintf(stderr,"Mat::sliceCopyMat() : eR(or eC)of srcSlice has to be larger than sR(or sC) \n");
        return -1;
    }else if( desSlice.eR < desSlice.sR || desSlice.eC < desSlice.sC){
        fprintf(stderr,"Mat::sliceCopyMat() : eR(or eC)of desSlice has to be larger than sR(or sC) \n");
        return -1;
    }else if( (desSlice.eR-desSlice.sR != srcSlice.eR-srcSlice.sR) || (desSlice.eC-desSlice.sC != srcSlice.eC-srcSlice.sC)) {
        fprintf(stderr,"Mat::sliceCopyMat() : The sizes between srcSlice and desSlice are not matched! \n");
        return -1;
    }

    elemptr sRowPtr, dRowPtr;
    uint32 m, n;
    uint32 src_col_start_p = srcSlice.sC*src.stepCol*src.byteStep;
    uint32 des_col_start_p = desSlice.sC*des.stepCol*des.byteStep;
    uint32 cp_byte_len = (srcSlice.eC - srcSlice.sC +1)*src.stepCol*src.byteStep;

    dRowPtr.uch_ptr = des.getDataPtr();
    for(m=srcSlice.sR, n=desSlice.sR; m <= (uint32)srcSlice.eR; ++m, ++n){
        sRowPtr = src.getRowElptr(m);
        dRowPtr = des.getRowElptr(n);
        memcpy(dRowPtr.uch_ptr+des_col_start_p, sRowPtr.uch_ptr+src_col_start_p, cp_byte_len);
    }

    return true;
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
    elemptr oth_ptr;
    oth_ptr.uch_ptr = other.getDataPtr();

    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] += oth_ptr.f64_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] += oth_ptr.f64_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] += oth_ptr.f64_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] += oth_ptr.f64_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] += oth_ptr.f64_ptr[k]; }
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] += oth_ptr.f32_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] += oth_ptr.f32_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] += oth_ptr.f32_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] += oth_ptr.f32_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] += oth_ptr.f32_ptr[k]; }
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] += oth_ptr.int_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] += oth_ptr.int_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] += oth_ptr.int_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] += oth_ptr.int_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] += oth_ptr.int_ptr[k]; }
        }
    }else if (othrdatT== DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] += oth_ptr.uch_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] += oth_ptr.uch_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] += oth_ptr.uch_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] += oth_ptr.uch_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] += oth_ptr.uch_ptr[k]; }
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] += oth_ptr.cmx_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] += oth_ptr.cmx_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] += oth_ptr.cmx_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] += oth_ptr.cmx_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] += oth_ptr.cmx_ptr[k]; }
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
    elemptr oth_ptr;
    oth_ptr.uch_ptr = other.getDataPtr();

    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] -= oth_ptr.f64_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] -= oth_ptr.f64_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] -= oth_ptr.f64_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] -= oth_ptr.f64_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] -= oth_ptr.f64_ptr[k]; }
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] -= oth_ptr.f32_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] -= oth_ptr.f32_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] -= oth_ptr.f32_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] -= oth_ptr.f32_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] -= oth_ptr.f32_ptr[k]; }
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] -= oth_ptr.int_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] -= oth_ptr.int_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] -= oth_ptr.int_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] -= oth_ptr.int_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] -= oth_ptr.int_ptr[k]; }
        }
    }else if (othrdatT== DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] -= oth_ptr.uch_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] -= oth_ptr.uch_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] -= oth_ptr.uch_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] -= oth_ptr.uch_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] -= oth_ptr.uch_ptr[k]; }
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] -= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] -= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] -= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] -= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] -= oth_ptr.cmx_ptr[k]; }
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
    elemptr oth_ptr;
    oth_ptr.uch_ptr = other.getDataPtr();

    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] *= oth_ptr.f64_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] *= oth_ptr.f64_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] *= oth_ptr.f64_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] *= oth_ptr.f64_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] *= oth_ptr.f64_ptr[k]; }
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] *= oth_ptr.f32_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] *= oth_ptr.f32_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] *= oth_ptr.f32_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] *= oth_ptr.f32_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] *= oth_ptr.f32_ptr[k]; }
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] *= oth_ptr.int_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] *= oth_ptr.int_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] *= oth_ptr.int_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] *= oth_ptr.int_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] *= oth_ptr.int_ptr[k]; }
        }
    }else if (othrdatT== DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] *= oth_ptr.uch_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] *= oth_ptr.uch_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] *= oth_ptr.uch_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] *= oth_ptr.uch_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] *= oth_ptr.uch_ptr[k]; }
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] *= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] *= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] *= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] *= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] *= oth_ptr.cmx_ptr[k]; }
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
    elemptr oth_ptr;
    oth_ptr.uch_ptr = other.getDataPtr();

    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] /= oth_ptr.f64_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] /= oth_ptr.f64_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] /= oth_ptr.f64_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] /= oth_ptr.f64_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] /= oth_ptr.f64_ptr[k]; }
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] /= oth_ptr.f32_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] /= oth_ptr.f32_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] /= oth_ptr.f32_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] /= oth_ptr.f32_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] /= oth_ptr.f32_ptr[k]; }
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] /= oth_ptr.int_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] /= oth_ptr.int_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] /= oth_ptr.int_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] /= oth_ptr.int_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] /= oth_ptr.int_ptr[k]; }
        }
    }else if (othrdatT== DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] /= oth_ptr.uch_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] /= oth_ptr.uch_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] /= oth_ptr.uch_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] /= oth_ptr.uch_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] /= oth_ptr.uch_ptr[k]; }
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : for(; k < length; ++k) { elptr.f64_ptr[k] /= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::FLOAT  : for(; k < length; ++k) { elptr.f32_ptr[k] /= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::INT    : for(; k < length; ++k) { elptr.int_ptr[k] /= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::UCHAR  : for(; k < length; ++k) { elptr.uch_ptr[k] /= oth_ptr.cmx_ptr[k]; } break;
        case DTYP::CMPLX  : for(; k < length; ++k) { elptr.cmx_ptr[k] /= oth_ptr.cmx_ptr[k]; }
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
        double *larg_mag_ch, tmp_mag;
        larg_mag_ch = new double[ch];
        for(; k<ch; ++k) { cl = elptr.cmx_ptr[k]; Aptrs.cmx_ptr[k]= cl; larg_mag_ch[k] = cl.square();}
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                tmp       = elptr.cmx_ptr[n];
                tmp_mag   = tmp.square();
                if( larg_mag_ch[k] < tmp_mag){
                    Aptrs.cmx_ptr[k] = tmp;
                    larg_mag_ch[k]      = tmp_mag;
                }
            }
        }
        delete [] larg_mag_ch;
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
        double *larg_mag_ch, tmp_mag;
        larg_mag_ch = new double[ch];
        for(; k<ch; ++k) { cl = elptr.cmx_ptr[k]; Aptrs.cmx_ptr[k]= cl; larg_mag_ch[k] = cl.square(); }
        for(m = ch ; m < length; m+=ch ){
            for(k=0, n=m ; k < ch; ++k, ++n){
                tmp       = elptr.cmx_ptr[n];
                tmp_mag   = tmp.square();
                if( larg_mag_ch[k] > tmp_mag){
                    Aptrs.cmx_ptr[k] = tmp;
                    larg_mag_ch[k]      = tmp_mag;
                }
            }
        }
        delete [] larg_mag_ch;
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

//Mat Mat::repeat_(const Mat& src, const uint32 rr, const uint32 rc, const uint32 rch){
//    if(src.isEmpty()) return Mat();
//
//    switch(src.getDatType()){
//    case DTYP::DOUBLE : return _repeat<double>(src, rr, rc, rch);
//    case DTYP::FLOAT  : return _repeat<float >(src, rr, rc, rch);
//    case DTYP::INT    : return _repeat<int32 >(src, rr, rc, rch);
//    case DTYP::UCHAR  : return _repeat<uchar >(src, rr, rc, rch);
//    case DTYP::CMPLX  : return _repeat<cmplx >(src, rr, rc, rch);
//    default           : return Mat();
//    }
//}
Mat Mat::repeat(const Mat& src, const uint32 rr, const uint32 rc, const uint32 rch){
    if(src.isEmpty()) return Mat();
    uint32 sr = src.getRow();
    uint32 sc = src.getCol();
    uint32 sch= src.getChannel();
    uint32 nr = sr * rr;
    uint32 nc = sc * rc;
    uint32 nch= sch*rch;

    uint32 x, y, z, i, k, n;
    uint32 ch_xpnd_sz = nch*sc*src.byteStep;
    uint32 cl_xpnd_sz = nch*nc*src.byteStep;
    Mat des(src.getDatType(), nr, nc, nch);
    elemptr sRow_ptr;
    elemptr tRow_ptr;
    elemptr ch_xpnd;

    for( y=0; y < sr ; ++y){
        sRow_ptr = src.getRowElptr(y);
        tRow_ptr = des.getRowElptr(y);
        ch_xpnd  = tRow_ptr;
        // make channel expanding array
        for(i=0; i < rch; ++i){
            if(src.getDatType() == DTYP::DOUBLE){
                for( x=0, z=i*sch, n=0; x < sc ; x+=sch, z+=nch){
                    for(k=0 ; k < sch; ++k, ++n)
                        ch_xpnd.f64_ptr[z+k] = sRow_ptr.f64_ptr[n];
                }
            }else if(src.getDatType()==DTYP::FLOAT){
                for( x=0, z=i*sch, n=0; x < sc ; x+=sch, z+=nch){
                    for(k=0 ; k < sch; ++k, ++n)
                        ch_xpnd.f32_ptr[z+k] = sRow_ptr.f32_ptr[n];
                }
            }else if(src.getDatType()==DTYP::INT){
                for( x=0, z=i*sch, n=0; x < sc ; x+=sch, z+=nch){
                    for(k=0 ; k < sch; ++k, ++n)
                        ch_xpnd.int_ptr[z+k] = sRow_ptr.int_ptr[n];
                }
            }else if(src.getDatType()==DTYP::UCHAR){
                for( x=0, z=i*sch, n=0; x < sc ; x+=sch, z+=nch){
                    for(k=0 ; k < sch; ++k, ++n)
                        ch_xpnd.uch_ptr[z+k] = sRow_ptr.uch_ptr[n];
                }
            }else if(src.getDatType()==DTYP::CMPLX){
                for( x=0, z=i*sch, n=0; x < sc ; x+=sch, z+=nch){
                    for(k=0 ; k < sch; ++k, ++n)
                        ch_xpnd.cmx_ptr[z+k] = sRow_ptr.cmx_ptr[n];
                }
            }
        }
        // replicating the expanded channel data into new column of des Mat.
        for(i=1, tRow_ptr.uch_ptr += ch_xpnd_sz; i < rc; ++i, tRow_ptr.uch_ptr += ch_xpnd_sz){
            memcpy(tRow_ptr.uch_ptr, ch_xpnd.uch_ptr, ch_xpnd_sz);
        }
    }
    // replicating the expanded column Mat into remaining rows of des Mat.
    for(y=0; y < sr ; ++y){
        sRow_ptr = des.getRowElptr(y);
        for( i=1, n=sr+y; i < rr; ++i, n+=sr){
            tRow_ptr = des.getRowElptr(n);
            memcpy(tRow_ptr.uch_ptr, sRow_ptr.uch_ptr, cl_xpnd_sz);
        }
    }
    return des;
}
} // end of jmat namespace
