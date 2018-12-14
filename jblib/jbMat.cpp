#include <stdio.h>
#include "jbMat.h"
#include <iostream>

namespace jmat {


int32 Mat::instant_count = 0;

//-- shallow copy version & using shared_ptr
void Mat::alloc(uint32 len){

    if(len < 0){
        mA = nullptr;
        dat_ptr = nullptr;
    }else{
        mA = shr_ptr (new uchar[len], std::default_delete<uchar[]>());
        dat_ptr = mA.get();
    }
    /*
    try{
        mA = (len==0)? nullptr : new double[len];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"memory allocation error: %s\n",ex.what());
        mA = nullptr;
    }
    */
}
void Mat::init(uint32 r, uint32 c, uint32 ch, DTYP dt){
    row = r;
    col = c;
    Nch = ch;
    lenRowCol = r*c;
    length = lenRowCol * ch;    
    datT = dt;
    switch(datT){
    case DTYP::UCHAR  : byteStep = 1; break;
    case DTYP::INT    : byteStep = 4; break;
    case DTYP::FLOAT  : byteStep = 4; break;
    case DTYP::DOUBLE : byteStep = 8; break;
    default           : byteStep = 1;
    }
    byteLen = length*byteStep;
    alloc(byteLen);
}
void Mat::initName(){
    obj_name =std::string( "Mat_") + std::to_string (instant_count);
    instant_count++;
}

Mat::Mat(DTYP dt):mA(nullptr){
    init(0,0,0, dt);
    initName();

}
Mat::Mat(DTYP dt, uint32 rc ):mA(nullptr){
    init(rc, rc, 1, dt);
    initName();
}
Mat::Mat(DTYP dt, uint32 r, uint32 c, uint32 ch):mA(nullptr){
    init(r, c, ch, dt);
    initName();
}

Mat::Mat(shr_ptr ma, DTYP dt, uint32 r, uint32 c, uint32 ch):mA(ma),row(r),col(c),Nch(ch),datT(dt){
    lenRowCol = r*c;
    length = lenRowCol * ch;
    switch(dt){
    case DTYP::UCHAR  : byteStep = 1; break;
    case DTYP::INT    : byteStep = 4; break;
    case DTYP::FLOAT  : byteStep = 4; break;
    case DTYP::DOUBLE : byteStep = 8; break;
    default           : byteStep = 1;
    }
    byteLen = length*byteStep;

    initName();
}

Mat::Mat(DTYP dt, uint32 r, uint32 c, uint32 ch, std::string name):mA(nullptr){
    init(r, c, ch, dt);
    obj_name = name;
}

Mat::Mat(const Mat& mat){
    row = mat.getRow();
    col = mat.getCol();
    Nch = mat.getChannel();
    length    = mat.getLength();
    lenRowCol = row*col;
    mA        = mat.getMat();
    datT      = mat.datT;
    byteStep  = mat.byteStep;
    byteLen   = mat.byteLen;
    sync_data_ptr();
    initName();

    /* -- deep copying
    alloc(length);

    double *pt_matdat = mat.getMat().get();
    double *pt_thisdat = mA.get();
    //for(int32 i=0; i < length; i++)
    //    pt_thisdat[i] = pt_matdat[i];
    std::copy(pt_matdat,pt_matdat+mat.length,pt_thisdat);
    */
    fprintf(stdout,"copy constructor\n");
}

Mat::Mat( std::initializer_list<double> list ){
    //-- Making a vector by column vector type    
    init(list.size(),1,1,DTYP::DOUBLE);

    double * dat_p = (double *)mA.get();
    if(mA!=nullptr){
        std::vector<double> v;
        v.insert(v.end(),list.begin(),list.end());

        for(int32 i=0;i<length;i++){
            dat_p[i] = v.at(i);
        }
    }
}
Mat::Mat( std::initializer_list<int32> list ){
    //-- Making a vector by column vector type
    init(list.size(),1,1,DTYP::INT);

    int32 *dat_p = (int32 *)mA.get();
    if(mA!=nullptr){
        std::vector<int32> v;
        v.insert(v.end(),list.begin(),list.end());

        for(int32 i=0;i<length;i++){
            dat_p[i] = v.at(i);
        }
    }
}
Mat::Mat( std::initializer_list<float> list ){
    //-- Making a vector by column vector type
    init(list.size(),1,1,DTYP::FLOAT);

    float *dat_p = (float *)mA.get();
    if(mA!=nullptr){
        std::vector<float> v;
        v.insert(v.end(),list.begin(),list.end());

        for(int32 i=0;i<length;i++){
            dat_p[i] = v.at(i);
        }
    }
}

Mat::~Mat(){}

void Mat::setRowCol(uint32 r, uint32 c, uint32 ch){
    int32 lenrc = r*c;
    int32 len   = lenrc*ch;

    if(length != len){
        mA.reset();
        alloc(len);
    }
    row = r;
    col = c;
    Nch = ch;
    lenRowCol = lenrc;
    length = len;
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

// call by reference
Mat& Mat::operator=(const Mat& other){

    fprintf(stdout,"Assign operator\n");
    row = other.getRow();
    col = other.getCol();
    Nch = other.getChannel();
    length = other.getLength();
    lenRowCol = row*col;
    byteStep  = getByteStep();
    byteLen   = length*byteStep;
    datT = getDatType();
    mA   = other.getMat();
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
    return multiplyMat(other);
}

Mat& Mat::operator*=(const double scalar){
    return minusScalar(scalar);
}

Mat& Mat::operator/=(const Mat& other){
    return divideMat(other);
}

Mat& Mat::operator/=(const double scalar){
    return divideScalar(scalar);
}


Mat Mat::operator+(const Mat& other) const{
    if( row != other.getRow() || col != other.getCol()|| Nch != other.getChannel() ){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return Mat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return Mat();
    }

    Mat sum = this->copy();
    sum.plusMat(other);

    return sum;
}

Mat Mat::operator+(const double scalar) const{
    if(isEmpty()){
        fprintf(stderr," operator+ with scalar : this mat is empty \n");
        return Mat();
    }
    Mat sum = this->copy();
    sum.plusScalar(scalar);

    return sum;
}

Mat Mat::operator-(const Mat& other) const{
    if( row != other.getRow() || col != other.getCol() || Nch != other.getChannel() ){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return Mat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return Mat();
    }

    Mat subtract = this->copy();
    subtract.minusMat(other);

    return subtract;
}
Mat Mat::operator-(const double scalar) const{
    if(isEmpty()){
        fprintf(stderr," operator- with scalar : this mat is empty \n");
        return Mat();
    }
    Mat subtract = this->copy();
    subtract.minusScalar(scalar);

    return subtract;
}

Mat Mat::operator*(const Mat& other) const{
    if( row != other.getRow() || col != other.getCol() || Nch != other.getChannel() ){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return Mat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return Mat();
    }
    Mat product = this->copy();
    product.multiplyMat(other);

    return product;
}
Mat Mat::operator*(const double scalar) const{
    if(isEmpty()){
        fprintf(stderr," operator* with scalar : this mat is empty \n");
        return Mat();
    }
    Mat product = this->copy();
    product.multiplyScalar(scalar);

    return product;
}

Mat Mat::operator/(const Mat& other) const{
    if( row != other.getRow() || col != other.getCol() || Nch != other.getChannel() ){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return Mat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return Mat();
    }
    Mat div = this->copy();
    div.divideMat(other);

    return div;
}
Mat Mat::operator/(const double scalar) const{
    if(isEmpty()){
        fprintf(stderr," operator/ with scalar : this mat is empty \n");
        return Mat();
    }
    Mat div = this->copy();
    div.divideScalar(scalar);

    return div;
}


uint32 Mat::reshape(uint32 r, uint32 c, uint32 ch){
    int32 rc   = r*c;
    int32 tlen = rc*ch;
    if( tlen != length){
        fprintf(stderr," reshape argument is not correct!\n");
        return -1;
    }
    row = r;
    col = c;
    Nch = ch;
    lenRowCol = rc;

    return 0;
}

void Mat::changeDType(const DTYP dt){
    if(dt == DTYP::DOUBLE){
        switch (datT) {
        case DTYP::FLOAT : _type_change<float,double>(); break;
        case DTYP::INT   : _type_change<int32  ,double>(); break;
        case DTYP::UCHAR : _type_change<uchar,double>(); break;
        default          : break;
        }
        datT = dt;
    }else if( dt== DTYP::FLOAT){
        switch (datT) {
        case DTYP::DOUBLE: _type_change<double,float>(); break;
        case DTYP::INT   : _type_change<int32   ,float>(); break;
        case DTYP::UCHAR : _type_change<uchar ,float>(); break;
        default          : break;
        }
        datT = dt;
    }else if( dt== DTYP::INT){
        switch (datT) {
        case DTYP::DOUBLE: _type_change<double,int32>(); break;
        case DTYP::FLOAT : _type_change<float ,int32>(); break;
        case DTYP::UCHAR : _type_change<uchar ,int32>(); break;
        default          : break;
        }
        datT = dt;
    }else if( dt== DTYP::UCHAR){
        switch (datT) {
        case DTYP::DOUBLE: _type_change<double,uchar>(); break;
        case DTYP::FLOAT : _type_change<float ,uchar>(); break;
        case DTYP::INT   : _type_change<int32   ,uchar>(); break;
        default          : break;
        }
        datT = dt;
    }else{
        fprintf(stderr, "argument dt of changeDType() is not match with any type\n");
    }
#ifdef _DEBUG
    switch(datT){
    case DTYP::DOUBLE: fprintf(stdout,"datT is changed to double\n"); break;
    case DTYP::FLOAT : fprintf(stdout,"datT is changed to float \n"); break;
    case DTYP::INT   : fprintf(stdout,"datT is changed to int32   \n"); break;
    case DTYP::UCHAR : fprintf(stdout,"datT is changed to uchar \n"); break;
    }
#endif
}

void Mat::transpose(){
    if(isEmpty()){
        fprintf(stderr," Transpose: This Mat is empty\n");
        return ;
    }

    uchar *tmA;
    uchar *mdat;

    try{
        tmA= new uchar[static_cast<unsigned long>(byteLen)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"Transpose Error: %s\n",ex.what());
        return;
    }

    mdat = mA.get();
    uint32 i, j, k, m, ch_offset, lhs_idx, rhs_idx;
    for(k=0; k < Nch; k++){       // channel step
        ch_offset = k* lenRowCol * byteStep;
        for(i=0; i < row; i++){   // row step
            for(j=0; j<col; j++){ // column step
                lhs_idx = (ch_offset + j*row +i)*byteStep;
                rhs_idx = (ch_offset + i*col +j)*byteStep;
                for(m=0; m< byteStep; m++){ // byte step
                    tmA[lhs_idx] = mdat[rhs_idx];
                }
            }
        }
    }
    mA.reset(tmA,std::default_delete<uchar[]>());
    sync_data_ptr();

    int32 row_tr = col;
    col = row;
    row = row_tr;
}

void Mat::printMat(const std::string objname) {

    if(!objname.empty())
        fprintf(stdout,"object : %s \n", objname.c_str());
    if(isEmpty()){
        fprintf(stdout,"Empty matrix\n");
        return ;
    }

    switch(datT){
    case DTYP::DOUBLE : _print((double *)dat_ptr); break;
    case DTYP::FLOAT  : _print((float  *)dat_ptr); break;
    case DTYP::INT    : _print((int32    *)dat_ptr); break;
    case DTYP::UCHAR  : _print((uchar  *)dat_ptr); break;
    default           : fprintf(stdout, "datT is not matching with any type\n");
    }
}

void Mat::printMat()  {
    printMat(obj_name);
}

Mat Mat::copy() const{
    Mat A(this->datT, row, col, Nch);

    uchar *pt_matdat  = A.getMat().get();
    uchar *pt_thisdat = this->mA.get();

    std::copy(pt_thisdat,pt_thisdat+length*byteStep,pt_matdat);

    return A;
}

Mat Mat::ones(uint32 r, uint32 c, uint32 ch, DTYP dt){
    if( r <= 0 || c <= 0 || ch <= 0){
        fprintf(stdout,"In ones method: arguments r , c and ch are to be larger than 0 ");
        return Mat();
    }

    Mat A(dt, r, c, ch);

    if(dt==DTYP::DOUBLE){
        double* pt_dat = (double *)A.getMat().get();
        for(uint32 i=0; i < r*c*ch; i++)
            pt_dat[i] = 1.0;
    }else if(dt==DTYP::FLOAT){
        float* pt_dat = (float *)A.getMat().get();
        for(uint32 i=0; i < r*c*ch; i++)
            pt_dat[i] = 1.0f;
    }else if(dt==DTYP::INT){
        int32* pt_dat = (int32 *)A.getMat().get();
        for(uint32 i=0; i < r*c*ch; i++)
            pt_dat[i] = 1;
    }else{
        uchar* pt_dat = (uchar *)A.getMat().get();
        for(uint32 i=0; i < r*c*ch; i++)
            pt_dat[i] = 1.0;
    }

    return A;
}


Mat Mat::zeros(uint32 r, uint32 c, uint32 ch, DTYP dt){
    if( r <= 0 || c <= 0 || ch <= 0){
        fprintf(stdout,"In zeros method: arguments r , c and ch are to be larger than 0 ");
        return Mat();
    }

    Mat A(dt, r, c, ch);

    uchar* pt_dat = A.getMat().get();
    for(uint32 i=0; i < r*c*ch*A.getByteStep(); i++){
        pt_dat[i] = 0;
    }

    return A;
}

Mat Mat::copyChannelN(const uint32 NoCh) const{
    Mat A(this->datT, row,col,1);

    int32 numCh = NoCh;
    if(numCh >= A.getChannel()){
        fprintf(stdout,"channel number of getNchannel() is out of bound\n The last channel is selected\n");
        numCh = A.getChannel()-1;
    }

    uchar *srcDat_pt = mA.get();
    uchar *tarDat_pt = A.mA.get();
    uint32 offset = numCh*lenRowCol*byteStep;
    for(uint32 i=offset; i < lenRowCol*byteStep; i++ )
        tarDat_pt[i] = srcDat_pt[i];

    return A;
}

void Mat::setChannelN(const Mat& src, const uint32 srcFromCh,const uint32 Channels, const uint32 tarToCh){
    if(src.getChannel() < srcFromCh+Channels ){
        fprintf(stdout,"setChannelN(): srcFromCh and Channels are not correct! \n");
        return ;
    }

    uint32 tar_ch;
    if(isEmpty()){
        init(src.getRow(),src.getCol(), Channels, src.getDatType() );
        tar_ch = 0;
    }else if( Nch < tarToCh+Channels){
        fprintf(stdout,"setChannelN(): tarToCh and Channels are not correct! \n");
        return ;
    }else if( src.getRow() != row && src.getCol() != col){
        fprintf(stdout,"*this and src are not equal size. *this:(%d, %d, %d) != src(%d, %d, %d)\n",row,col,Nch,src.getRow(),src.getCol(),src.getChannel());
        return ;
    }else
        tar_ch = tarToCh;

    sync_data_ptr();
    uchar *srcdat_ptr = src.getMat().get();
    uint32 src_idx = srcFromCh*lenRowCol*byteStep;
    uint32 tar_idx = tar_ch*lenRowCol*byteStep;
    uint32 chStep  = lenRowCol*byteStep;

    for(uint32 j=0; j < Channels ; j++){
        for(uint32 i=0; i < lenRowCol*byteStep; i++ ){
            dat_ptr[tar_idx++] = srcdat_ptr[src_idx++];
        }
    }
}

void Mat::setName(std::string name){
    obj_name = name;
}

Mat Mat::copySubMat(const uint32 startRow, const uint32 endRow, const uint32 startCol, const uint32 endCol) const {
    if( startRow > row || endRow > row || startCol > col || endCol > col){
        fprintf(stdout,"copySubMat() : one or more arguments are out of bound from *this mat \n");
        return Mat();
    }

    int32 new_row = endRow-startRow+1;
    int32 new_col = endCol-startCol+1;
    Mat A( datT, new_row, new_col, Nch);
    uchar *tardat_ptr = A.getMat().get();

    uint32 ch_offset = 0 ;
    uint32 k = 0;
    uint32 r, c, ch;
    uint32 offset;
    uint32 lenRCByteStep = lenRowCol*byteStep;
    uint32 rowByteStep   = col*byteStep;
    uint32 startColByte  = startCol*byteStep;
    uint32 endColByte    = endCol*byteStep;
    uint32 colstart, colend;
    uint32 rowstart      = startRow*rowByteStep;
    uint32 rowend        = endRow*rowByteStep;
    for( ch=0; ch < Nch; ch++){
        for( r = rowstart; r <= rowend; r+=rowByteStep ){
            offset   = ch_offset + r ;
            colstart = offset + startColByte;
            colend   = offset + endColByte;
            for( c = colstart; c <= colend ; c++){
                tardat_ptr[k++] = dat_ptr[ c ];
            }
        }
        ch_offset += lenRCByteStep;
    }
    return A;
}

Mat& Mat::plusMat(const Mat& other){
    if(isEmpty() || other.isEmpty()) {
        fprintf(stdout, "plusMat method : either of this or other is empty\n");
        return *this;
    }

    DTYP othrdatT = other.getDatType();
    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : _plus_mat<double>( getDataPtr<double>(), other.getDataPtr<double>(), length); break;
        case DTYP::FLOAT  : _plus_mat<float >( getDataPtr<float >(), other.getDataPtr<double>(), length); break;
        case DTYP::INT    : _plus_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<double>(), length); break;
        case DTYP::UCHAR  : _plus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        default           : _plus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : _plus_mat<double>( getDataPtr<double>(), other.getDataPtr<float >(), length); break;
        case DTYP::FLOAT  : _plus_mat<float >( getDataPtr<float >(), other.getDataPtr<float >(), length); break;
        case DTYP::INT    : _plus_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<float >(), length); break;
        case DTYP::UCHAR  : _plus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        default           : _plus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : _plus_mat<double>( getDataPtr<double>(), other.getDataPtr<int32   >(), length); break;
        case DTYP::FLOAT  : _plus_mat<float >( getDataPtr<float >(), other.getDataPtr<int32   >(), length); break;
        case DTYP::INT    : _plus_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<int32   >(), length); break;
        case DTYP::UCHAR  : _plus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<int32   >(), length); break;
        default           : _plus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<int32   >(), length); break;
        }
    }else{
        switch ( datT ){
        case DTYP::DOUBLE : _plus_mat<double>( getDataPtr<double>(), other.getDataPtr<uchar >(), length); break;
        case DTYP::FLOAT  : _plus_mat<float >( getDataPtr<float >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::INT    : _plus_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::UCHAR  : _plus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        default           : _plus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        }
    }
    return *this;
}

Mat& Mat::minusMat(const Mat& other){
    if(isEmpty() || other.isEmpty()) {
        fprintf(stdout, "plusMat method : either of this or other is empty\n");
        return *this;
    }
    DTYP othrdatT = other.getDatType();
    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : _minus_mat<double>( getDataPtr<double>(), other.getDataPtr<double>(), length); break;
        case DTYP::FLOAT  : _minus_mat<float >( getDataPtr<float >(), other.getDataPtr<double>(), length); break;
        case DTYP::INT    : _minus_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<double>(), length); break;
        case DTYP::UCHAR  : _minus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        default           : _minus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : _minus_mat<double>( getDataPtr<double>(), other.getDataPtr<float >(), length); break;
        case DTYP::FLOAT  : _minus_mat<float >( getDataPtr<float >(), other.getDataPtr<float >(), length); break;
        case DTYP::INT    : _minus_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<float >(), length); break;
        case DTYP::UCHAR  : _minus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        default           : _minus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : _minus_mat<double>( getDataPtr<double>(), other.getDataPtr<int32   >(), length); break;
        case DTYP::FLOAT  : _minus_mat<float >( getDataPtr<float >(), other.getDataPtr<int32   >(), length); break;
        case DTYP::INT    : _minus_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<int32   >(), length); break;
        case DTYP::UCHAR  : _minus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<int32   >(), length); break;
        default           : _minus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<int32   >(), length); break;
        }
    }else{
        switch ( datT ){
        case DTYP::DOUBLE : _minus_mat<double>( getDataPtr<double>(), other.getDataPtr<uchar >(), length); break;
        case DTYP::FLOAT  : _minus_mat<float >( getDataPtr<float >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::INT    : _minus_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::UCHAR  : _minus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        default           : _minus_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        }
    }
    return *this;
}

Mat& Mat::multiplyMat(const Mat& other){
    if(isEmpty() || other.isEmpty()) {
        fprintf(stdout, "plusMat method : either of this or other is empty\n");
        return *this;
    }
    DTYP othrdatT = other.getDatType();
    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : _multiply_mat<double>( getDataPtr<double>(), other.getDataPtr<double>(), length); break;
        case DTYP::FLOAT  : _multiply_mat<float >( getDataPtr<float >(), other.getDataPtr<double>(), length); break;
        case DTYP::INT    : _multiply_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<double>(), length); break;
        case DTYP::UCHAR  : _multiply_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        default           : _multiply_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : _multiply_mat<double>( getDataPtr<double>(), other.getDataPtr<float >(), length); break;
        case DTYP::FLOAT  : _multiply_mat<float >( getDataPtr<float >(), other.getDataPtr<float >(), length); break;
        case DTYP::INT    : _multiply_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<float >(), length); break;
        case DTYP::UCHAR  : _multiply_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        default           : _multiply_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : _multiply_mat<double>( getDataPtr<double>(), other.getDataPtr<int32   >(), length); break;
        case DTYP::FLOAT  : _multiply_mat<float >( getDataPtr<float >(), other.getDataPtr<int32   >(), length); break;
        case DTYP::INT    : _multiply_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<int32   >(), length); break;
        case DTYP::UCHAR  : _multiply_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<int32   >(), length); break;
        default           : _multiply_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<int32   >(), length); break;
        }
    }else{
        switch ( datT ){
        case DTYP::DOUBLE : _multiply_mat<double>( getDataPtr<double>(), other.getDataPtr<uchar >(), length); break;
        case DTYP::FLOAT  : _multiply_mat<float >( getDataPtr<float >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::INT    : _multiply_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::UCHAR  : _multiply_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        default           : _multiply_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        }
    }
    return *this;
}

Mat& Mat::divideMat(const Mat& other){
    if(isEmpty() || other.isEmpty()) {
        fprintf(stdout, "plusMat method : either of this or other is empty\n");
        return *this;
    }
    DTYP othrdatT = other.getDatType();
    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : _divide_mat<double>( getDataPtr<double>(), other.getDataPtr<double>(), length); break;
        case DTYP::FLOAT  : _divide_mat<float >( getDataPtr<float >(), other.getDataPtr<double>(), length); break;
        case DTYP::INT    : _divide_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<double>(), length); break;
        case DTYP::UCHAR  : _divide_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        default           : _divide_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : _divide_mat<double>( getDataPtr<double>(), other.getDataPtr<float >(), length); break;
        case DTYP::FLOAT  : _divide_mat<float >( getDataPtr<float >(), other.getDataPtr<float >(), length); break;
        case DTYP::INT    : _divide_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<float >(), length); break;
        case DTYP::UCHAR  : _divide_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        default           : _divide_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : _divide_mat<double>( getDataPtr<double>(), other.getDataPtr<int32   >(), length); break;
        case DTYP::FLOAT  : _divide_mat<float >( getDataPtr<float >(), other.getDataPtr<int32   >(), length); break;
        case DTYP::INT    : _divide_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<int32   >(), length); break;
        case DTYP::UCHAR  : _divide_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<int32   >(), length); break;
        default           : _divide_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<int32   >(), length); break;
        }
    }else{
        switch ( datT ){
        case DTYP::DOUBLE : _divide_mat<double>( getDataPtr<double>(), other.getDataPtr<uchar >(), length); break;
        case DTYP::FLOAT  : _divide_mat<float >( getDataPtr<float >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::INT    : _divide_mat<int32   >( getDataPtr<int32   >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::UCHAR  : _divide_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        default           : _divide_mat<uchar >( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        }
    }
    return *this;
}

Mat& Mat::plusScalar(const double scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _plus_scalar<double,double>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _plus_scalar<float ,double>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _plus_scalar<int32   ,double>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _plus_scalar<uchar ,double>(getDataPtr<uchar >(), scalar, length); break;
    default           : _plus_scalar<uchar ,double>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::plusScalar(const float scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _plus_scalar<double,float>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _plus_scalar<float ,float>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _plus_scalar<int32   ,float>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _plus_scalar<uchar ,float>(getDataPtr<uchar >(), scalar, length); break;
    default           : _plus_scalar<uchar ,float>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::plusScalar(const int32 scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _plus_scalar<double,int32>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _plus_scalar<float ,int32>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _plus_scalar<int32   ,int32>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _plus_scalar<uchar ,int32>(getDataPtr<uchar >(), scalar, length); break;
    default           : _plus_scalar<uchar ,int32>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}
Mat& Mat::plusScalar(const uchar scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _plus_scalar<double,uchar>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _plus_scalar<float ,uchar>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _plus_scalar<int32   ,uchar>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _plus_scalar<uchar ,uchar>(getDataPtr<uchar >(), scalar, length); break;
    default           : _plus_scalar<uchar ,uchar>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::minusScalar(const double scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _minus_scalar<double,double>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _minus_scalar<float ,double>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _minus_scalar<int32   ,double>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _minus_scalar<uchar ,double>(getDataPtr<uchar >(), scalar, length); break;
    default           : _minus_scalar<uchar ,double>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::minusScalar(const float scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _minus_scalar<double,float>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _minus_scalar<float ,float>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _minus_scalar<int32   ,float>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _minus_scalar<uchar ,float>(getDataPtr<uchar >(), scalar, length); break;
    default           : _minus_scalar<uchar ,float>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::minusScalar(const int32 scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _minus_scalar<double,int32>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _minus_scalar<float ,int32>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _minus_scalar<int32   ,int32>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _minus_scalar<uchar ,int32>(getDataPtr<uchar >(), scalar, length); break;
    default           : _minus_scalar<uchar ,int32>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}
Mat& Mat::minusScalar(const uchar scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _minus_scalar<double,uchar>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _minus_scalar<float ,uchar>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _minus_scalar<int32   ,uchar>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _minus_scalar<uchar ,uchar>(getDataPtr<uchar >(), scalar, length); break;
    default           : _minus_scalar<uchar ,uchar>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::multiplyScalar(const double scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _multiply_scalar<double,double>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _multiply_scalar<float ,double>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _multiply_scalar<int32   ,double>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _multiply_scalar<uchar ,double>(getDataPtr<uchar >(), scalar, length); break;
    default           : _multiply_scalar<uchar ,double>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::multiplyScalar(const float scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _multiply_scalar<double,float>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _multiply_scalar<float ,float>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _multiply_scalar<int32   ,float>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _multiply_scalar<uchar ,float>(getDataPtr<uchar >(), scalar, length); break;
    default           : _multiply_scalar<uchar ,float>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::multiplyScalar(const int32 scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _multiply_scalar<double,int32>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _multiply_scalar<float ,int32>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _multiply_scalar<int32   ,int32>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _multiply_scalar<uchar ,int32>(getDataPtr<uchar >(), scalar, length); break;
    default           : _multiply_scalar<uchar ,int32>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}
Mat& Mat::multiplyScalar(const uchar scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _multiply_scalar<double,uchar>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _multiply_scalar<float ,uchar>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _multiply_scalar<int32   ,uchar>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _multiply_scalar<uchar ,uchar>(getDataPtr<uchar >(), scalar, length); break;
    default           : _multiply_scalar<uchar ,uchar>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::divideScalar(const double scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _divide_scalar<double,double>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _divide_scalar<float ,double>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _divide_scalar<int32   ,double>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _divide_scalar<uchar ,double>(getDataPtr<uchar >(), scalar, length); break;
    default           : _divide_scalar<uchar ,double>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::divideScalar(const float scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _divide_scalar<double,float>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _divide_scalar<float ,float>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _divide_scalar<int32   ,float>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _divide_scalar<uchar ,float>(getDataPtr<uchar >(), scalar, length); break;
    default           : _divide_scalar<uchar ,float>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

Mat& Mat::divideScalar(const int32 scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _divide_scalar<double,int32>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _divide_scalar<float ,int32>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _divide_scalar<int32   ,int32>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _divide_scalar<uchar ,int32>(getDataPtr<uchar >(), scalar, length); break;
    default           : _divide_scalar<uchar ,int32>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}
Mat& Mat::divideScalar(const uchar scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _divide_scalar<double,uchar>(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _divide_scalar<float ,uchar>(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _divide_scalar<int32   ,uchar>(getDataPtr<int32   >(), scalar, length); break;
    case DTYP::UCHAR  : _divide_scalar<uchar ,uchar>(getDataPtr<uchar >(), scalar, length); break;
    default           : _divide_scalar<uchar ,uchar>(getDataPtr<uchar >(), scalar, length);
    }
    return *this;
}

} // end of jmat namespace
