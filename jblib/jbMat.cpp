#include <stdio.h>
#include "jbMat.h"
#include <iostream>
#include <utility>

namespace jmat {


int32 Mat::instant_count = 0;

//-- shallow copy version & using shared_ptr
void Mat::alloc(uint32 len){

    if(len == 0){
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
    case DTYP::CMPLX  : byteStep = sizeof(cmplx); break;
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
    row = mat.row;
    col = mat.col;
    Nch = mat.Nch;
    length    = mat.length;
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
    //for(int32 i=0; i < length; i++)
    //    pt_thisdat[i] = pt_matdat[i];
    std::copy(pt_matdat,pt_matdat+mat.length,pt_thisdat);
    */
#ifdef _DEBUG_MODE_
    fprintf(stdout,"copy constructor\n");
#endif
}

Mat::Mat( Mat&& mat){
    row = mat.row;
    col = mat.col;
    Nch = mat.Nch;
    length    = mat.length;
    lenRowCol = mat.lenRowCol;
    mA        = mat.mA;
    datT      = mat.datT;
    byteStep  = mat.byteStep;
    byteLen   = mat.byteLen;
    sync_data_ptr();
    initName();

#ifdef _DEBUG_MODE_
    fprintf(stdout,"move constructor\n");
#endif
}

Mat::Mat( std::initializer_list<double> list ){
    //-- Making a vector by column vector type    
    init(list.size(),1,1,DTYP::DOUBLE);

    double * dat_p = (double *)mA.get();
    if(mA!=nullptr){
        std::vector<double> v;
        v.insert(v.end(),list.begin(),list.end());

        for(uint32 i=0;i<length;i++){
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

        for(uint32 i=0;i<length;i++){
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

        for(uint32 i=0;i<length;i++){
            dat_p[i] = v.at(i);
        }
    }
}

Mat::~Mat(){}

void Mat::setRowCol(uint32 r, uint32 c, uint32 ch){
    uint32 lenrc = r*c;
    uint32 len   = lenrc*ch;

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

// for assigning lvalue of right-hand side
Mat& Mat::operator=(const Mat& other){
//The parameter to the ‘operator=()’ is passed by value which calls copy constructor
//to create an object local to the ‘operator=()’.
//Than the value of the temp object is swapped with ‘*this’ object
#ifdef _DEBUG_MODE_
    fprintf(stdout,"rvalue referance Assign operator\n");
#endif
    if( &other == this) return *this;
    row      = other.row;
    col      = other.col;
    Nch      = other.Nch;
    length   = other.length;
    lenRowCol= other.lenRowCol;
    byteStep = other.byteStep;
    byteLen  = other.byteLen;
    datT     = other.datT;
    mA       = other.mA;
    sync_data_ptr();
    return *this;
}


// for rvalue assignment of right-hand side
Mat& Mat::operator=(Mat&& other){

#ifdef _DEBUG_MODE_
    fprintf(stdout,"rvalue referance Assign operator\n");
#endif
    if( &other == this) return *this;
    std::swap(row,other.row);
    std::swap(col,other.col);
    std::swap(Nch,other.Nch);
    std::swap(length,other.length);
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
    return multiplyMat(other);
}

Mat& Mat::operator*=(const double scalar){
    return multiplyScalar(scalar);
}

Mat& Mat::operator/=(const Mat& other){
    return divideMat(other);
}

Mat& Mat::operator/=(const double scalar){
    return dividedByScalar(scalar);
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

Mat operator+(const Mat& A, const double scalar){
    if(A.isEmpty()){
        fprintf(stderr," operator+ with scalar : Mat A is empty \n");
        return Mat();
    }
    Mat sum = A.copy();
    sum.plusScalar(scalar);

    return sum;
}
Mat operator+(const double scalar , const Mat& A){
    if(A.isEmpty()){
        fprintf(stderr," operator+ with scalar : Mat A is empty \n");
        return Mat();
    }
    Mat sum = A.copy();
    sum.plusScalar(scalar);

    return sum;
}

Mat operator-(const Mat& A, const double scalar) {
    if(A.isEmpty()){
        fprintf(stderr," operator- with scalar : this mat is empty \n");
        return Mat();
    }
    Mat subtract = A.copy();
    subtract.minusScalar(scalar);

    return subtract;
}
Mat operator-(const double scalar,const Mat& A) {
    if(A.isEmpty()){
        fprintf(stderr," operator- with scalar : this mat is empty \n");
        return Mat();
    }
    Mat subtract = A.copy();
    subtract.minusScalar(scalar);

    return subtract;
}

Mat operator*(const Mat& A, const double scalar){
    if(A.isEmpty()){
        fprintf(stderr," operator* with scalar : this mat is empty \n");
        return Mat();
    }
    Mat product = A.copy();
    product.multiplyScalar(scalar);

    return product;
}
Mat operator*(const double scalar, const Mat& A){
    if(A.isEmpty()){
        fprintf(stderr," operator* with scalar : this mat is empty \n");
        return Mat();
    }
    Mat product = A.copy();
    product.multiplyScalar(scalar);

    return product;
}

Mat operator/(const Mat& lhs, const double scalar){
    if(lhs.isEmpty()){
        fprintf(stderr," operator/ with scalar : this mat is empty \n");
        return Mat();
    }
    Mat div = lhs.copy();
    div.dividedByScalar(scalar);

    return div;
}
Mat operator/(const double scalar, const Mat& rhs){
    if(rhs.isEmpty()){
        fprintf(stderr," operator/ with scalar : this mat is empty \n");
        return Mat();
    }
    Mat div = rhs.copy();
    div.dividingScalar(scalar);

    return div;
}

int32 Mat::reshape(uint32 r, uint32 c, uint32 ch){
    uint32 rc   = r*c;
    uint32 tlen = rc*ch;
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
        case DTYP::INT   : _type_change<int32,double>(); break;
        case DTYP::UCHAR : _type_change<uchar,double>(); break;
        default          : break;
        }
        datT = dt;
    }else if( dt== DTYP::FLOAT){
        switch (datT) {
        case DTYP::DOUBLE: _type_change<double,float>(); break;
        case DTYP::INT   : _type_change<int32 ,float>(); break;
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
        case DTYP::INT   : _type_change<int32 ,uchar>(); break;
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
                memcpy(&tmA[lhs_idx], &mdat[rhs_idx], byteStep);
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
    case DTYP::INT    : _print((int32  *)dat_ptr); break;
    case DTYP::UCHAR  : _print((uchar  *)dat_ptr); break;
    case DTYP::CMPLX  : _print((cmplx  *)dat_ptr); break;
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
    uint32 len = A.getLength();
    if(dt==DTYP::DOUBLE){
        double* pt_dat = A.getDataPtr<double>();
        for(uint32 i=0; i < len; i++)
            pt_dat[i] = 1.0;
    }else if(dt==DTYP::FLOAT){
        float* pt_dat = A.getDataPtr<float >();
        for(uint32 i=0; i < len; i++)
            pt_dat[i] = 1.0f;
    }else if(dt==DTYP::INT){
        int32* pt_dat = A.getDataPtr<int32 >();
        for(uint32 i=0; i < len; i++)
            pt_dat[i] = 1;
    }else if(dt==DTYP::UCHAR){
        uchar* pt_dat = A.getDataPtr<uchar >();
        for(uint32 i=0; i < len; i++)
            pt_dat[i] = 1;
    }else if(dt==DTYP::CMPLX){
        cmplx* pt_dat = A.getDataPtr<cmplx >();
        for(uint32 i=0; i < len; i++)
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
        fprintf(stdout, "plusMat method : either of this or other is empty\n");
        return *this;
    }

    DTYP othrdatT = other.getDatType();
    if(  othrdatT == DTYP::DOUBLE){
        switch ( datT ){
        case DTYP::DOUBLE : _plus_mat( getDataPtr<double>(), other.getDataPtr<double>(), length); break;
        case DTYP::FLOAT  : _plus_mat( getDataPtr<float >(), other.getDataPtr<double>(), length); break;
        case DTYP::INT    : _plus_mat( getDataPtr<int32 >(), other.getDataPtr<double>(), length); break;
        case DTYP::UCHAR  : _plus_mat( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        case DTYP::CMPLX  : _plus_mat( getDataPtr<cmplx >(), other.getDataPtr<double>(), length); break;
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : _plus_mat( getDataPtr<double>(), other.getDataPtr<float >(), length); break;
        case DTYP::FLOAT  : _plus_mat( getDataPtr<float >(), other.getDataPtr<float >(), length); break;
        case DTYP::INT    : _plus_mat( getDataPtr<int32 >(), other.getDataPtr<float >(), length); break;
        case DTYP::UCHAR  : _plus_mat( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        case DTYP::CMPLX  : _plus_mat( getDataPtr<cmplx >(), other.getDataPtr<float >(), length); break;
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : _plus_mat( getDataPtr<double>(), other.getDataPtr<int32 >(), length); break;
        case DTYP::FLOAT  : _plus_mat( getDataPtr<float >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::INT    : _plus_mat( getDataPtr<int32 >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::UCHAR  : _plus_mat( getDataPtr<uchar >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::CMPLX  : _plus_mat( getDataPtr<cmplx >(), other.getDataPtr<int32 >(), length); break;
        }
    }else if (othrdatT== DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : _plus_mat( getDataPtr<double>(), other.getDataPtr<uchar >(), length); break;
        case DTYP::FLOAT  : _plus_mat( getDataPtr<float >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::INT    : _plus_mat( getDataPtr<int32 >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::UCHAR  : _plus_mat( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::CMPLX  : _plus_mat( getDataPtr<cmplx >(), other.getDataPtr<uchar >(), length); break;
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : _plus_mat( getDataPtr<double>(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::FLOAT  : _plus_mat( getDataPtr<float >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::INT    : _plus_mat( getDataPtr<int32 >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::UCHAR  : _plus_mat( getDataPtr<uchar >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::CMPLX  : _plus_mat( getDataPtr<cmplx >(), other.getDataPtr<cmplx >(), length); break;
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
        case DTYP::DOUBLE : _minus_mat( getDataPtr<double>(), other.getDataPtr<double>(), length); break;
        case DTYP::FLOAT  : _minus_mat( getDataPtr<float >(), other.getDataPtr<double>(), length); break;
        case DTYP::INT    : _minus_mat( getDataPtr<int32 >(), other.getDataPtr<double>(), length); break;
        case DTYP::UCHAR  : _minus_mat( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        case DTYP::CMPLX  : _minus_mat( getDataPtr<cmplx >(), other.getDataPtr<double>(), length); break;
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : _minus_mat( getDataPtr<double>(), other.getDataPtr<float >(), length); break;
        case DTYP::FLOAT  : _minus_mat( getDataPtr<float >(), other.getDataPtr<float >(), length); break;
        case DTYP::INT    : _minus_mat( getDataPtr<int32 >(), other.getDataPtr<float >(), length); break;
        case DTYP::UCHAR  : _minus_mat( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        case DTYP::CMPLX  : _minus_mat( getDataPtr<cmplx >(), other.getDataPtr<float >(), length); break;
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : _minus_mat( getDataPtr<double>(), other.getDataPtr<int32 >(), length); break;
        case DTYP::FLOAT  : _minus_mat( getDataPtr<float >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::INT    : _minus_mat( getDataPtr<int32 >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::UCHAR  : _minus_mat( getDataPtr<uchar >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::CMPLX  : _minus_mat( getDataPtr<cmplx >(), other.getDataPtr<int32 >(), length); break;
        }
    }else if( othrdatT == DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : _minus_mat( getDataPtr<double>(), other.getDataPtr<uchar >(), length); break;
        case DTYP::FLOAT  : _minus_mat( getDataPtr<float >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::INT    : _minus_mat( getDataPtr<int32 >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::UCHAR  : _minus_mat( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::CMPLX  : _minus_mat( getDataPtr<cmplx >(), other.getDataPtr<uchar >(), length); break;
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : _minus_mat( getDataPtr<double>(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::FLOAT  : _minus_mat( getDataPtr<float >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::INT    : _minus_mat( getDataPtr<int32 >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::UCHAR  : _minus_mat( getDataPtr<uchar >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::CMPLX  : _minus_mat( getDataPtr<cmplx >(), other.getDataPtr<cmplx >(), length); break;
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
        case DTYP::DOUBLE : _multiply_mat( getDataPtr<double>(), other.getDataPtr<double>(), length); break;
        case DTYP::FLOAT  : _multiply_mat( getDataPtr<float >(), other.getDataPtr<double>(), length); break;
        case DTYP::INT    : _multiply_mat( getDataPtr<int32 >(), other.getDataPtr<double>(), length); break;
        case DTYP::UCHAR  : _multiply_mat( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        case DTYP::CMPLX  : _multiply_mat( getDataPtr<cmplx >(), other.getDataPtr<double>(), length); break;
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : _multiply_mat( getDataPtr<double>(), other.getDataPtr<float >(), length); break;
        case DTYP::FLOAT  : _multiply_mat( getDataPtr<float >(), other.getDataPtr<float >(), length); break;
        case DTYP::INT    : _multiply_mat( getDataPtr<int32 >(), other.getDataPtr<float >(), length); break;
        case DTYP::UCHAR  : _multiply_mat( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        case DTYP::CMPLX  : _multiply_mat( getDataPtr<cmplx >(), other.getDataPtr<float >(), length); break;
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : _multiply_mat( getDataPtr<double>(), other.getDataPtr<int32 >(), length); break;
        case DTYP::FLOAT  : _multiply_mat( getDataPtr<float >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::INT    : _multiply_mat( getDataPtr<int32 >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::UCHAR  : _multiply_mat( getDataPtr<uchar >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::CMPLX  : _multiply_mat( getDataPtr<cmplx >(), other.getDataPtr<int32 >(), length); break;
        }
    }else if( othrdatT == DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : _multiply_mat( getDataPtr<double>(), other.getDataPtr<uchar >(), length); break;
        case DTYP::FLOAT  : _multiply_mat( getDataPtr<float >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::INT    : _multiply_mat( getDataPtr<int32 >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::UCHAR  : _multiply_mat( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::CMPLX  : _multiply_mat( getDataPtr<cmplx >(), other.getDataPtr<uchar >(), length); break;
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : _multiply_mat( getDataPtr<double>(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::FLOAT  : _multiply_mat( getDataPtr<float >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::INT    : _multiply_mat( getDataPtr<int32 >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::UCHAR  : _multiply_mat( getDataPtr<uchar >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::CMPLX  : _multiply_mat( getDataPtr<cmplx >(), other.getDataPtr<cmplx >(), length); break;
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
        case DTYP::DOUBLE : _divide_mat( getDataPtr<double>(), other.getDataPtr<double>(), length); break;
        case DTYP::FLOAT  : _divide_mat( getDataPtr<float >(), other.getDataPtr<double>(), length); break;
        case DTYP::INT    : _divide_mat( getDataPtr<int32 >(), other.getDataPtr<double>(), length); break;
        case DTYP::UCHAR  : _divide_mat( getDataPtr<uchar >(), other.getDataPtr<double>(), length); break;
        case DTYP::CMPLX  : _divide_mat( getDataPtr<cmplx >(), other.getDataPtr<double>(), length); break;
        }
    }else if( othrdatT == DTYP::FLOAT){
        switch ( datT ){
        case DTYP::DOUBLE : _divide_mat( getDataPtr<double>(), other.getDataPtr<float >(), length); break;
        case DTYP::FLOAT  : _divide_mat( getDataPtr<float >(), other.getDataPtr<float >(), length); break;
        case DTYP::INT    : _divide_mat( getDataPtr<int32 >(), other.getDataPtr<float >(), length); break;
        case DTYP::UCHAR  : _divide_mat( getDataPtr<uchar >(), other.getDataPtr<float >(), length); break;
        case DTYP::CMPLX  : _divide_mat( getDataPtr<cmplx >(), other.getDataPtr<float >(), length); break;
        }
    }else if( othrdatT == DTYP::INT){
        switch ( datT ){
        case DTYP::DOUBLE : _divide_mat( getDataPtr<double>(), other.getDataPtr<int32 >(), length); break;
        case DTYP::FLOAT  : _divide_mat( getDataPtr<float >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::INT    : _divide_mat( getDataPtr<int32 >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::UCHAR  : _divide_mat( getDataPtr<uchar >(), other.getDataPtr<int32 >(), length); break;
        case DTYP::CMPLX  : _divide_mat( getDataPtr<cmplx >(), other.getDataPtr<int32 >(), length); break;
        }
    }else if( othrdatT == DTYP::UCHAR){
        switch ( datT ){
        case DTYP::DOUBLE : _divide_mat( getDataPtr<double>(), other.getDataPtr<uchar >(), length); break;
        case DTYP::FLOAT  : _divide_mat( getDataPtr<float >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::INT    : _divide_mat( getDataPtr<int32 >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::UCHAR  : _divide_mat( getDataPtr<uchar >(), other.getDataPtr<uchar >(), length); break;
        case DTYP::CMPLX  : _divide_mat( getDataPtr<cmplx >(), other.getDataPtr<uchar >(), length); break;
        }
    }else if( othrdatT == DTYP::CMPLX){
        switch ( datT ){
        case DTYP::DOUBLE : _divide_mat( getDataPtr<double>(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::FLOAT  : _divide_mat( getDataPtr<float >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::INT    : _divide_mat( getDataPtr<int32 >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::UCHAR  : _divide_mat( getDataPtr<uchar >(), other.getDataPtr<cmplx >(), length); break;
        case DTYP::CMPLX  : _divide_mat( getDataPtr<cmplx >(), other.getDataPtr<cmplx >(), length); break;
        }
    }
    return *this;
}

Mat& Mat::plusScalar(const double scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _plus_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _plus_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _plus_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _plus_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _plus_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::plusScalar(const float scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _plus_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _plus_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _plus_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _plus_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _plus_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::plusScalar(const int32 scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _plus_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _plus_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _plus_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _plus_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _plus_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}
Mat& Mat::plusScalar(const uchar scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _plus_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _plus_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _plus_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _plus_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _plus_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::minusScalar(const double scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _minus_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _minus_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _minus_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _minus_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _minus_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::minusScalar(const float scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _minus_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _minus_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _minus_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _minus_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _minus_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::minusScalar(const int32 scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _minus_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _minus_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _minus_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _minus_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _minus_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}
Mat& Mat::minusScalar(const uchar scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _minus_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _minus_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _minus_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _minus_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _minus_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::multiplyScalar(const double scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _multiply_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _multiply_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _multiply_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _multiply_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _multiply_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::multiplyScalar(const float scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _multiply_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _multiply_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _multiply_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _multiply_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _multiply_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::multiplyScalar(const int32 scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _multiply_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _multiply_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _multiply_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _multiply_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _multiply_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}
Mat& Mat::multiplyScalar(const uchar scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _multiply_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _multiply_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _multiply_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _multiply_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _multiply_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::dividedByScalar(const double scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _divided_by_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _divided_by_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _divided_by_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _divided_by_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _divided_by_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::dividedByScalar(const float scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _divided_by_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _divided_by_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _divided_by_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _divided_by_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _divided_by_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::dividedByScalar(const int32 scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _divided_by_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _divided_by_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _divided_by_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _divided_by_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _divided_by_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}
Mat& Mat::dividedByScalar(const uchar scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _divided_by_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _divided_by_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _divided_by_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _divided_by_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _divided_by_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::dividingScalar(const double scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _dividing_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _dividing_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _dividing_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _dividing_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _dividing_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::dividingScalar(const float scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _dividing_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _dividing_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _dividing_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _dividing_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _dividing_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat& Mat::dividingScalar(const int32 scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _dividing_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _dividing_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _dividing_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _dividing_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _dividing_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}
Mat& Mat::dividingScalar(const uchar scalar){
    if(isEmpty()) return *this;

    switch(datT){
    case DTYP::DOUBLE : _dividing_scalar(getDataPtr<double>(), scalar, length); break;
    case DTYP::FLOAT  : _dividing_scalar(getDataPtr<float >(), scalar, length); break;
    case DTYP::INT    : _dividing_scalar(getDataPtr<int32 >(), scalar, length); break;
    case DTYP::UCHAR  : _dividing_scalar(getDataPtr<uchar >(), scalar, length); break;
    case DTYP::CMPLX  : _dividing_scalar(getDataPtr<cmplx >(), scalar, length);
    }
    return *this;
}

Mat Mat::max(){
    if(isEmpty()) return Mat();
    switch(datT){
    case DTYP::DOUBLE : return _max<double>();
    case DTYP::FLOAT  : return _max<float >();
    case DTYP::INT    : return _max<int32 >();
    case DTYP::UCHAR  : return _max<uchar >();
    default           : return _max<cmplx >(); //case DTYP::CMPLX  : return _max<cmplx >();
    }
}

Mat Mat::min(){
    if(isEmpty()) return Mat();

    switch(datT){
    case DTYP::DOUBLE : return _min<double>();
    case DTYP::FLOAT  : return _min<float >();
    case DTYP::INT    : return _min<int32 >();
    case DTYP::UCHAR  : return _min<uchar >();
    default           : return _min<cmplx >(); //case DTYP::CMPLX  : return _min<cmplx >();
    }
}

Mat Mat::mean(){
    if(isEmpty()) return Mat();

    switch(datT){
    case DTYP::DOUBLE : return _mean<double>();
    case DTYP::FLOAT  : return _mean<float >();
    case DTYP::INT    : return _mean<int32 >();
    case DTYP::UCHAR  : return _mean<uchar >();
    default           : return _mean<cmplx >(); //case DTYP::CMPLX  : return _mean<cmplx >();
    }
}

Mat Mat::std(){
    if(isEmpty()) return Mat();

    switch(datT){
    case DTYP::DOUBLE : return _std<double>();
    case DTYP::FLOAT  : return _std<float >();
    case DTYP::INT    : return _std<int32 >();
    case DTYP::UCHAR  : return _std<uchar >();
    default           : return _std<cmplx >(); //case DTYP::CMPLX  : return _std<cmplx >();
    }
}

Mat Mat::sum(){
    if(isEmpty()) return Mat();

    switch(datT){
    case DTYP::DOUBLE : return _sum<double>();
    case DTYP::FLOAT  : return _sum<float >();
    case DTYP::INT    : return _sum<int32 >();
    case DTYP::UCHAR  : return _sum<uchar >();
    default           : return _sum<cmplx >(); //case DTYP::CMPLX  : return _sum<cmplx >();
    }
}

Mat Mat::repeat(const Mat& src, const uint32 rr, const uint32 rc, const uint32 rch){
    if(src.isEmpty()) return Mat();

    if(src.getChannel() > 1){
        fprintf(stderr,"src Mat argument of repeat method should have only 1 channel\n");
        return Mat();
    }

    switch(src.getDatType()){
    case DTYP::DOUBLE : return _repeat<double>(src, rr, rc, rch);
    case DTYP::FLOAT  : return _repeat<float >(src, rr, rc, rch);
    case DTYP::INT    : return _repeat<int32 >(src, rr, rc, rch);
    case DTYP::UCHAR  : return _repeat<uchar >(src, rr, rc, rch);
    default           : return _repeat<cmplx >(src, rr, rc, rch); // case DTYP::CMPLX
    }
}

} // end of jmat namespace
