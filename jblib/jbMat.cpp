#include <stdio.h>
#include "jbMat.h"
#include <float.h>
#include <iostream>

#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

int jbMat::instant_count = 0;

//-- shallow copy version & using shared_ptr
void jbMat::alloc(int len){

    if(len < 0){
        mA = nullptr;
        dat_ptr = nullptr;
    }else{
        mA = ptr_double (new double[len], std::default_delete<double[]>());
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
void jbMat::init(int r, int c, int ch){
    row = r;
    col = c;
    Nch = ch;
    lenRowCol = r*c;
    length = lenRowCol * ch;
    alloc(length);
}
void jbMat::initName(){
    obj_name =std::string( "jbMat_") + std::to_string (instant_count);
    instant_count++;
}

jbMat::jbMat():mA(nullptr){
    init(0,0,0);
    initName();

}
jbMat::jbMat(int rc ):mA(nullptr){
    init(rc, rc, 1);
    initName();
}
jbMat::jbMat(int r, int c, int ch):mA(nullptr){
    init(r, c, ch);
    initName();
}

jbMat::jbMat(int r, int c, int ch, std::string name):mA(nullptr){
    init(r, c, ch);
    obj_name = name;
}

jbMat::jbMat(const jbMat& mat){
    row = mat.getRow();
    col = mat.getCol();
    Nch = mat.getChannel();
    length    = mat.getLength();
    lenRowCol = row*col;
    mA      = mat.getMat();
    sync_data_ptr();
    initName();

    /* -- deep copying
    alloc(length);

    double *pt_matdat = mat.getMat().get();
    double *pt_thisdat = mA.get();
    //for(int i=0; i < length; i++)
    //    pt_thisdat[i] = pt_matdat[i];
    std::copy(pt_matdat,pt_matdat+mat.length,pt_thisdat);
    */
    fprintf(stdout,"copy constructor\n");
}

jbMat::jbMat( std::initializer_list<double> list ){
    //-- Making a vector by column vector type    
    init(list.size(),1,1);

    dat_ptr = mA.get();
    if(mA!=nullptr){
        std::vector<double> v;
        v.insert(v.end(),list.begin(),list.end());

        for(int i=0;i<length;i++){
            dat_ptr[i] = v.at(i);
        }
    }
}
jbMat::jbMat( std::initializer_list<int> list ){
    //-- Making a vector by column vector type
    init(list.size(),1,1);

    dat_ptr = mA.get();
    if(mA!=nullptr){
        std::vector<int> v;
        v.insert(v.end(),list.begin(),list.end());

        for(int i=0;i<length;i++){
            dat_ptr[i] = v.at(i);
        }
    }
}
jbMat::jbMat( std::initializer_list<float> list ){
    //-- Making a vector by column vector type
    init(list.size(),1,1);

    dat_ptr = mA.get();
    if(mA!=nullptr){
        std::vector<float> v;
        v.insert(v.end(),list.begin(),list.end());

        for(int i=0;i<length;i++){
            dat_ptr[i] = v.at(i);
        }
    }
}

jbMat::~jbMat(){}

void jbMat::setRowCol(int r, int c, int ch){
    int lenrc = r*c;
    int len   = lenrc*ch;

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
jbMat& jbMat::operator=( jbMat other){

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
jbMat& jbMat::operator=(const jbMat& other){

    fprintf(stdout,"Assign operator\n");
    row = other.getRow();
    col = other.getCol();
    Nch = other.getChannel();
    length = other.getLength();
    lenRowCol = row*col;
    mA  = other.getMat();
    sync_data_ptr();

    return *this;
}

jbMat jbMat::operator+(const jbMat& other){
    if(this->row != other.getRow() || this->col != other.getCol()|| this->Nch != other.getChannel() ){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return jbMat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return *this;
    }
    jbMat sum(this->row, this->col, this->Nch);
    for(int k=0; k < this->length ; k++)
        sum[k] = (*this)[k] + other[k];

    return sum;
}
jbMat jbMat::operator+(const double scalar){
    if(this->isEmpty()){
        jbMat tmp(1,1,1);
        tmp[0] = scalar;
        return tmp;
    }else{
        jbMat sum(this->row, this->col, this->Nch);
        for(int k=0; k < this->length ; k++)
            sum[k] = (*this)[k] + scalar;

        return sum;
    }
}
jbMat& jbMat::operator+=(const jbMat& other){

    if(this->row != other.getRow() || this->col != other.getCol() || this->Nch != other.getChannel()){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return *this;
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return *this;
    }

    for(int k=0; k < this->length ; k++)
        (*this)[k] += other[k];

    return *this;
}
jbMat& jbMat::operator+=(const double scalar){
    if(this->isEmpty()){
        this->setRowCol(1,1,1);
        (*this)[0] = scalar;
        return *this;
    }else{
        for(int k=0; k < this->length ; k++)
            (*this)[k] += scalar;

        return *this;
    }
}
jbMat jbMat::operator-(const jbMat& other){
    if(this->row != other.getRow() || this->col != other.getCol() || this->Nch != other.getChannel() ){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return jbMat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return *this;
    }
    jbMat sum(this->row, this->col, this->Nch);
    for(int k=0; k < this->length ; k++)
        sum[k] = (*this)[k] - other[k];

    return sum;
}
jbMat jbMat::operator-(const double scalar){
    if(this->isEmpty()){
        jbMat tmp(1,1,1);
        tmp[0] = -scalar;
        return tmp;
    }else{
        jbMat sum(this->row, this->col, this->Nch);
        for(int k=0; k < this->length ; k++)
            sum[k] = (*this)[k] - scalar;

        return sum;
    }
}
jbMat& jbMat::operator-=(const jbMat& other){

    if(this->row != other.getRow() || this->col != other.getCol() || this->Nch != other.getChannel()){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return *this;
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return *this;
    }

    for(int k=0; k < this->length ; k++)
        (*this)[k] -= other[k];

    return *this;
}
jbMat& jbMat::operator-=(const double scalar){
    if(this->isEmpty()){
        this->setRowCol(1,1,1);
        (*this)[0] = scalar;
        return *this;
    }else{
        for(int k=0; k < this->length ; k++)
            (*this)[k] -= scalar;

        return *this;
    }
}
jbMat jbMat::operator*(const jbMat& other){
    if(this->row != other.getRow() || this->col != other.getCol() || this->Nch != other.getChannel() ){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return jbMat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return *this;
    }
    jbMat sum(this->row, this->col, this->Nch);
    for(int k=0; k < this->length ; k++)
        sum[k] = (*this)[k] * other[k];

    return sum;
}
jbMat jbMat::operator*(const double scalar){
    if(this->isEmpty()){
        jbMat tmp(1,1,1);
        tmp[0] = -scalar;
        return tmp;
    }else{
        jbMat sum(this->row, this->col, this->Nch);
        for(int k=0; k < this->length ; k++)
            sum[k] = (*this)[k] * scalar;

        return sum;
    }
}
jbMat& jbMat::operator*=(const jbMat& other){

    if(this->row != other.getRow() || this->col != other.getCol() || this->Nch != other.getChannel() ){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return *this;
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return *this;
    }

    for(int k=0; k < this->length ; k++)
        (*this)[k] *= other[k];

    return *this;
}
jbMat& jbMat::operator*=(const double scalar){
    if(this->isEmpty()){
        this->setRowCol(1,1,1);
        (*this)[0] = scalar;
        return *this;
    }else{
        for(int k=0; k < this->length ; k++)
            (*this)[k] *= scalar;

        return *this;
    }
}
jbMat jbMat::operator/(const jbMat& other){
    if(this->row != other.getRow() || this->col != other.getCol() || this->Nch != other.getChannel() ){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return jbMat();
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return *this;
    }
    jbMat sum(this->row, this->col, this->Nch);
    for(int k=0; k < this->length ; k++)
        sum[k] = (*this)[k] / other[k];

    return sum;
}
jbMat jbMat::operator/(const double scalar){
    if(this->isEmpty()){
        jbMat tmp(1,1,1);
        tmp[0] = -scalar;
        return tmp;
    }else{
        jbMat sum(this->row, this->col, this->Nch);
        for(int k=0; k < this->length ; k++)
            sum[k] = (*this)[k] / scalar;

        return sum;
    }
}
jbMat& jbMat::operator/=(const jbMat& other){

    if(this->row != other.getRow() || this->col != other.getCol() || this->Nch != other.getChannel()){
        fprintf(stdout,"Stop adding two operands because the both operands are not same size\n");
        return *this;
    }else if(this->isEmpty()){
        fprintf(stdout,"Stop adding two operands because the both operands are empty\n");
        return *this;
    }

    for(int k=0; k < this->length ; k++)
        (*this)[k] /= other[k];

    return *this;
}
jbMat& jbMat::operator/=(const double scalar){
    if(this->isEmpty()){
        this->setRowCol(1,1,1);
        (*this)[0] = scalar;
        return *this;
    }else{
        for(int k=0; k < this->length ; k++)
            (*this)[k] /= scalar;

        return *this;
    }
}
double& jbMat::operator[] (int i) const{
    if(isEmpty()) return *mA.get(); //*mA;

    if(i >= length){
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        i = length-1;
    }

    return (mA.get())[i];
}
double& jbMat::operator() (int i) const{
    if(isEmpty()) return *mA.get(); //*mA;
    if(i >= length){
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        i = length-1;
    }
    return (mA.get())[i];
}
double& jbMat::operator() (int r, int c) const{
    if(isEmpty()) return *mA.get(); //*mA;
    int idx = r*col + c;
    if(idx >= length) {
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        idx = length-1;
    }

    return (mA.get())[idx];
}

double& jbMat::operator() (int r, int c, int ch) const{
    if(isEmpty()) return *mA.get(); //*mA;
    int idx = ch*lenRowCol + (r*col + c);
    if(idx >= length) {
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        idx = length-1;
    }

    return (mA.get())[idx];
}

int jbMat::reshape(int r, int c, int ch){
    int rc   = r*c;
    int tlen = rc*ch;
    if( tlen != length){
        fprintf(stderr," reshape argument is not correct!\n");
        return -1;
    }
    /*
    int itch, itrc, k;
    if( ch != Nch ){
        double *tmA;
        double *mdat=mA.get();
        try{
            tmA = new double[static_cast<unsigned long>(length)];
        }catch(std::bad_alloc& ex){
            fprintf(stderr,"memory allocation error in jbMat::reshape() : %s\n",ex.what());
            tmA = nullptr;
            return -1;
        }
        k=0;
        for (int ich=0; ich < Nch; ich++){
            for(int irc=0; irc < row*col ; irc++){
                itch = k / rc;
                itrc = (k - itch*rc)*ch;
                tmA[itch + itrc] = mdat[ich+irc*Nch];
                k++;
            }
        }
        mA.reset(tmA,std::default_delete<double[]>());

    }
    */
    row = r;
    col = c;
    Nch = ch;
    lenRowCol = rc;

    return 0;
}

void jbMat::transpose(){
    if(isEmpty()){
        fprintf(stderr," Transpose: This jbMat is empty\n");
        return ;
    }

    double *tmA;
    try{
        tmA= new double[static_cast<unsigned long>(length)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"Transpose Error: %s\n",ex.what());
        return;
    }
    double *mdat = mA.get();
    int i, j, k, ch_offset;
    for(k=0; k < Nch; k++){
        ch_offset = k* lenRowCol;
        for(i=0; i < row; i++){
            for(j=0; j<col; j++)
                tmA[ch_offset + j*row+i] = mdat[ch_offset + i*col+j];
        }
    }
    mA.reset(tmA,std::default_delete<double[]>());
    sync_data_ptr();

    int row_tr = col;
    col = row;
    row = row_tr;
}

void jbMat::printMat() const {
    printMat(obj_name);
}
void jbMat::printMat(const std::string objname) const
{
    const int bufsz = 2049;
    char buf[bufsz]="\0";
    char tmp[bufsz];
    int i,j;

    const double neg_max_double = -DBL_EPSILON ;
    const double pos_min_double = DBL_EPSILON ;

    double val;
    int k;
    double* mdat = mA.get();
    /*
    for( k=0; k < Nch; k++){      
        fprintf(stdout,"channel: %d \n",k);

        for( i=0; i< row*col*Nch; i+= col*Nch){
            snprintf(buf,bufsz,"[");
            for( j=0; j< col*Nch; j+= Nch){
                val = mdat[i+j+k];
                if( val >= neg_max_double && val <= pos_min_double)
                    val = 0.0;
                snprintf(tmp,bufsz," %.4f ",val);
                strncat(buf,tmp,bufsz);
            }
            strncat(buf,"]",1);
            fprintf(stdout,"%s\n",buf);
        }
    } */
    if(!objname.empty())
        fprintf(stdout,"object : %s \n", objname.c_str());

    int ch_offset;
    for( k=0; k < Nch; k++){
        fprintf(stdout,"channel: %d \n",k);
        ch_offset = k*lenRowCol;
        for( i=0; i < lenRowCol; i+=col){ // rows
            snprintf(buf,bufsz,"[");
            for( j=0; j < col; j++){ // columns
                val = mdat[i+j+ch_offset];
                if( val >= neg_max_double && val <= pos_min_double)
                    val = 0.0;
                snprintf(tmp,bufsz," %.4f ",val);
                strncat(buf,tmp,bufsz);
            }
            strncat(buf,"]",1);
            fprintf(stdout,"%s\n",buf);
        }
    }
}

jbMat jbMat::copy() const{
    jbMat A(row, col, Nch);

    double *pt_matdat  = A.getMat().get();
    double *pt_thisdat = this->mA.get();

    std::copy(pt_thisdat,pt_thisdat+length,pt_matdat);

    return A;
}

jbMat jbMat::ones(int r, int c, int ch){
    if( r <= 0 || c <= 0 || ch <= 0){
        fprintf(stdout,"In ones method: arguments r , c and ch are to be larger than 0 ");
        return jbMat();
    }

    jbMat A(r, c, ch);

    double* pt_dat = A.getMat().get();
    for(int i=0; i < r*c*ch; i++){
        pt_dat[i] = 1.0;
    }

    return A;
}


jbMat jbMat::zeros(int r, int c, int ch){
    if( r <= 0 || c <= 0 || ch <= 0){
        fprintf(stdout,"In zeros method: arguments r , c and ch are to be larger than 0 ");
        return jbMat();
    }

    jbMat A(r, c, ch);

    double* pt_dat = A.getMat().get();
    for(int i=0; i < r*c*ch; i ++){
        pt_dat[i] = 0.0;
    }

    return A;
}

/*
jbMat jbMat::getChannelN(const unsigned int NoCh){
    jbMat A(*this);
    int numCh = NoCh;
    if(numCh >= A.getChannel()){
        fprintf(stdout,"channel number of getNchannel() is out of bound\n The last channel is selected\n");
        numCh = A.getChannel()-1;
    }

    A.mA = A.mA+A.lenRowCol*numCh;
    return A;
}*/

jbMat jbMat::copyChannelN(const unsigned int NoCh) const{
    jbMat A(row,col,1);

    int numCh = NoCh;
    if(numCh >= A.getChannel()){
        fprintf(stdout,"channel number of getNchannel() is out of bound\n The last channel is selected\n");
        numCh = A.getChannel()-1;
    }

    double *srcDat_pt = mA.get();
    double *tarDat_pt = A.mA.get();
    int offset = numCh*lenRowCol;
    for(int i=0; i < lenRowCol; i++ )
        tarDat_pt[i] = srcDat_pt[offset+i];

    return A;
}

void jbMat::setChannelN(const jbMat& src, const unsigned int srcCh, const unsigned int tarCh){
    if(src.getChannel()-1 < srcCh){
        fprintf(stdout,"setChannelN(): src argument has less channel than argument srcCh\n");
        return ;
    }

    int tar_ch;
    if(isEmpty()){
        init(src.getRow(),src.getCol(), 1);
        tar_ch = 0;
    }else if( src.getRow() != row && src.getCol() != col && src.getChannel() != Nch){
        fprintf(stdout,"*this and src are not equal size. *this:(%d, %d, %d) != src(%d, %d, %d)\n",row,col,Nch,src.getRow(),src.getCol(),src.getChannel());
        return ;
    }else
        tar_ch = tarCh;

    sync_data_ptr();
    double *srcdat_ptr = src.getMat().get();
    int src_offset = srcCh*lenRowCol;
    int tar_offset = tar_ch*lenRowCol;

    for(int i=0; i < lenRowCol; i++ )
        dat_ptr[tar_offset+i] = srcdat_ptr[src_offset+i];

}

void jbMat::setChannelN(const jbMat& src, const unsigned int srcFromCh,const unsigned int Channels, const unsigned int tarToCh){
    if(src.getChannel()-1 < srcFromCh+Channels ){
        fprintf(stdout,"setChannelN(): srcFromCh and Channels are not correct! \n");
        return ;
    }

    int tar_ch;
    if(isEmpty()){
        init(src.getRow(),src.getCol(), Channels );
        tar_ch = 0;
    }else if( src.getRow() != row && src.getCol() != col && src.getChannel() != Nch){
        fprintf(stdout,"*this and src are not equal size. *this:(%d, %d, %d) != src(%d, %d, %d)\n",row,col,Nch,src.getRow(),src.getCol(),src.getChannel());
        return ;
    }else
        tar_ch = tarToCh;

    sync_data_ptr();
    double *srcdat_ptr = src.getMat().get();
    int src_offset = srcFromCh*lenRowCol;
    int tar_offset = tar_ch*lenRowCol;

    for(int j=0; j < Channels ; j++){
        for(int i=0; i < lenRowCol; i++ )
            dat_ptr[tar_offset+i] = srcdat_ptr[src_offset+i];
        src_offset += lenRowCol;
        tar_offset += lenRowCol;
    }

}

void jbMat::setName(std::string name){
    obj_name = name;
}
