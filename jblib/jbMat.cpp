#include "jbMat.h"
#include <stdio.h>
#include <float.h>
#include <iostream>
#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

//-- shallow copy version & using shared_ptr

void jbMat::alloc(int len){

    mA = ptr_double (new double[len], std::default_delete<double[]>());
    /*
    try{
        mA = (len==0)? nullptr : new double[len];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"memory allocation error: %s\n",ex.what());
        mA = nullptr;
    }
    */
}

jbMat::jbMat():mA(nullptr),row(0),col(0),Nch(0){length =0; lenRowCol=0; }
jbMat::jbMat(int rc ):mA(nullptr){
    jbMat(rc, rc, 1);
}
jbMat::jbMat(int r, int c, int ch):mA(nullptr),row(r),col(c),Nch(ch){
    lenRowCol = r*c;
    length    = lenRowCol*ch;

    alloc(length);
}

jbMat::jbMat(const jbMat& mat){
    row = mat.row;
    col = mat.col;
    Nch = mat.Nch;
    length    = mat.length;
    lenRowCol = mat.lenRowCol;
    mA = mat.getMat();

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
    length = list.size();
    row    = length;
    col    = 1;
    Nch    = 1;
    lenRowCol = length;

    alloc(length);
    double* pt_mdat = mA.get();
    if(mA!=nullptr){
        std::vector<double> v;
        v.insert(v.end(),list.begin(),list.end());

        for(int i=0;i<length;i++){
            pt_mdat[i] = v.at(i);
        }
    }
}
jbMat::jbMat( std::initializer_list<int> list ){
    //-- Making a vector by column vector type
    length = list.size();
    row    = length;
    col    = 1;
    Nch    = 1;
    lenRowCol= length;

    alloc(length);
    double* pt_mdat = mA.get();
    if(mA!=nullptr){
        std::vector<int> v;
        v.insert(v.end(),list.begin(),list.end());

        for(int i=0;i<length;i++){
            pt_mdat[i] = v.at(i);
        }
    }
}
jbMat::jbMat( std::initializer_list<float> list ){
    //-- Making a vector by column vector type
    length = list.size();
    row    = length;
    col    = 1;
    Nch    = 1;
    lenRowCol= length;

    alloc(length);
    double* pt_mdat = mA.get();
    if(mA!=nullptr){
        std::vector<float> v;
        v.insert(v.end(),list.begin(),list.end());

        for(int i=0;i<length;i++){
            pt_mdat[i] = v.at(i);
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

}


//-- overloading operators : it calls copy constructor
jbMat& jbMat::operator=(jbMat other){
/*
The parameter to the ‘operator=()’ is passed by value which calls copy constructor
to create an object local to the ‘operator=()’.
Than the value of the temp object is swapped with ‘*this’ object
*/
    fprintf(stdout,"Assign operator\n");
    row = other.row;
    col = other.col;
    Nch = other.Nch;
    length = other.getLength();
    lenRowCol = row*col;
    mA  = other.getMat();
/*
    std::swap(row,other.row);
    std::swap(col,other.col);
    std::swap(Nch,other.Nch);

    std::swap(length,other.length);
    std::swap(lenRowCol,other.lenRowCol);
    std::swap(mA,other.mA); // swapping mA pointer
*/

    /*
    if(this != &other){
        length = row*col;
        double *tmA = new double[length];
        std::copy(other_ma,other_ma+length, tmA);
        mA = tmA;
    } */

    return *this;
}

jbMat jbMat::operator+(const jbMat& other){
    if(this->row != other.row || this->col != other.col || this->Nch != other.Nch){
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

    if(this->row != other.row || this->col != other.col || this->Nch != other.Nch){
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
    if(this->row != other.row || this->col != other.col || this->Nch != other.Nch){
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

    if(this->row != other.row || this->col != other.col || this->Nch != other.Nch){
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
    if(this->row != other.row || this->col != other.col || this->Nch != other.Nch){
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

    if(this->row != other.row || this->col != other.col || this->Nch != other.Nch){
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
    if(this->row != other.row || this->col != other.col || this->Nch != other.Nch){
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

    if(this->row != other.row || this->col != other.col || this->Nch != other.Nch){
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
    int rc    = row * col;
    int chidx = i / rc; // chidx is used to start offset of data buffer
    int rcidx = chidx + (i - rc*chidx)*Nch;
    if(i >= length){
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        i = length-1;
    }

    //return *(mA+rcidx);
    return (mA.get())[rcidx];
}
double& jbMat::operator() (int i) const{
    if(isEmpty()) return *mA.get(); //*mA;
    int rc    = row * col;
    int chidx = i / rc; // chidx is used to start offset of data buffer
    int rcidx = chidx + (i - rc*chidx)*Nch;
    if(i >= length){
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        i = length-1;
    }
    //return *(mA+rcidx);
    return (mA.get())[rcidx];
}
double& jbMat::operator() (int r, int c) const{
    if(isEmpty()) return *mA.get(); //*mA;
    int idx = r*col + c;
    if(idx >= length) {
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        idx = length-1;
    }
    //return *(mA+idx);
    return (mA.get())[idx];
}

double& jbMat::operator() (int r, int c, int ch) const{
    if(isEmpty()) return *mA.get(); //*mA;
    int idx = ch + (r*col + c)*Nch;
    if(idx >= length) {
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        idx = length-1;
    }
    //return *(mA+idx);
    return (mA.get())[idx];
}

int jbMat::reshape(int r, int c, int ch){
    int rc   = r*c;
    int tlen = rc*ch;
    if( tlen != length){
        fprintf(stderr," reshape argument is not correct!\n");
        return -1;
    }
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
    int i, j, k;
    for(k=0; k < Nch; k++){
        for(i=0; i < row; i++){
            for(j=0; j<col; j++)
                tmA[k+(j*row+i)*Nch] = mdat[k+(i*col+j)*Nch];
        }
    }
    mA.reset(tmA,std::default_delete<double[]>());

    int rows = col;
    col = row;
    row = rows;
}

void jbMat::printMat() const
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

    }
}

jbMat jbMat::copy() {
    jbMat A(row, col, Nch);

    double *pt_matdat  = A.getMat().get();
    double *pt_thisdat = this->mA.get();

    std::copy(pt_thisdat,pt_thisdat+length,pt_matdat);

    return A;
}

jbMat jbMat::ones(int r, int c, int ch){
    if( r < 0 || c < 0 || ch < 0){
        fprintf(stdout,"In ones method: r , c and ch are to be larger than 0 ");
        return jbMat();
    }

    jbMat A(r, c, ch);

    double* pt_dat = A.getMat().get();
    for(int i=0; i < r*c*ch; i ++){
        pt_dat[i] = 1.0;
    }

    return A;
}


jbMat jbMat::zeros(int r, int c, int ch){
    if( r < 0 || c < 0 || ch < 0){
        fprintf(stdout,"In ones method: r , c and ch are to be larger than 0 ");
        return jbMat();
    }

    jbMat A(r, c, ch);

    double* pt_dat = A.getMat().get();
    for(int i=0; i < r*c*ch; i ++){
        pt_dat[i] = 0.0;
    }

    return A;
}
