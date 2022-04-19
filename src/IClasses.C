
#include "IClasses.h"

using namespace std;

const int OFFSET = 1;
  
  /*  the indexing and storage allocation for vectors and matrices follow the 
        prescriptions of Numerical Recipes in C, by 
            Press, Teukolsky, Vetterling, and Flannery  */

IVector::IVector() { 
    therror(" to construct an IVector, use the constructor with parameters");
}

IVector::IVector(int xnl, int xnh){
       init(xnl,xnh);
       zero();
}

IVector::IVector(const IVector & A){
       init(A.nl, A.nh);
       load(A);
}

IVector::~IVector(){ 
	destroy();
}

IVector & IVector::operator = (const IVector & A) {
	if (A.nl != nl || A.nh != nh){
	  therror(" attempt to assign from an incompatible IVector");
        }
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = A.ptr[i];
	return * this;
}


int IVector::l() const {return nl;}
int IVector::h() const {return nh;} 

void IVector::fill(int n, ...) {
	if (n != (1+ nh-nl)) therror(" wrong size argument for fillVector ");
        va_list a;
        va_start(a,n);
	int i;
        for (i = nl; i<= nh; i++) ptr[i] = va_arg(a,int);
        va_end(a);
}

void IVector::zero() { 
	int i;
	for (i=nl; i<= nh; i++)  ptr[i] = 0;
}


IVector IVector::operator - () const {
        IVector A(nl,nh);
        int i;
        for (i = nl; i<= nh; i++) A.ptr[i] = -ptr[i]; 
        return A;
}

IVector & IVector::operator *= (int  A) {
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] * A;
	return * this;
}

IVector & IVector::operator += (const IVector & A) {
        if (A.nl != nl || A.nh != nh){
	  therror(" attempt to add incompatible IVectors ");}
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] +  A.ptr[i];
	return * this;
}
 
IVector & IVector::operator -= (const IVector & A) {
	if (A.nl != nl || A.nh != nh){
	  therror(" attempt to subtract incompatible IVectors ");}
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] - A.ptr[i];
	return * this;
}

 
void IVector::init(int xnl,int xnh){
      nl = xnl; nh = xnh;
      if ((nh-nl)<0) 
         therror(" the number of components of a IVector should be positive ");
      ptr = new int[nh-nl+1 + OFFSET];
      if (!ptr) therror("allocation failure in IVector()");
      ptr -= (nl - OFFSET); 
}

void IVector::load(const IVector & A) {
	int i;
	for (i=nl; i<= nh; i++)  ptr[i] = A.ptr[i];
}

void IVector::destroy() {
	ptr += (nl-OFFSET);
	delete [] ptr;
}	



std::ostream & operator << (std::ostream & os, const IVector & A){
    int i;
    os << "<";
    for (i=A.nl;i<= A.nh;i++) os << "\t" << A.ptr[i];
    os << ">" << endl;
    return os;
}


int operator ==  (const IVector &A,const IVector &B) {
  if ( A.nl != B.nl  || A.nh != B.nh) return 0;
  if ( A.ptr == B.ptr) return 1;
  int i;
  for (i=A.nl; i<= A.nh; i++) {  if (A.ptr[i] != B.ptr[i]) return 0;};
  return 1;
}


int operator !=  (const IVector &A,const IVector &B) {
  return (A == B ? 0 : 1);
}

IVector operator + (const IVector & A, const IVector & B){
  if (A.nl != B.nl || A.nh != B.nh){
    therror(" attempt to add incompatible IVectors ");}
  IVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] + B.ptr[i];
  return C; 
}
    
IVector operator - (const IVector & A, const IVector & B){
  if (A.nl != B.nl || A.nh != B.nh){
    therror(" attempt to subtract incompatible IVectors ");}
  IVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] - B.ptr[i];
  return C; 
}

IVector operator * (int A, const IVector & B){
  IVector C(B.nl,B.nh);
  int i;
  for (i=B.nl; i<= B.nh; i++) C.ptr[i] = A * B.ptr[i];
  return C; 
}

IVector operator * (const IVector & A, int  B){
  IVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] * B;
  return C; 
}

int  operator * (const IVector & A, const IVector & B){
  if (A.nl != B.nl || A.nh != B.nh){
    therror(" attempt to dot incompatible IVectors ");}
  int s = 0;
  int i;
  for (i=A.nl; i<= A.nh;i++) s += A.ptr[i] * B.ptr[i];
  return s;
}



IMatrix::IMatrix() { 
      therror("to construct an IMatrix, use the constructor with parameters");
}

IMatrix::IMatrix(int xml, int xmh, int xnl, int xnh){
       init(xml,xmh,xnl,xnh);
       zero();
}

IMatrix::IMatrix(const IMatrix & A){
       init(A.ml,A.mh,A.nl, A.nh);
       load(A);
}

IMatrix::IMatrix(int xnl, int xnh, int a) {
       init(xnl,xnh,xnl,xnh);
       int i,j;
       for (i=ml; i<= mh; i++) {
	 for (j=nl;j<=nh; j++) ptr[i][j] = 0;
	 ptr[i][i] =  a;
       }
}

IMatrix::~IMatrix(){ 
	destroy();
}

IMatrix &  IMatrix::operator = (const IMatrix & A) {
	if (A.ml != ml || A.mh != mh ||A.nl != nl || A.nh != nh){
	  therror(" attempt to assign from an incompatible IMatrix");
        }
	int i,j;
	for (i = ml; i<=mh; i++) 
           for (j = nl; j<=nh; j++)  ptr[i][j] = A.ptr[i][j];
	return * this;
}

int IMatrix::rl() const {return ml;}
int IMatrix::rh() const {return mh;}
int IMatrix::cl() const {return nl;}
int IMatrix::ch() const {return nh;} 

void IMatrix::fillRow(int m, int n, ...) {
	if (n != (1+ nh-nl)) therror(" wrong size argument for fillVector ");
        va_list a;
        va_start(a,n);
	int i;
        for (i = nl; i<= nh; i++) ptr[m][i] = va_arg(a,int);
        va_end(a);
}

void IMatrix::readRow(int n,const IVector & A){
        if (A.nl != nl || A.nh != nh)
         therror(" attempt to fill a IMatrix with an incompatible IVector");
        int i;
        for (i= nl; i<= nh; i++)  ptr[n][i] = A.ptr[i];
}

void IMatrix::readColumn(int n, const IVector & A){
        if (A.nl != ml || A.nh != mh)
         therror(" attempt to fill a IMatrix with an incompatible IVector");
        int i;
        for (i= ml; i<= mh; i++)  ptr[i][n] = A.ptr[i];
}

void IMatrix::zero() {
	int i,j;
	for (i=ml; i<= mh; i++)
	  for (j=nl; j<= nh; j++)  ptr[i][j] = 0;
}

IVector IMatrix::row(int i) const{
	IVector A(nl,nh);
	int j;
        for (j=nl;j<=nh;j++) A[j] = ptr[i][j];
	return A;
}

IVector IMatrix::column(int j) const{
	IVector A(ml,mh);
	int i;
        for (i=ml;i<=mh;i++) A[i] = ptr[i][j];
	return A;
}

IMatrix IMatrix::operator - () const {
        IMatrix A(ml,mh,nl,nh);
	int i,j;
	for (i = ml; i<=mh; i++)
	  for (j=nl;j<=nh;j++) A.ptr[i][j] = - ptr[i][j];
	return A;
}

IMatrix & IMatrix::operator *= (int A) {
	int i,j;
	for (i = ml; i<=mh; i++)
	  for (j=nl;j<=nh;j++) ptr[i][j] = ptr[i][j] * A;
	return * this; 
}

IMatrix & IMatrix::operator += (const IMatrix & A) {
        if (A.ml != ml || A.mh != mh ||  A.nl != nl || A.nh != nh){
	  therror(" attempt to add incompatible IMatrices ");}
	int i,j;
	for (i = ml; i<=mh; i++) 
	  for (j= nl;j<= nh; j++) ptr[i][j] = ptr[i][j] +  A.ptr[i][j];
	return * this;
}

IMatrix & IMatrix::operator -= (const IMatrix & A) {
        if (A.ml != ml || A.mh != mh ||  A.nl != nl || A.nh != nh){
	  therror(" attempt to subtract incompatible IMatrices ");}
	int i,j;
	for (i = ml; i<=mh; i++) 
	  for (j= nl;j<= nh; j++) ptr[i][j] = ptr[i][j] -  A.ptr[i][j];
	return * this;
}

void IMatrix::init(int xml, int xmh, int xnl,int xnh){
	ml = xml; mh = xmh; nl = xnl; nh = xnh;
        if ((mh-ml)<0) 
           therror(" the number of rows of a IMatrix should be positive ");
        if ((nh-nl)<0) 
           therror(" the number of columns of a IMatrix should be positive ");
	int nrow = mh-ml+1;
        int ncol = nh-nl+1;
	ptr = new int * [nrow + OFFSET];
	if (!ptr)  therror ("allocation failure in IMatrix()");
	ptr += OFFSET;
        ptr -= ml;    /*  now ptr has the value s.t. ptr + ml is the correct
                              starting address   */
	ptr[ml] = new int[nrow*ncol + OFFSET];
	if (!ptr[ml]) therror("allocation failure, stage 2, in IMatrix()");
	ptr[ml] += OFFSET;
	ptr[ml] -= nl;
             /* now evaluate the pointers to rows  */
	for (int i = ml+1; i<= mh; i++) ptr[i] = ptr[i-1] + ncol;
} 

void IMatrix::load(const IMatrix & A) {
	int i,j;
	for (i=ml; i<= mh; i++)
	  for (j=nl; j<= nh; j++)  ptr[i][j] = A.ptr[i][j];
}

void IMatrix::destroy() {
	int * w = ptr[ml] + nl - OFFSET;
	delete [] w;
	ptr += (ml - OFFSET);
	delete [] ptr;
}

ostream & operator << (ostream & os, const IMatrix & A){
    int i,j;
    for (i=A.ml;i<=A.mh;i++) {
       os << endl <<  "<";
       for (j=A.nl;j<= A.nh;j++) os << "\t" << A.ptr[i][j];
    os << "\t>" << endl;
    }
    os << endl;
    return os;
}

int operator ==  (const IMatrix &A,const IMatrix &B) {
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) return 0;
  if ( A.ptr == B.ptr) return 1;
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) if (A.ptr[i][j] != B.ptr[i][j]) return 0;
  return 1;
}

int operator !=  (const IMatrix &A,const IMatrix &B) {
  return (A == B ? 0 : 1);
}


IMatrix operator + (const IMatrix & A, const IMatrix & B){
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) 
    therror(" attempt to add incompatible IMatrices ");
  IMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j]+ B.ptr[i][j];
  return C; 
}
    
IMatrix operator - (const IMatrix & A, const IMatrix & B){
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) 
    therror(" attempt to subtract incompatible IMatrices ");
  IMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j] - B.ptr[i][j];
  return C;
} 

IMatrix operator * (int A, const IMatrix & B){
  IMatrix C(B.ml,B.mh,B.nl,B.nh);
  int i,j;
  for (i=B.ml; i<= B.mh; i++)  
    for (j=B.nl;j<=B.nh; j++) C.ptr[i][j] = A * B.ptr[i][j];
  return C; 
}

IMatrix operator * (const IMatrix & A, int B){
  IMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j] * B;
  return C; 
}

IVector operator * (const IVector & A, const IMatrix & B){
  if (A.nl != B.ml || A.nh != B.mh)
    therror(" attempt to multiply incompatible IVector and IMatrix ");
  IVector C(B.nl,B.nh);
  int i,j;
  for (j=B.nl; j<= B.nh; j++) {
    C.ptr[j] = 0;
    for (i=B.ml;i<=B.mh; i++) C.ptr[j] += A.ptr[i] * B.ptr[i][j];
  }
  return C; 
}

IVector operator * (const IMatrix & A, const IVector & B){
  if (A.nl != B.nl || A.nh != B.nh)
    therror(" attempt to multiply incompatible IMatrix and IVector ");
  IVector C(A.ml,A.mh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++) {
    C.ptr[i] = 0;
    for (j=B.nl;j<=B.nh; j++) C.ptr[i] += A.ptr[i][j] * B.ptr[j];
  }
  return C; 
}

IMatrix operator * (const IMatrix & A, const IMatrix & B){
  if (A.nl != B.ml || A.nh != B.mh)
    therror(" attempt to multiply incompatible DMatrices ");
  IMatrix C(A.ml,A.mh,B.nl,B.nh);
  int i,j,k;
  for (i=A.ml; i<= A.mh; i++) {
    for (j=B.nl; j<= B.nh; j++) {
      C.ptr[i][j] = 0;
      for (k=A.nl;k<=A.nh;k++) C.ptr[i][j] += A.ptr[i][k] * B.ptr[k][j];
    }
  }
  return C; 
}

IMatrix dyad(const IVector & A, const IVector & B){
  IMatrix  C(A.nl,A.nh,B.nl,B.nh);
  int i,j;
  for (i=A.nl; i<= A.nh; i++)
    for (j=B.nl;j<=B.nh; j++) C.ptr[i][j] = A.ptr[i] * B.ptr[j];
  return C; 
}


/*  methods for inverse,  determinant, etc.     */

IMatrix IMatrix::transpose() const {
  IMatrix A(nl,nh,ml,mh);
  int i,j;
  for (i=ml; i<= mh; i++)
    for (j=nl; j<= nh; j++) A.ptr[j][i] = ptr[i][j];
  return A;
}

int IMatrix::trace() const {
  if (mh != nh || ml != nl)
         therror(" attempt to take the trace of an asymmetric IMatrix ");
  int i;   int tr = 0;
  for (i = nl; i <= nh; i++)  tr += ptr[i][i];
  return tr;
}

void IMatrix::rowint(int i, int j) {
  int k;
  int el;
  for (k= nl; k<= nh; k++) {
    el = ptr[i][k];
    ptr[i][k] = ptr[j][k];
    ptr[j][k] = el;
  }
}

void IMatrix::colint(int i, int j) {
  int k;
  int el;
  for (k= ml; k<= mh; k++) {
    el = ptr[k][i];
    ptr[k][i] = ptr[k][j];
    ptr[k][j] = el;
  }
}


