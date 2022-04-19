
#ifndef CCLASSES_H
#define CCLASSES_H

/*   classes for complex vectors and matrices   */

#include "thbasics.h"
#include "complex.h"
#include "IClasses.h"
#include "DClasses.h"

typedef  complex*   cprt;  

class CMatrix;

class CVector{
 
 friend class CMatrix;

 protected:

     int nl, nh;
     complex * ptr;

 public:
    
    CVector();
    CVector(int nl, int nh);
    CVector(const CVector & A);
    CVector(const DVector & A);
    CVector(const IVector & A);

    ~CVector();
 
    CVector & operator = (const CVector & A);

    inline complex &  operator [](int i){  return ptr[i]; }
    int l() const;
    int h() const;

    void fill(int n, ...);
    void zero();
  
    friend std::ostream & operator << (std::ostream & os, const CVector & A);

    friend int operator == (const CVector & A, const CVector & B);
    friend int operator != (const CVector & A, const CVector & B);
    friend int eqacc (const CVector & A, const CVector & B, double acc);
        /*  tests equality  to accuracy acc    */
    
    friend CVector operator + (const CVector & A, const CVector & B);
    friend CVector operator - (const CVector & A, const CVector & B);
    friend CVector operator * (complex a, const CVector & B);
    friend CVector operator * (const CVector & A, complex b);
    friend CVector operator / (const CVector & A, complex b);
    friend complex  operator * (const CVector & A, const CVector & B);
    friend complex  norm(const CVector & A, const CVector & B);
    friend CVector operator * (const CMatrix & A, const CVector & B);
    friend CVector operator * (const CVector & A, const CMatrix & B);
    friend CMatrix dyad (const CVector & A, const CVector & B);

    CVector operator - () const;
    CVector & operator *= (complex a);
    CVector & operator /= (complex a);
    CVector & operator += (const CVector & A);
    CVector & operator -= (const CVector & A);
    CVector adjoint() const;

 private:
 
    void init(int nl, int nh);
    void load(const CVector & A);
    void loadI(const IVector & A);
    void loadD(const DVector & A);
    void destroy();

};

std::ostream & operator << (std::ostream & os, const CVector & A);

int operator == (const CVector & A, const CVector & B);
int operator != (const CVector & A, const CVector & B);
int eqacc (const CVector & A, const CVector & B, double acc);
    
CVector operator + (const CVector & A, const CVector & B);
CVector operator - (const CVector & A, const CVector & B);
CVector operator * (complex a, const CVector & B);
CVector operator * (const CVector & A, complex b);
CVector operator / (const CVector & A, complex b);
complex  operator * (const CVector & A, const CVector & B);
complex  norm(const CVector & A, const CVector & B);

  
class CMatrix{
 
 friend class CVector;

 protected:

     int ml, mh, nl, nh;
     complex * * ptr;

 public:
    
    CMatrix();
    CMatrix(int ml, int mh, int nl, int nh);
    CMatrix(const CMatrix & A);
    CMatrix(const DMatrix & A);
    CMatrix(const IMatrix & A);
    CMatrix(int nl, int nh, complex a);

    ~CMatrix();
 
    CMatrix & operator = (const CMatrix & A);
    inline cprt  &  operator [](int i){return ptr[i];}
    int rl() const;
    int rh() const;
    int cl() const;
    int ch() const;  

    void fillRow(int m, int n, ...);
        /*  fills row m with n elements   */
    void readRow(int n, const CVector & A);
    void readColumn(int n, const CVector & A);
    void zero();
  
    friend std::ostream & operator << (std::ostream & os, const CMatrix & A);

    friend int operator == (const CMatrix & A, const CMatrix & B);
    friend int operator != (const CMatrix & A, const CMatrix & B);
    friend int eqacc (const CMatrix & A, const CMatrix & B, double acc);
        /*  tests equality  to accuracy acc    */
    
    friend CMatrix operator + (const CMatrix & A, const CMatrix & B);
    friend CMatrix operator - (const CMatrix & A, const CMatrix & B);
    friend CMatrix operator * (complex a, const CMatrix & B);
    friend CMatrix operator * (const CMatrix & A, complex b);
    friend CMatrix operator / (const CMatrix & A, complex b);
    friend CVector operator * (const CMatrix & A, const CVector & B);
    friend CVector operator * (const CVector & A, const CMatrix & B);
    friend CMatrix operator * (const CMatrix & A, const CMatrix & B);
    friend CMatrix dyad (const CVector & A, const CVector & B);

    CVector row(int i) const;
    CVector column(int j) const;

    CMatrix operator - () const;
    CMatrix & operator *= (complex a);
    CMatrix & operator /= (complex a);
    CMatrix & operator += (const CMatrix & A);
    CMatrix & operator -= (const CMatrix & A);

    CMatrix transpose() const;
    CMatrix adjoint() const;
    complex trace() const;
    complex det() const;
    double logdet() const;
    CMatrix inverse() const;

    friend CMatrix iproduct(const CMatrix & A, const CMatrix & B);
       /*  iproduct returns   A^{-1} * B   */

    void rowint(int i, int j);
    void colint(int i, int j);

 private:

    void init(int ml, int mh, int nl, int nh);
    void load(const CMatrix & A);
    void loadI(const IMatrix & A);
    void destroy();
    double maxRow(int i) const;
    double maxCol(int i) const;
    CMatrix LUdecomp(int * index, complex & d) const;
    void LUbacksub(CMatrix & B, int * index);

 public:
     /*  helper function for diagonalize() below  */

     inline void ROTATE(complex s,complex  ta,
                                int i,int j,int k,int l);
     inline void ROTATE2(complex s,complex  ta,
                                int i,int j,int k,int l);
     inline void ROTATE3(complex s,complex  ta,
                                int i,int j,int k,int l);
 
};
    
    
std::ostream & operator << (std::ostream & os, const CMatrix & A);

int operator == (const CMatrix & A, const CMatrix & B);
int operator != (const CMatrix & A, const CMatrix & B);
int eqacc (const CMatrix & A, const CMatrix & B, double acc);
    
CMatrix operator + (const CMatrix & A, const CMatrix & B);
CMatrix operator - (const CMatrix & A, const CMatrix & B);
CMatrix operator * (complex a, const CMatrix & B);
CMatrix operator * (const CMatrix & A, complex b);
CMatrix operator / (const CMatrix & A, complex b);
CVector operator * (const CMatrix & A, const CVector & B);
CVector operator * (const CVector & A, const CMatrix & B);
CMatrix operator * (const CMatrix & A, const CMatrix & B);
CMatrix dyad (const CVector & A, const CVector & B);
    
CMatrix iproduct(const CMatrix & A, const CMatrix & B);


/*  Hermitian matrix   */
class HMatrix: public CMatrix {
  public:
     HMatrix(): CMatrix() {} 
     HMatrix(int xml, int xmh,int xnl, int xnh): 
              CMatrix(xml,xmh,xnl,xnh) {}
     HMatrix(const CMatrix & A): CMatrix(A){}
     HMatrix(int xnl, int xnh, complex a): CMatrix(xnl,xnh,a) {}

     CMatrix diagonalize(DVector & lambda, double acc) const;
        /*  returns eigenvectors  V[i][j] = i th component of j th vector 
             the eigenvalues are returned in lambda[i], which must be 
	       preassigned with the correct size ; 
                 acc is the accuracy with which the diagonalization is done */
  private:
     void check() const;
};



#endif      /*  CCLASSES_H     */
