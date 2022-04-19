
#ifndef DCLASSES_H
#define DCLASSES_H

/*   classes for double vectors and matrices   */

#include "thbasics.h"
#include "IClasses.h"

typedef  double*   dprt;  

class DMatrix;

class DVector{
 
 friend class DMatrix;
 friend class CVector;

 protected:

     int nl, nh;
     double * ptr;

 public:
    
    DVector();
    DVector(int nl, int nh);
    DVector(const DVector & A);
    DVector(const IVector & A);

    ~DVector();
 
    DVector & operator = (const DVector & A);

    inline double &  operator [](int i){  return ptr[i]; }
    int l() const;
    int h() const;

    void fill(int n, ...);
    void zero();
  
    friend std::ostream & operator << (std::ostream & os, const DVector & A);

    friend int operator == (const DVector & A, const DVector & B);
    friend int operator != (const DVector & A, const DVector & B);
    friend int eqacc (const DVector & A, const DVector & B,double acc);
        /*  tests equality  to accuracy acc    */
    
    friend DVector operator + (const DVector & A, const DVector & B);
    friend DVector operator - (const DVector & A, const DVector & B);
    friend DVector operator * (double a, const DVector & B);
    friend DVector operator * (const DVector & A, double b);
    friend DVector operator / (const DVector & A, double b);
    friend double  operator * (const DVector & A, const DVector & B);
    friend DVector operator * (const DMatrix & A, const DVector & B);
    friend DVector operator * (const DVector & A, const DMatrix & B);
    friend DMatrix dyad (const DVector & A, const DVector & B);

    DVector operator - () const;
    DVector & operator *= (double a);
    DVector & operator /= (double a);
    DVector & operator += (const DVector & A);
    DVector & operator -= (const DVector & A);

 private:
 
    void init(int nl, int nh);
    void load(const DVector & A);
    void loadI(const IVector & A);
    void destroy();

};

std::ostream & operator << (std::ostream & os, const DVector & A);

int operator == (const DVector & A, const DVector & B);
int operator != (const DVector & A, const DVector & B);
int eqacc (const DVector & A, const DVector & B,double acc);
    
DVector operator + (const DVector & A, const DVector & B);
DVector operator - (const DVector & A, const DVector & B);
DVector operator * (double a, const DVector & B);
DVector operator * (const DVector & A, double b);
DVector operator / (const DVector & A, double b);
double  operator * (const DVector & A, const DVector & B);


  
class DMatrix{
 
 friend class DVector;
 friend class CMatrix;

 protected:

     int ml, mh, nl, nh;
     double * * ptr;

 public:
    
    DMatrix();
    DMatrix(int ml, int mh, int nl, int nh);
    DMatrix(const DMatrix & A);
    DMatrix(const IMatrix & A);
    DMatrix(int nl, int nh, double a);

    ~DMatrix();
 
    DMatrix & operator = (const DMatrix & A);
    inline dprt  &  operator [](int i){return ptr[i];}
    int rl() const;
    int rh() const;
    int cl() const;
    int ch() const;  

    void fillRow(int m, int n, ...);
        /*  fills row m with n elements   */
    void readRow(int n, const DVector & A);
    void readColumn(int n, const DVector & A);
    void zero();
  
    friend std::ostream & operator << (std::ostream & os, const DMatrix & A);

    friend int operator == (const DMatrix & A, const DMatrix & B);
    friend int operator != (const DMatrix & A, const DMatrix & B);
    friend int eqacc (const DMatrix & A, const DMatrix & B,double acc);
        /*  tests equality  to accuracy acc    */
    
    friend DMatrix operator + (const DMatrix & A, const DMatrix & B);
    friend DMatrix operator - (const DMatrix & A, const DMatrix & B);
    friend DMatrix operator * (double a, const DMatrix & B);
    friend DMatrix operator * (const DMatrix & A, double b);
    friend DMatrix operator / (const DMatrix & A, double b);
    friend DVector operator * (const DMatrix & A, const DVector & B);
    friend DVector operator * (const DVector & A, const DMatrix & B);
    friend DMatrix operator * (const DMatrix & A, const DMatrix & B);
    friend DMatrix dyad (const DVector & A, const DVector & B);

    DVector row(int i) const;
    DVector column(int j) const;

    DMatrix operator - () const;
    DMatrix & operator *= (double a);
    DMatrix & operator /= (double a);
    DMatrix & operator += (const DMatrix & A);
    DMatrix & operator -= (const DMatrix & A);

    DMatrix transpose() const;
    double trace() const;
    double det() const;
    double logdet() const;
    DMatrix inverse() const;

    friend DMatrix iproduct(const DMatrix & A, const DMatrix & B);
       /*  iproduct returns   A^{-1} * B   */

    void rowint(int i, int j);
    void colint(int i, int j);

 private:

    void init(int ml, int mh, int nl, int nh);
    void load(const DMatrix & A);
    void loadI(const IMatrix & A);
    void destroy();
    double maxRow(int i) const;
    double maxCol(int i) const;
    DMatrix LUdecomp(int * index, double & d) const;
    void LUbacksub(DMatrix & B, int * index);

 public:
     /*  helper function for diagonalize() below  */
    inline void ROTATE(double s,double ta,
                                int i,int j,int k,int l);
};
    
    
std::ostream & operator << (std::ostream & os, const DMatrix & A);

int operator == (const DMatrix & A, const DMatrix & B);
int operator != (const DMatrix & A, const DMatrix & B);
int eqacc (const DMatrix & A, const DMatrix & B,double acc);
    
DMatrix operator + (const DMatrix & A, const DMatrix & B);
DMatrix operator - (const DMatrix & A, const DMatrix & B);
DMatrix operator * (double a, const DMatrix & B);
DMatrix operator * (const DMatrix & A, double b);
DMatrix operator / (const DMatrix & A, double b);
DVector operator * (const DMatrix & A, const DVector & B);
DVector operator * (const DVector & A, const DMatrix & B);
DMatrix operator * (const DMatrix & A, const DMatrix & B);
DMatrix dyad (const DVector & A, const DVector & B);
    
DMatrix iproduct(const DMatrix & A, const DMatrix & B);

/*    real symmetric matrix  */
class RSMatrix: public DMatrix {
  public:
     RSMatrix(): DMatrix() {} 
     RSMatrix(int xml, int xmh,int xnl, int xnh): 
              DMatrix(xml,xmh,xnl,xnh) {}
     RSMatrix(const DMatrix & A): DMatrix(A){}
     RSMatrix(int xnl, int xnh, double a): DMatrix(xnl,xnh,a) {}

     DMatrix diagonalize(DVector & lambda, double acc) const;
        /*  returns eigenvectors  V[i][j] = i th component of j th vector 
             the eigenvalues are returned in lambda[i], which must be 
	       preassigned with the correct size ; 
                 acc is the accuracy with which the diagonalization is done */
  protected:
     void check() const;
};


#endif      /*  DCLASSES_H     */
