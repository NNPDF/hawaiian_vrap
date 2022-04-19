
#ifndef ICLASSES_H
#define ICLASSES_H

/*   classes for int vectors and matrices   */

#include "thbasics.h"
#include <iosfwd>

typedef  int*   iprt;  

class IMatrix;

class IVector{
 
 friend class IMatrix;
 friend class DVector; 
 friend class CVector;

 protected:

     int nl, nh;
     int * ptr;

 public:
    
    IVector();
    IVector(int nl, int nh);
    IVector(const IVector & A);

    ~IVector();
 
    IVector & operator = (const IVector & A);

    inline int &  operator [](int i){  return ptr[i]; }
    int l() const;
    int h() const;

    void fill(int n, ...);
    void zero();
  
    friend std::ostream & operator << (std::ostream & os, const IVector & A);

    friend int operator == (const IVector & A, const IVector & B);
    friend int operator != (const IVector & A, const IVector & B);
    
    friend IVector operator + (const IVector & A, const IVector & B);
    friend IVector operator - (const IVector & A, const IVector & B);
    friend IVector operator * (int a, const IVector & B);
    friend IVector operator * (const IVector & A, int b);
    friend int  operator * (const IVector & A, const IVector & B);
    friend IVector operator * (const IMatrix & A, const IVector & B);
    friend IVector operator * (const IVector & A, const IMatrix & B);
    friend IMatrix dyad (const IVector & A, const IVector & B);

    IVector operator - () const;
    IVector & operator *= (int a);
    IVector & operator += (const IVector & A);
    IVector & operator -= (const IVector & A);

 private:
 
    void init(int nl, int nh);
    void load(const IVector & A);
    void destroy();

};

std::ostream & operator << (std::ostream & os, const IVector & A);

int operator == (const IVector & A, const IVector & B);
int operator != (const IVector & A, const IVector & B);
    
IVector operator + (const IVector & A, const IVector & B);
IVector operator - (const IVector & A, const IVector & B);
IVector operator * (int a, const IVector & B);
IVector operator * (const IVector & A, int b);
int  operator * (const IVector & A, const IVector & B);


  
class IMatrix{
 
 friend class IVector;
 friend class DMatrix;
 friend class CMatrix;

 protected:

     int ml, mh, nl, nh;
     int * * ptr;

 public:
    
    IMatrix();
    IMatrix(int ml, int mh, int nl, int nh);
    IMatrix(const IMatrix & A);
    IMatrix(int nl, int nh, int a);

    ~IMatrix();
 
    IMatrix & operator = (const IMatrix & A);
    inline iprt  &  operator [](int i){return ptr[i];}
    int rl() const;
    int rh() const;
    int cl() const;
    int ch() const;  

    void fillRow(int m, int n, ...);
        /*  fills row m with n elements   */
    void readRow(int n, const IVector & A);
    void readColumn(int n, const IVector & A);
    void zero();
  
    friend std::ostream & operator << (std::ostream & os, const IMatrix & A);

    friend int operator == (const IMatrix & A, const IMatrix & B);
    friend int operator != (const IMatrix & A, const IMatrix & B);
    
    friend IMatrix operator + (const IMatrix & A, const IMatrix & B);
    friend IMatrix operator - (const IMatrix & A, const IMatrix & B);
    friend IMatrix operator * (int a, const IMatrix & B);
    friend IMatrix operator * (const IMatrix & A, int b);
    friend IVector operator * (const IMatrix & A, const IVector & B);
    friend IVector operator * (const IVector & A, const IMatrix & B);
    friend IMatrix operator * (const IMatrix & A, const IMatrix & B);
    friend IMatrix dyad (const IVector & A, const IVector & B);

    IVector row(int i) const;
    IVector column(int j) const;

    IMatrix  operator - () const;
    IMatrix & operator *= (int a);
    IMatrix & operator += (const IMatrix & A);
    IMatrix & operator -= (const IMatrix & A);

    IMatrix transpose() const;
    int trace() const;

    void rowint(int i, int j);
    void colint(int i, int j);

 private:

    void init(int ml, int mh, int nl, int nh);
    void load(const IMatrix & A);
    void destroy();
};
    
    
std::ostream & operator << (std::ostream & os, const IMatrix & A);

int operator == (const IMatrix & A, const IMatrix & B);
int operator != (const IMatrix & A, const IMatrix & B);
    
IMatrix operator + (const IMatrix & A, const IMatrix & B);
IMatrix operator - (const IMatrix & A, const IMatrix & B);
IMatrix operator * (int a, const IMatrix & B);
IMatrix operator * (const IMatrix & A, int b);
IVector operator * (const IMatrix & A, const IVector & B);
IVector operator * (const IVector & A, const IMatrix & B);
IMatrix operator * (const IMatrix & A, const IMatrix & B);
IMatrix dyad (const IVector & A, const IVector & B);
    

#endif      /*  ICLASSES_H     */
