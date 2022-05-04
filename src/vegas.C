
#include "NAClasses.h"
#include "random.h"

using namespace std;

      /* It is necessary to supply a random number generator for vegas.
           This generator should follow our conventions that ran() gives
          the next random number and ran(n) initializes the sequence 
            when n is a negative integer.  Here it the default:  */

inline double ran(int n=0) {return ran1(n);}


/* Here are the member functions of the VegasGrid class  */

VegasGrid::VegasGrid(DVector Xl, DVector Xh, int ncall, int ranseed):
  xl(Xl), xh(Xh), dx(Xh-Xl),xi(Xl.l(),Xl.h(),0,VNDMX){	
          mds = 1;  nprt = -1;
	  nl = xl.l();  
	  nh = xl.h();
	  ndim = 1+ nh-nl;
          ran(-ranseed);
          reset(ncall);
}

void VegasGrid::reset(int ncall){
          it = 0;
          int i,j,k,n;
	  nd = VNDMX;
	  ng = 1;
          for (j=nl; j<=nh; j++)  xi[j][1] = 1.0;
          if (mds) {
              ng = pow(ncall/2.0 + 0.25, 1.0/ndim);
	      mds = 1;
              if ((2*ng-VNDMX) >= 0 ) {
                 mds = -1;
                 n = ng/VNDMX + 1;
                 nd = ng/n;
                 ng = n*nd;
              }
	   }
	   for (k=1, i=1; i <= ndim; i++) k *= ng;
           npg = MAX(ncall/k,2);
           calls = npg * k;
           dxg = 1.0/ng;
           dv2g = 1.0/(npg-1.0);
           xnd = nd;
           dxg = xnd * dxg;
           xjac = 1.0/calls;
           for (j=nl; j<= nh; j++)   xjac *= dx[j];
           DVector r(1,nd);
           for (j=1; j<= nd; j++) r[j] = 1.0;
           for (j=nl; j<= nh; j++) rebin(1.0/xnd,r,j);
	   si = schi = swgt = 0.0;
} 

VegasGrid::VegasGrid(const VegasGrid & A):
          xl(A.xl),xh(A.xh),dx(A.dx),xi(A.xi){
    nl = A.nl;
    nh = A.nh;
    ndim = A.ndim;
    nd = A.nd;
    xnd = A.xnd;
    npg = A.npg;
    ng = A.ng;
    calls = A.calls;
    xjac = A.xjac;
    dxg = A.dxg;
    dv2g = A.dv2g;
    si = A.si;
    swgt = A.swgt;
    schi = A.schi;
    it = A.it;
    mds = A.mds;
    nprt = A.nprt;
}

VegasGrid & VegasGrid::operator = (const VegasGrid & A){
    xl = A.xl;
    xh = A.xh;
    dx = A.dx;
    xi = A.xi;
    nl = A.nl;
    nh = A.nh;
    ndim = A.ndim;
    nd = A.nd;
    xnd = A.xnd;
    npg = A.npg;
    ng = A.ng;
    calls = A.calls;
    xjac = A.xjac;
    dxg = A.dxg;
    dv2g = A.dv2g;
    si = A.si;
    swgt = A.swgt;
    schi = A.schi;
    it = A.it;
    mds = A.mds;
    nprt = A.nprt;
    return * this;
}



/* this function clears the answers accumulated by a VegasGrid  */
    
void VegasGrid::clear() {
     si = 0.0;
     swgt = 0.0;
     schi = 0.0;
     it = 0;
   }

/* this function changes the number of calls associated with a VegasGrid */

void VegasGrid::morecalls(int nc){
    int i,k;
    cout << ng << " " << npg << "  " << calls << endl;
    int nx = nc/calls;
    int ny  = pow(nx + 0.25, 1.0/ndim);
    ng *= ny;
    for (k=1, i=1; i <= ndim; i++) k *= ny;
    int nz = nx/k;
    npg *= nz;
    calls *= k*nz;
    cout << ng << " " << npg << "  " << calls << endl;
    xjac *= 1.0/(k*nz);
    dxg = xnd/ng;
    dv2g = 1.0/(npg-1);
}
     

/*    this function modifies the grid for random sampling   */

VegasGrid VegasGrid::sampling(int ncalls) const {
    VegasGrid A(*this);
    A.npg = ncalls;
    A.ng = 1;
    A.calls = ncalls;
    A.xjac *= calls/A.calls;
    A.dxg = nd;
    A.dv2g = 1.0/(ncalls-1);
    A.si = 0.0;
    A.swgt = 0.0;
    A.schi = 0.0;
    A.it = 0;
    A.mds = 0;
    return A;
}              

void VegasGrid::refine(DMatrix d) {
  //   cout << d ;                // print grid, for debugging purposes
       int i,j;
       double rc,xo,xn;
       DVector dt(1,nd), r(1,nd);
       for (j=nl; j <= nh; j++) {
	 xo = d[j][1];
	 xn = d[j][2];
	 d[j][1] = (xo + xn)/2.0;
	 dt[j] = d[j][1];
	 for (i=2; i < nd; i++) {
            rc = xo+xn;
	    xo = xn;
	    xn = d[j][i+1];
	    d[j][i] = (rc+xn)/3.0;
	    dt[j] += d[j][i];
         }
	 d[j][nd] = (xo+xn)/2.0;
         dt[j] += d[j][nd];
       }
       for (j=nl; j<= nh; j++) {
	 rc = 0.0;
	 for (i=1; i<=nd; i++) {
	    if (d[j][i] < TINY) d[j][i] = TINY;
            r[i] = pow((1.0-d[j][i]/dt[j])/(log(dt[j])-log (d[j][i])),
                                              VEGASALPH);
	    rc += r[i];
	 }
         rebin(rc/xnd,r, j);
       }
    }

void VegasGrid::rebin(double rc, DVector & r, int j) {
     int i,k = 0;
     double dr =0.0, xn = 0.0, xo;
     DVector xinew(1,nd);
   
     for (i=1;i<nd;i++) {
       while (rc > dr) {
         dr += r[++k];
         xo = xn;
         xn = xi[j][k];
       }
       dr -= rc;
       xinew[i] = xn - (xn-xo)*dr/r[k];
     }
     for (i=1; i<nd; i++) xi[j][i] = xinew[i];
     xi[j][nd] = 1.0;
   } 

DVector VegasGrid::choose(double & xvol,VegasCounter & VC) { 
     xvol = xjac;
     DVector x(nl,nh);
     int j;
     for (j=nl; j <= nh; j++){
        double xn = (VC.kg[j] - ran()) * dxg + 1.0;
        int ixn = xn;
        int iaj  = MAX(MIN(ixn,nd),1);
        double xo = xi[j][iaj] - xi[j][iaj-1];
        double rc = xi[j][iaj-1] + (xn - iaj)*xo;
        x[j] = xl[j] + rc * dx[j];
        xvol *= xnd * xo; 
     }
     VC.move();
     return x;
   }
  
DVector VegasGrid::choose(double & xvol,VegasCounter & VC,IVector & ia) { 
     xvol = xjac;
     DVector x(nl,nh);
     int j;
     for (j=nl; j <= nh; j++){
        double xn = (VC.kg[j] - ran()) * dxg + 1.0;
        int ixn = xn;
        int iaj  = MAX(MIN(ixn,nd),1);
        double xo = xi[j][iaj] - xi[j][iaj-1];
        double rc = xi[j][iaj-1] + (xn - iaj)*xo;
        x[j] = xl[j] + rc * dx[j];
        xvol *= xnd * xo; 
        ia[j] = iaj; 
     }
     VC.move();
     return x;
}


void VegasGrid::record(double fb, double f2b,double & ti, double & tsi){
        double f2bb = sqrt(f2b*npg);
        f2bb = (f2bb-fb)*(f2bb+fb);
        if (f2bb <= 0.0) f2b = TINY;
        ti += fb;
        tsi += f2bb;
}

void VegasGrid::finish(double ti, double & tsi, double & integral, 
        double & sd, double & chi2){
     tsi *= dv2g;
     double wgt = 1.0/tsi;
     si += wgt*ti;
     schi += wgt*ti*ti;
     swgt += wgt;
     integral = si/swgt;
     chi2 = (schi-si*integral)/(it-0.9999);
     if (chi2 < 0.0) chi2 = 0.0;
     sd = sqrt(1.0/swgt);
}
    
     
void VegasGrid::integralinfo(int itt,double ti,double tsi,double integral,
                   double sd,double chi2){       
      cout << "iteration no. " << itt << " : integral = " << ti << " +/- " <<
                  sqrt(tsi) << endl;
      cout << "all iterations   it = "<<it << " : integral = "
        << integral << " +/- " << sd << "  chi2/IT n = " << chi2 
                                << endl<< endl;
   }
   
           


/* Here are the member functions of the VegasCounter class  */


VegasCounter::VegasCounter(const VegasGrid & grid): kg(grid.nl,grid.nh){
    step = 0;
    done = 0;
    nl = grid.nl;
    nh = grid.nh;
    ng = grid.ng;
    npg = grid.npg;
    kct = 1;
    for ( int i = nl; i <= nh; i++)  kg[i] = 1;
}

void VegasCounter::reset(){
    step = 0;
    done = 0;
    kct = 1;
    for ( int i = nl; i <= nh; i++)  kg[i] = 1;
}

void VegasCounter::resetstep(){
    step = 0;
    kct = 1;
}

void VegasCounter::move(){  
    kct += 1;
    if (kct <= npg) return;
    kct = 1;
    step = 1;
    int j;
    for (j= nh; j >= nl; j--) {
      kg[j] += 1; 
      if (kg[j] <= ng ) break;
      kg[j] = 1;
    }
    if (j < nl) done = 1;
}
  


         /*  adaptive Monte Carlo integration of the volume under 
                surface(V)  in the area specified by the vectors xl, xu.
            The function returns the values of the integral, an estimate
          of the standard deviation (in sd) and a chi^2 per degree of 
           freedom (in chi2).  */

double Surface::vegas(VegasGrid & grid, int its, double & sd, double & chi2) {
  int itt;
  int nl = grid.nl; int nh = grid.nh; int nd = grid.nd;
  int j,k;
  double f, fb, f2, f2b, ti, tsi, integral; 
  double wgt = 0.0;
  DVector x(nl,nh);
  IVector ia(nl,nh);
  DMatrix d(nl,nh,1,nd);
  VegasCounter VC(grid);
  for (itt = 1; itt <= its; itt++) { 
     grid.it++;
     ti = 0.0;  tsi = 0.0;
     VC.reset();
     while (VC.done == 0){
        VC.resetstep();
        fb = 0.0;  f2b = 0.0;
        while (VC.step == 0){
           x = grid.choose(wgt,VC,ia);
 	   f = wgt*surface(x);
	   f2 = f*f;
	   fb += f;
	   f2b += f2;
           if (grid.mds >= 0){
             for (j=nl; j<= nh; j++) d[j][ia[j]] += f2;
           }
        }
        grid.record(fb,f2b,ti,tsi);
        if (grid.mds < 0) {
           for (j=nl; j<= nh; j++) d[j][ia[j]] += f2b;
        }
     }
     grid.finish(ti,tsi,integral,sd,chi2);
     if (grid.nprt > 0) grid.integralinfo(itt,ti,tsi,integral,sd,chi2);
     grid.refine(d);
  }
  return integral;
}


double Surface::vegasint(VegasGrid & grid, int its,
           double & sd, double & chi2) {
  int itt;
  double f, fb, f2, f2b, integral,ti,tsi; 
  double wgt = 0.0;
  DVector x(grid.nl,grid.nh);
  VegasCounter VC(grid);
  for (itt = 1; itt <= its; itt++) { 
    // Enable pineappl filling only for the last iteration
    if (itt == its) vegas_piner->enable(true);
    grid.it++;
    ti = 0.0; tsi = 0.0;
    VC.reset();
    while (VC.done == 0){
        VC.resetstep();
        fb = f2b = 0.0;
        while (VC.step == 0){
            x = grid.choose(wgt,VC);
            // I no longer feel shame
            vegas_piner->vegas_wgt = wgt;
            f = wgt*surface(x);
            f2 = f*f;
            fb += f;
            f2b += f2;
        }
        grid.record(fb,f2b,ti,tsi);
    }
  /*   compute the final results for this iteration  */
    grid.finish(ti,tsi,integral,sd,chi2);
    if (grid.nprt > 0)  grid.integralinfo(itt,ti,tsi,integral,sd,chi2);
  }
  return integral;
}

double Surface::sampleint(VegasGrid & grid, int its,
           double & sd, double & chi2) {
  if (grid.mds != 0 || grid.ng != 1) 
             therror("  the VegasGrid is not set up for random sampling ");
  int n;
  double f, fb, f2, f2b, integral; 
  double wgt = 0.0;
  DVector x(grid.nl,grid.nh);
  VegasCounter VC(grid);
  for (n = 1; n<= its; n++){
     grid.it++;
     VC.reset();
     fb = f2b = 0.0;
     while (VC.done == 0){
       x = grid.choose(wgt,VC); 
       f = wgt*surface(x);
       f2 = f*f;
       fb += f;
       f2b += f2;
     }
  /*   compute the final results   */
     f2b = sqrt(f2b*grid.calls);
     f2b = (f2b - fb)*(f2b + fb);
     grid.finish(fb,f2b,integral,sd,chi2);
     if (grid.nprt > 0)  grid.integralinfo(n,fb,f2b,integral,sd,chi2);
  }
  return integral;
}

DVector Surface::sample(VegasGrid & grid, double & weight){
  if (grid.mds != 0 || grid.ng != 1) 
             therror("  the VegasGrid is not set up for random sampling ");
  int n;
  DVector x(grid.nl,grid.nh);
  double f, fb, f2, f2b, integral; 
  double wgt = 0.0;
  VegasCounter VC(grid);
  grid.it++;
  fb = f2b = 0.0;
  x = grid.choose(wgt,VC); 
  f = wgt*surface(x);
  f2 = f*f;
  fb += f;
  f2b += f2;
  weight = f;
  return x;
}
