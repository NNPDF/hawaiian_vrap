#pragma once

#include "LHApdf.h"
#include "pdf.h"
#include "Vlumifns_LHApdf.h"
#include <pineappl_capi.h>
#include <vector>
#include<bits/stdc++.h>

namespace pinerap {

typedef double (*LuminosityFunction)(pdfArray const &, pdfArray const &,
                                     collider);

void reconstruct_lumi(LuminosityFunction lumi, collider c,
                      std::vector<int32_t> &pdg_ids,
                      std::vector<double> &factors);

class CheffPanopoulos {
    /*
     * Flavour up vrap with a bit of pineappl
     *
     * The grid construction works as follows:
     *  vrap is able to compute one pair of Y and M at once
     *  due to different tricks we use, we can equally store only one pair of Y of M at once
     *  since we cannot fill 2D anyway we are going to create one grid per pair with dummy binning (given next_grid_index)
     *  
     *  Every time that `create_grid` is called, `next_grid_index` will be checked.
     *  If next_grid_index > 1 then the current `grid` will be merged into `mother_grid` and a new grid will be created
     *  for vrap to fill.
     *
     *  At the end of the calculation the `mother_grid` will be saved.
     *
     *
     */
  public:
    CheffPanopoulos();
    void create_grid(int max_orders, double q2, collider coll);
    void fill_grid(int order, LuminosityFunction lumi_function, double x1,
                   double x2, double weight);
    void set_prefactor(const double);
    void enable(const bool state);
    void rebin(const std::vector<std::pair<double, double>>);
    void save();
    double vegas_wgt = 1.0;

  private:
    double constant_q2, prefactor = 1.0;
    int next_grid_index = 0;
    // This flag is used to disable the pineappl filling until the last iteration of a Vegas call
    std::vector<LuminosityFunction> luminosities;
    bool is_enabled = true;
    pineappl_grid *grid, *mother_grid; 
};

typedef double (unlog_f1)(double, double, double, double); // z, nf, muF/Q, muR/Q
typedef double (unlog_f0)(double, double, double); // nf, muF/Q, muR/Q
typedef double (unlog_f2)(double, double); // z, mu/Q
void unlog_muFmuR0(double Nf, unlog_f0 fun, double* logterms);
void unlog_muFmuR(double, double, unlog_f1, double, double*);
void unlog_muF(double, unlog_f2, double, double*);


typedef double (unlog_r1)(double, double, double, double, double); // ys, z, nf, muF/Q, muR/Q
typedef double (unlog_nf)(double, double, double, double); // ys, z, Nf, mu/Q
typedef double (unlog_r2)(double, double, double); // ys, z, mu/Q
void r_unlog_muFmuR(double, double, double, unlog_r1, double, double*);
void r_unlog_muF(double, double, double, unlog_nf, double, double*);
void r_unlog_muF2(double, double, unlog_r2, double, double*);
} // namespace pinerap
