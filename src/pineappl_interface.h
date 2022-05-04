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
     */
  public:
    CheffPanopoulos();
    void create_grid(int max_orders, double q2);
    void fill_grid(int order, LuminosityFunction lumi_function, double x1,
                   double x2, double weight);
    void set_prefactor(double);
    void enable(bool state);
    void save();
    double vegas_wgt = 1.0;

  private:
    double constant_q2, prefactor = 1.0;
    // This flag is used to disable the pineappl filling until the last iteration of a Vegas call
    std::vector<LuminosityFunction> luminosities;
    bool is_enabled = true;
    pineappl_grid *grid;
};
} // namespace pinerap
