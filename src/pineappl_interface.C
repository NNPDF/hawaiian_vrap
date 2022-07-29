#include "pineappl_interface.h"
#include "QCDbasics.h"

#include <cstddef>
#include <iostream>
#include <utility>

using namespace pinerap;

void pinerap::reconstruct_lumi(LuminosityFunction lumi, collider c,
                      std::vector<int32_t> &pdg_ids,
                      std::vector<double> &factors) {
    pdfArray f1;
    pdfArray f2;
    pdfArray sampler;

    std::vector<std::pair<double pdfArray::*, int32_t>> partons = {
        {&pdfArray::tbar, -6},  {&pdfArray::bbar, -5}, {&pdfArray::cbar, -4},
        {&pdfArray::sbar, -3},  {&pdfArray::ubar, -2}, {&pdfArray::dbar, -1},
        {&pdfArray::gluon, 21}, {&pdfArray::d, 1},     {&pdfArray::u, 2},
        {&pdfArray::s, 3},      {&pdfArray::c, 4},     {&pdfArray::b, 5},
        {&pdfArray::t, 6},      {&pdfArray::x, 9999},
    };

    // Create a sampler to automatically skip partons that the PDF does not include
    LHAComputePdf(0.1, 10.0, sampler);

    // loop over partons of the first initial state
    for (auto const a : partons) {
        if (sampler.*a.first == 0.0) continue;

        // loop over partons of the second initial state
        for (auto const b : partons) {
            if (sampler.*b.first == 0.0) continue;

            f1 = pdfArray();
            f2 = pdfArray();

            f1.*a.first = 1.0;
            f2.*b.first = 1.0;

            // process is always DY (1)
            double const result = lumi(f1, f2, c);

            // skip zero partonic combinations
            if (result == 0.0) {
                continue;
            }

            pdg_ids.push_back(a.second);
            pdg_ids.push_back(b.second);
            factors.push_back(result);
        }
    }
}

/* Creates a vector with all possible luminosities
 */
CheffPanopoulos::CheffPanopoulos() {
    luminosities.push_back(&qqbar_lumi_dy);
    // appear at NLO
    luminosities.push_back(&qg_lumi);
    luminosities.push_back(&gq_lumi);
    // appear at NNLO
    luminosities.push_back(&qqbar_ax_lumi);
    luminosities.push_back(&qqbar_BC_lumi);
    luminosities.push_back(&qqbar_lumi_g);
    luminosities.push_back(&qq_11_lumi);
    luminosities.push_back(&qq_12_lumi);
    luminosities.push_back(&qq_12_lumi);
    luminosities.push_back(&qq_12_ax_lumi);
    luminosities.push_back(&qq_22_lumi);
    luminosities.push_back(&qq_CE1_lumi);
    luminosities.push_back(&qq_CE2_lumi);
    luminosities.push_back(&qq_CF_lumi);
    luminosities.push_back(&gg_lumi);
}

/*
 * Initializes the pineappl grid with all the necessary parameters
 */
void CheffPanopoulos::create_grid(int max_orders, double q2, collider coll) {
    // Check whether this is a new grid that needs to be absorbed by mother_grid
    if (next_grid_index > 1) pineappl_grid_merge_and_delete(mother_grid, grid);

    auto *lumi = pineappl_lumi_new();

    // At LO we only have one lumi channel
    std::vector<int32_t> pdg_ids;
    std::vector<double> factors;
    // Take the collider type from the runcard
    for (auto lumi_function: luminosities) {
        reconstruct_lumi(lumi_function, coll, pdg_ids, factors);
        if (!pdg_ids.empty()) { 
            pineappl_lumi_add(lumi, factors.size(), pdg_ids.data(), factors.data());
        }
        pdg_ids.clear();
        factors.clear();
    }


    // Leading Order
    // ------------------- (as, a, muR, muF)
    std::vector<uint32_t> orders{0, 2, 0, 0};

    if (order_flag > 0) { // Add NLO orders
        orders.insert(orders.end(), {1, 2, 0, 0}); // nlo 
        orders.insert(orders.end(), {1, 2, 1, 0}); // muR
        orders.insert(orders.end(), {1, 2, 0, 1}); // muF
    }
    if (order_flag > 1) { // Add NNLO orders
        orders.insert(orders.end(), {2, 2, 0, 0}); // nnlo
        orders.insert(orders.end(), {2, 2, 1, 0}); // muR
        orders.insert(orders.end(), {2, 2, 0, 1}); // muF
        orders.insert(orders.end(), {2, 2, 1, 1}); // muR x muF
        orders.insert(orders.end(), {2, 2, 2, 0}); // muR^2
        orders.insert(orders.end(), {2, 2, 0, 2}); // muF^2 
    }
    int how_many_orders = orders.size() / 4;

    double grid_limit = 1.0*next_grid_index;
    double bins[] = {grid_limit, grid_limit+1.0};
    int how_many_bins = 1;

    auto *keyval = pineappl_keyval_new();
    pineappl_keyval_set_double(keyval, "q2_min", 2.5);
    grid = pineappl_grid_new(lumi, how_many_orders, orders.data(), how_many_bins, bins,
                             keyval);

    if (next_grid_index == 0) mother_grid = grid;
    next_grid_index += 1;

    constant_q2 = q2;
    vegas_wgt = 1.0;

    pineappl_keyval_delete(keyval);
    pineappl_lumi_delete(lumi);
}

/*
 * Given an order:
 *  0 - LO
 *  1 - NLO
 *  2 - NNLO
 * and a luminosity channel
 * fills the grid
 */
void CheffPanopoulos::fill_grid(int order, LuminosityFunction lumi_function, double x1,
                                double x2, double weight) {

    if (is_enabled) {
        auto lumi_index = std::find(luminosities.begin(), luminosities.end(), lumi_function);
        int lumi_channel = lumi_index - luminosities.begin();
        double res = weight*prefactor*vegas_wgt;
        if (order > 0) res /= PI;
        if (order > 3) res /= PI;
        pineappl_grid_fill(grid, x1, x2, constant_q2, order, next_grid_index-0.5, lumi_channel, res);
    }
}

void CheffPanopoulos::save() {
    char const *filename = "test.pineappl.lz4";
    pineappl_grid_write(mother_grid, filename);
    pineappl_grid_delete(mother_grid);
}

void CheffPanopoulos::rebin(const std::vector<std::pair<double, double>> qy_bins) {
    // Merge the last bin into mother grid
    if (next_grid_index > 1) pineappl_grid_merge_and_delete(mother_grid, grid);
    const int nbins = qy_bins.size();
    std::vector<double> normalization(nbins, 1.0);
    std::vector<double> limits;
    limits.reserve(2*nbins*2);
    for(auto const qy: qy_bins) {
        const double Q = qy.first;
        const double y = qy.second;
        limits.push_back(Q);
        limits.push_back(Q);
        limits.push_back(y);
        limits.push_back(y);
    }
    pineappl_grid_set_remapper(mother_grid, 2, normalization.data(), limits.data());
}

void CheffPanopoulos::set_prefactor(const double val){
    prefactor = val;
}

void CheffPanopoulos::enable(const bool state) {
    is_enabled = state;
}

// unlogging
void pinerap::unlog_muFmuR0(double Nf, unlog_f0 fun, double* logterms) {
    const double e = std::exp(0.5);
    const double ooe = std::exp(-0.5);
    //                           muF muQ
    const double c0 = fun(Nf, 1., 1.);
    const double c1pc4 = fun(Nf, 1., e) -c0;
    const double mc1pc4 = fun(Nf, 1., ooe) -c0;
    const double c2pc5 = fun(Nf, e, 1.) -c0;
    const double mc2pc5 = fun(Nf, ooe, 1.) -c0;
    const double c3 = fun(Nf, e, e) - c0 - c1pc4 - c2pc5;

    logterms[0] += (c1pc4-mc1pc4)/2.0;
    logterms[1] += (c2pc5-mc2pc5)/2.0;
    logterms[2] += c3;
    logterms[3] += (c1pc4+mc1pc4)/2.0;
    logterms[4] += (c2pc5+mc2pc5)/2.0;
}
void pinerap::unlog_muFmuR(double z, double y, unlog_f1 fun, double factor, double* logterms) {
    // The function will produce
    // fun(... muF, muR) = c0 + c1*log(muR) + c2*log(muF) + c3*log(muF*muR) + c4*log(muR)**2 + c5*log(muF)**2
    // We will use the following values to find all coefficients
    // log(1) = 0
    // log(e) = 1
    // log(1/e) = -1
    const double e = std::exp(0.5);
    const double ooe = std::exp(-0.5);
    //                           muF muQ
    const double c0 = fun(z, y, 1., 1.);
    const double c1pc4 = fun(z, y, 1., e) -c0;
    const double mc1pc4 = fun(z, y, 1., ooe) -c0;
    const double c2pc5 = fun(z, y, e, 1.) -c0;
    const double mc2pc5 = fun(z, y, ooe, 1.) -c0;
    const double c3 = fun(z, y, e, e) - c0 - c1pc4 - c2pc5;

    logterms[0] += factor*(c1pc4-mc1pc4)/2.0;
    logterms[1] += factor*(c2pc5-mc2pc5)/2.0;
    logterms[2] += factor*c3;
    logterms[3] += factor*(c1pc4+mc1pc4)/2.0;
    logterms[4] += factor*(c2pc5+mc2pc5)/2.0;
}
void pinerap::unlog_muF(double z, unlog_f2 fun, double factor, double* logterms) {
    // The function will produce at max:
    // fun(... muF) = c0 + c1*log(muF) + c2*log(muF)**2
    // We will use the following values to find all coefficients
    // log(1) = 0
    // log(e) = 1
    // log(1/e) = -1
    const double e = std::exp(0.5);
    const double ooe = std::exp(-0.5);
    //                           muF muQ
    const double c0 = fun(z, 1.);
    const double c1pc2 = fun(z, e) -c0;
    const double mc1pc2 = fun(z, ooe) -c0;

    logterms[1] += factor*(c1pc2-mc1pc2)/2.0;
    logterms[4] += factor*(c1pc2+mc1pc2)/2.0;
}
void pinerap::r_unlog_muFmuR(double ys, double z, double nf, unlog_r1 fun, double factor, double* logterms) {
    const double e = std::exp(0.5);
    const double ooe = std::exp(-0.5);
    //                           muF muQ
    const double c0 = fun(ys,z, nf, 1., 1.);
    const double c1pc4 = fun(ys,z, nf, 1., e) -c0;
    const double mc1pc4 = fun(ys,z, nf, 1., ooe) -c0;
    const double c2pc5 = fun(ys,z, nf, e, 1.) -c0;
    const double mc2pc5 = fun(ys,z, nf, ooe, 1.) -c0;
    const double c3 = fun(ys,z, nf, e, e) - c0 - c1pc4 - c2pc5;

    logterms[0] += factor*(c1pc4-mc1pc4)/2.0;
    logterms[1] += factor*(c2pc5-mc2pc5)/2.0;
    logterms[2] += factor*c3;
    logterms[3] += factor*(c1pc4+mc1pc4)/2.0;
    logterms[4] += factor*(c2pc5+mc2pc5)/2.0;
}
void pinerap::r_unlog_muF2(double ys, double z, unlog_r2 fun, double factor, double* logterms) {
    const double e = std::exp(0.5);
    const double ooe = std::exp(-0.5);
    const double c0 = fun(ys,z, 1.);
    const double c1pc2 = fun(ys,z, e) -c0;
    const double mc1pc2 = fun(ys,z, ooe) -c0;

    logterms[1] += factor*(c1pc2-mc1pc2);
    logterms[4] += factor*(c1pc2+mc1pc2);
}
void pinerap::r_unlog_muF(double ys, double z, double Nf, unlog_nf fun, double factor, double* logterms) {
    const double e = std::exp(0.5);
    const double ooe = std::exp(-0.5);
    const double c0 = fun(ys,z,Nf, 1.);
    const double c1pc2 = fun(ys,z,Nf, e) -c0;
    const double mc1pc2 = fun(ys,z,Nf, ooe) -c0;

    logterms[1] += factor*(c1pc2-mc1pc2)/2.0;
    logterms[4] += factor*(c1pc2+mc1pc2)/2.0;
}
