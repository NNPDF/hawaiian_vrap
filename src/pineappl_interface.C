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

    // TODO: add missing partons
    std::vector<std::pair<double pdfArray::*, int32_t>> partons = {
        {&pdfArray::tbar, -6},  {&pdfArray::bbar, -5}, {&pdfArray::cbar, -4},
        {&pdfArray::sbar, -3},  {&pdfArray::ubar, -2}, {&pdfArray::dbar, -1},
        {&pdfArray::gluon, 21}, {&pdfArray::d, 1},     {&pdfArray::u, 2},
        {&pdfArray::s, 3},      {&pdfArray::c, 4},     {&pdfArray::b, 5},
        {&pdfArray::t, 6},      {&pdfArray::x, 9999}, // TODO: what is `x`?
    };

    // loop over partons of the first initial state
    for (auto const a : partons) {
        // loop over partons of the second initial state
        for (auto const b : partons) {
            // TODO: check that this zero-initializes
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
    luminosities.push_back(&qg_lumi);
    luminosities.push_back(&gq_lumi);
}

/*
 * Initializes the pineappl grid with all the necessary parameters
 */
void CheffPanopoulos::create_grid(int max_orders, double q2) {
    // Check whether this is a new grid that needs to be absorbed by mother_grid
    if (next_grid_index > 1) pineappl_grid_merge_and_delete(mother_grid, grid);

    auto *lumi = pineappl_lumi_new();

    // At LO we only have one lumi channel
    std::vector<int32_t> pdg_ids;
    std::vector<double> factors;
    // we only care about collider piso for now (=3) take it from the runcard in
    // the future
    for (auto lumi_function: luminosities) {
        reconstruct_lumi(lumi_function, piso, pdg_ids, factors);
        pineappl_lumi_add(lumi, factors.size(), pdg_ids.data(), factors.data());
        pdg_ids.clear();
        factors.clear();
    }


    // Only LO for now
    // ------------------- (as, a, muR, muF)
    std::vector<uint32_t> orders{0, 2, 0, 0};

    if (order_flag > 0) { // Add NLO orders
        orders.insert(orders.end(), {1, 2, 0, 0});
        orders.insert(orders.end(), {1, 2, 1, 0});
        orders.insert(orders.end(), {1, 2, 0, 1});
    }
    if (order_flag > 1) { // Add NLO orders
        orders.insert(orders.end(), {2, 2, 0, 0});
        orders.insert(orders.end(), {2, 2, 1, 0});
        orders.insert(orders.end(), {2, 2, 0, 1});
        orders.insert(orders.end(), {2, 2, 1, 1});
        orders.insert(orders.end(), {2, 2, 1, 1});
        orders.insert(orders.end(), {2, 2, 2, 0});
        orders.insert(orders.end(), {2, 2, 0, 2});
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

    pineappl_keyval_delete(keyval);
    pineappl_lumi_delete(lumi);
}

/*
 * Given an order:
 *  0 - LO
 *  1 - NLO
 *  2 - NNLO
 * and a luminosity channel
 *  TBD
 * fills the grid
 */
void CheffPanopoulos::fill_grid(int order, LuminosityFunction lumi_function, double x1,
                                double x2, double weight) {
    if (is_enabled) {
        auto lumi_index = std::find(luminosities.begin(), luminosities.end(), lumi_function);
        int lumi_channel = lumi_index - luminosities.begin();
        double res = weight*prefactor*vegas_wgt;
        if (order > 0) res /= PI;
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

    std::cout << "Rebinning\n";
    const int nbins = qy_bins.size();
    const double normalization[nbins] = {1.0};
    double limits[2*nbins*2];
    int i = 0;
    for(auto const qy: qy_bins) {
        const double Q = qy.first;
        const double y = qy.second;
        limits[i++] = Q;
        limits[i++] = Q;
        limits[i++] = y;
        limits[i++] = y;
    }
    pineappl_grid_set_remapper(mother_grid, 2, normalization, limits);
}

void CheffPanopoulos::set_prefactor(const double val){
    prefactor = val;
}

void CheffPanopoulos::enable(const bool state) {
    is_enabled = state;
}
