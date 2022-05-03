#include "pineappl_interface.h"

#include <cstddef>
#include <iostream>
#include <utility>

using namespace pinerap;

void reconstruct_lumi(
    double (*lumi) (pdfArray const&, pdfArray const&, collider),
    collider c,
    std::vector<int32_t>& pdg_ids,
    std::vector<double>& factors
) {
    pdfArray f1;
    pdfArray f2;

    // TODO: add missing partons
    std::vector<std::pair<double pdfArray::*, int32_t>> partons = {
        { &pdfArray::tbar, -6 },
        { &pdfArray::bbar, -5 },
        { &pdfArray::cbar, -4 },
        { &pdfArray::sbar, -3 },
        { &pdfArray::ubar, -2 },
        { &pdfArray::dbar, -1 },
        { &pdfArray::gluon, 21 },
        { &pdfArray::d, 1 },
        { &pdfArray::u, 2 },
        { &pdfArray::s, 3 },
        { &pdfArray::c, 4 },
        { &pdfArray::b, 5 },
        { &pdfArray::t, 6 },
        { &pdfArray::x, 9999 }, // TODO: what is `x`?
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

            double const result = lumi(f1, f2, c);

            pdg_ids.push_back(a.second);
            pdg_ids.push_back(b.second);
            factors.push_back(result);
        }
    }
}

/*
 * Initializes the pineappl grid with all the necessary parameters
 */
void CheffPanopoulos::create_grid(int max_orders, double q2) {
    auto* lumi = pineappl_lumi_new();

    // At LO we only have q Qb where q=Q (Z) or q != Q (W)
    // (diagonal ckm)

    // TODO: start with Z only
    int pdg_ids[] = {
        1, -1, -1, 1,
        2, -2, -2, 2,
        3, -3, -3, 3,
        4, -4, -4, 4,
    };
    double ckm_factors[8] = { 1.0 };
    pineappl_lumi_add(lumi, 8, pdg_ids, ckm_factors);

    // TODO only LO for now
    uint32_t orders[] = { 0, 2, 0, 0 };
    int how_many_orders = 1;
    
    // TODO
    double bins[] = {0.0, 1.0};
    int how_many_bins = 1;

    auto* keyval = pineappl_keyval_new();
    pineappl_keyval_set_double(keyval, "q2_min", 10); 
    grid = pineappl_grid_new(lumi, how_many_orders, orders, how_many_bins, bins, keyval);
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
void CheffPanopoulos::fill_grid(int order, int lumi_channel, double x1, double x2, double weight) {
    pineappl_grid_fill(grid, x1, x2, constant_q2, order, 0.5, lumi_channel, weight);

}

void CheffPanopoulos::save() {
    char const* filename = "test.pineappl.lz4";
    pineappl_grid_write(grid, filename);
    pineappl_grid_delete(grid);
}
