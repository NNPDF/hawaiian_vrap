#include "pineappl_interface.h"

#include <cstddef>
#include <iostream>

using namespace pinerap;
/*
 * Initializes the pineappl grid with all the necessary parameters
 */
void CheffPanopoulos::create_grid(int max_orders, double q2) {
    auto* lumi = pineappl_lumi_new();

    // At LO we only have q Qb where q=Q (Z) or q != Q (W)
    // (diagonal ckm)

    // TODO: start with Z only
    int pdg_ids[] = {
        1, -1,
        2, -2,
        3, -3,
        4, -4,
        5, -5
    };
    double ckm_factors[5] = { 1.0 };
    pineappl_lumi_add(lumi, 5, pdg_ids, ckm_factors);

    // TODO only LO for now
    uint32_t orders[] = { 0, 2, 0, 0 };
    int how_many_orders = 1;
    
    // TODO
    double bins[] = {0.0, 1.0};
    int how_many_bins = 1;

    grid = pineappl_grid_new(lumi, how_many_orders, orders, how_many_bins, bins, nullptr);
    constant_q2 = q2;
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
