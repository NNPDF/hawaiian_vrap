#pragma once

#include <pineappl_capi.h>


namespace pinerap {

    class CheffPanopoulos {
        /*
        * Flavour up vrap with a bit of pineappl
        */
        public:
            // the constructor creates the grid with all different orders and channels
            void create_grid(int max_orders, double q2);
            void fill_grid(int order, int lumi_channel, double x1, double x2, double weight);
            void save();

        private:
            double constant_q2;
            pineappl_grid* grid;
    };
}
