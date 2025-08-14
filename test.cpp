#include "tree_binning.h"
#include <iostream>

int main(){
    std::vector<std::array<double,2>> sig = {
        {1,0.5},{1.2,0.4},{2,1.1},{2.2,1.3}
    };
    std::vector<std::array<double,2>> bkg = {
        {0.5,0.7},{1.1,0.9},{1.8,0.6},{2.5,1.8}
    };

    // a 2d binning with maximum 4 leaves and a simple FOM for splitting signal and background that maximizes statistical power
    TreeBinning b(dim_t<2>{}, 4, [](double s, double b){ return s*s/(s+b); });
    b.fit(sig, bkg);

    histogram_t<2> signal_h = b.signal_hist();
    histogram_t<2> bkg_h = b.bkg_hist();
    std::cout << "Bounds (x) : " << signal_h.xmin[0] << ", " << signal_h.xmax[0] << std::endl;
    std::cout << "Bounds (y) : " << signal_h.xmin[1] << ", " << signal_h.xmax[1] << std::endl;
    for (auto it = signal_h.sorted_bins.begin(); it != signal_h.sorted_bins.end(); ++it) {
        auto& s_leaf = *it;
        auto& b_leaf = bkg_h.sorted_bins.at(std::distance(signal_h.sorted_bins.begin(), it));
        std::cout << "x in [" << s_leaf.xmin[0] << "," <<s_leaf.xmax[0] << "), "
                  << "y in [" << s_leaf.xmin[1] << "," <<s_leaf.xmax[1] << ") "
                  << "(S = "  << s_leaf.content <<", B = " <<b_leaf.content << ")" << std::endl;
    }

    return 0;
}
