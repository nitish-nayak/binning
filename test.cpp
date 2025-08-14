#include "tree_binning.h"
#include <iostream>

int main(){
    std::vector<std::array<double,2>> sig = {
        {1,0.5},{1.2,0.4},{2,1.1},{2.2,1.3}
    };
    std::vector<std::array<double,2>> bkg = {
        {0.5,0.7},{1.1,0.9},{1.8,0.6},{2.5,1.8}
    };

    // a 2d binning with a simple FOM for splitting signal and background that maximizes statistical power
    TreeBinning b(dim_t<2>{}, 4, [](double s, double b){ return s*s/(s+b); });
    b.fit(sig, bkg);

    for (auto& leaf : b.leaves()) {
        std::cout << "x in [" << leaf.xmin[0] << ","<<leaf.xmax[0]<<") "
                  << "y in [" << leaf.xmin[1] << ","<<leaf.xmax[1]<<") "
                  << "(S="<<leaf.S<<", B="<<leaf.B<<")\n";
    }
}
