#include "../tree_binning.h"
#include "root.h"
#include "TCanvas.h"

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

    root_histogram_t<2> signal_root_h(signal_h);
    root_histogram_t<2> bkg_root_h(bkg_h);

    TH2Poly* hsignal = signal_root_h.ToTH2Poly("signal2d");
    TH1D* hbkg = bkg_root_h.ToTH1DCollapsed("bkg1d");

    std::cout << "Signal, Bkg Counts : " << hsignal->Integral() << ", " << hbkg->Integral() << std::endl;

    TCanvas* c = new TCanvas();
    c->Divide(2, 1);

    c->cd(1);
    hsignal->Draw("colz");
    c->cd(2);
    hbkg->Draw("hist");
    c->cd();
    c->Print("roottest.pdf");

    return 0;
}
