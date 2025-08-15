#include "../tree_binning.h"

#include "TH1D.h"
#include "TH2Poly.h"

template<size_t D>
class root_histogram_t {
public:
    root_histogram_t(histogram_t<D> hist):
        f_hist(std::make_shared<histogram_t<D>>(hist)) {}

    // fill a TH2Poly object
    TH2Poly* ToTH2Poly(const char* name="h") const {
      static_assert(D == 2, "I need a 2-dim object to convert to a TH2Poly!");
      // ROOT has a nice TH2Poly object just for our purpose
      TH2Poly* ret = new TH2Poly(name, "", f_hist->xmin[0], f_hist->xmax[0], f_hist->xmin[1], f_hist->xmax[1]);
      for(int i = 0; i < f_hist->nbins; i++){
        auto& bin = (f_hist->sorted_bins).at(i);
        int bin_i = ret->AddBin(bin.xmin[0], bin.xmin[1], bin.xmax[0], bin.xmax[1]);
        ret->SetBinContent(bin_i, bin.content);
      }
      return ret;
    }

    // collapse everything to a TH1D
    TH1D* ToTH1DCollapsed(const char* name="h") const {
      static_assert(D > 0, "I need atleast a 1-dim object to collapse to a TH1!");
      // collapse incoming histogram to a 1D histogram just based on bin indices
      TH1D* ret = new TH1D(name, "", f_hist->nbins, 0, f_hist->nbins);
      for(int i = 0; i < f_hist->nbins; i++){
        ret->SetBinContent(i+1, (f_hist->sorted_bins).at(i).content);
      }
      return ret;
    }

    // if I have a 1D histogram, then convert that to a TH1D
    TH1D* ToTH1D(const char* name="h") const {
      static_assert(D == 1, "I need a 1-dim object to convert to a TH1!");
      // first get an array to the bin edges
      std::vector<double> edges;
      for(int i = 0; i < f_hist->nbins; i++){
        edges.push_back(f_hist->sorted_bins.at(i).xmin[0]);
        // fill in the last boundary as well
        if(i == f_hist->nbins - 1)
          edges.push_back(f_hist->sorted_bins.at(i).xmax[0]);
      }
      TH1D* ret = new TH1D(name, "", f_hist->nbins, &edges[0]);
      for(int i = 0; i < f_hist->nbins; i++){
        ret->SetBinContent(i+1, (f_hist->sorted_bins).at(i).content);
      }
      return ret;
    }

private:
    std::shared_ptr<histogram_t<D>> f_hist;
};
