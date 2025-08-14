#ifndef TREE_BINNING_H
#define TREE_BINNING_H

#include <array>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <memory>
#include <cmath>
#include <limits>
#include <queue>
#include <numeric>
#include <functional>
#include <type_traits>

template <size_t D>
struct event_t {
  std::array<double, D> x;
  double weight;
  bool is_signal;
  constexpr size_t dim() { return D; }
};

template <size_t D>
struct node_t {
  // range of events that the node is looking at for splitting
  size_t idx0, idx1;
  std::array<double, D> xmin, xmax;
  double S = 0, B = 0;

  node_t(size_t i0, size_t i1,
       const std::array<double, D>& lo,
       const std::array<double, D>& hi):
    idx0(i0), idx1(i1), xmin(lo), xmax(hi) {}

  std::unique_ptr<node_t<D>> left, right;

  // scan across events based on given range
  void extract_counts(const std::vector<event_t<D>>& events){
    S = 0, B = 0;
    for (size_t k = idx0; k < idx1; ++k) {
      auto& e = events[k];
      (e.is_signal ? S : B) += e.weight;
    }
  }
  double x_split;
  size_t dim_split;
};

template <size_t D>
struct node_info_t {
  node_t<D>* node;
  double gain;
  // order splits based on their gain
  bool operator<(node_info_t const& o) const {
    return gain < o.gain;
  }
};

// this will read nicer
template<size_t D>
using dim_t = std::integral_constant<size_t, D>;

/// Greedy, priority-queue‐driven KD‐tree binning
/// the nice priority-queue and partition idea was suggested by chatgpt
template <size_t D, typename... Args>
class TreeBinning {
private:
  int f_maxleaves;
  std::function<double(double, double, Args...)> funcFOM;
  std::tuple<Args...> f_extraArgs;

public:
  // if I want to deduce Args from the arguments defined for func, I have to deduce everything else including the dimension which is of non-type here
  // which means the arguments to the ctor needed for the deduction needs to be known at compile time (dim, func)
  // since C++17 I can do this cleanly like below using integral_constants
  // otherwise I have to do weird metaprogramming things like pass in a C-array of size D like (&)[D] and guide the deduction like we do after this
  // ofcourse I also want to enforce funcFOM to take in the first two arguments as double (for signal and bkg counts)
  // and keep Args for the rest
  template<typename F>
  TreeBinning(dim_t<D> /* dim */, int max_leaves, F&& func):
    f_maxleaves(max_leaves), funcFOM(std::forward<F>(func)) {}


  struct histogram_t {
    std::array<double, D> xmin, xmax;
    double S, B;
  };

  const std::vector<histogram_t>& leaves() const { return f_bins; }

  /// Inputs: lists of D-vectors for signal & background,
  /// optional same-length weight arrays.
  void fit(Args&& ... args,
           const std::vector<std::array<double,D>>& sig_x,
           const std::vector<std::array<double,D>>& bkg_x,
           const std::vector<double>& sig_w = {},
           const std::vector<double>& bkg_w = {})
  {
    f_extraArgs = std::make_tuple(std::forward<Args>(args)...);
    // build events
    double totalS = 0, totalB = 0;
    f_events.clear();
    for(size_t i = 0; i < sig_x.size(); ++i) {
      double sig_wgt = sig_w.empty() ? 1.0 : sig_w[i];
      totalS += sig_wgt;
      f_events.push_back({sig_x[i], sig_wgt, true});
    }
    for(size_t i = 0; i < bkg_x.size(); ++i) {
      double bkg_wgt = bkg_w.empty() ? 1.0 : bkg_w[i];
      totalB += bkg_wgt;
      f_events.push_back({bkg_x[i], bkg_wgt, false});
    }

    // root covers [0,f_events.size()), infinite bounds
    std::array<double,D> lo, hi;
    lo.fill(-std::numeric_limits<double>::infinity());
    hi.fill( std::numeric_limits<double>::infinity());
    root = std::make_unique<node_t<D>>(0, f_events.size(), lo, hi);
    root->S = totalS;
    root->B = totalB;

    // greedy splitting
    grow_tree();

    // collect leaves
    collect_leaves(root.get(), f_bins);
  }

private:
  std::unique_ptr<node_t<D>> root;
  std::vector<event_t<D>> f_events;
  std::vector<histogram_t> f_bins;
  // methods to grow the tree
  void grow_tree() {
    std::priority_queue<node_info_t<D>> pq;
    split_node(root.get(), pq);

    int leaves = 1;
    while (leaves < f_maxleaves && !pq.empty()) {
      // pick the best leaf to split
      auto best = pq.top(); pq.pop();
      node_t<D>* n = best.node;

      // Partition f_events[idx0,i1) in place by n->x_split on n->dim_split
      auto begin = f_events.begin() + n->idx0;
      auto end   = f_events.begin() + n->idx1;
      auto split_it  = std::partition(begin, end,
          [&](event_t<D> const& e) {
            return e.x[n->dim_split] < n->x_split;
      });
      size_t i_split = std::distance(f_events.begin(), split_it);

      // Left child
      auto loL = n->xmin, hiL = n->xmax;
      hiL[n->dim_split] = n->x_split;
      n->left = std::make_unique<node_t<D>>(n->idx0, i_split, loL, hiL);
      (n->left)->extract_counts(f_events);

      // Right child
      auto loR = n->xmin, hiR = n->xmax;
      loR[n->dim_split] = n->x_split;
      n->right = std::make_unique<node_t<D>>(i_split, n->idx1, loR, hiR);
      (n->right)->extract_counts(f_events);

      ++leaves;
      split_node(n->left.get(),  pq);
      split_node(n->right.get(), pq);
    }
  }

  // Find the best split for node n and push to pq
  void split_node(node_t<D>* n, std::priority_queue<node_info_t<D>>& pq)
  {
    size_t N = n->idx1 - n->idx0;
    if (N < 2) return;
    // start at the current node and get its fom
    double parent_fom = std::apply(
        [&](auto&... args){ return funcFOM(n->S, n->B, args...); },
        f_extraArgs);

    // local offsets 0..N-1 for sorting
    std::vector<size_t> local(N);
    std::iota(local.begin(), local.end(), 0);
    double split_gain = 0;

    for (size_t dim = 0; dim < D; ++dim) {
      // sort offsets by the variable
      std::sort(local.begin(), local.end(),
          [&](size_t a, size_t b) {
            return f_events[n->idx0 + a].x[dim] < f_events[n->idx0 + b].x[dim];
      });

      double s_left = 0, b_left = 0;

      for (size_t j = 1; j < N; ++j) {
        // scan across our f_events vector that's sorted in place already
        auto& min_event = f_events[n->idx0 + local[j-1]];
        if(min_event.is_signal)
          s_left += min_event.weight;
        else
          b_left += min_event.weight;
        auto& max_event = f_events[n->idx0 + local[j]];
        if (max_event.x[dim] == min_event.x[dim]) continue;

        // evaluate information gain at current position
        double xval    = 0.5*(min_event.x[dim] + max_event.x[dim]);
        double s_right = n->S - s_left;
        double b_right = n->B - b_left;

        double fom_l = std::apply(
            [&](auto&&... args){ return funcFOM(s_left, b_left, args...); },
            f_extraArgs);
        double fom_r = std::apply(
            [&](auto&&... args){ return funcFOM(s_right, b_right, args...); },
            f_extraArgs);
        double gain = (fom_l + fom_r) - parent_fom;

        if (gain > split_gain) {
          split_gain = gain;
          n->x_split = xval;
          n->dim_split = dim;
        }
      }
    }
    // push into the queue
    if (split_gain > 0) {
      pq.push({n, split_gain});
    }
  }

  // Recursively collect leaf hyperrectangles
  void collect_leaves(node_t<D>* n, std::vector<histogram_t>& out) const {
    if (!n->left) {
      out.push_back({n->xmin, n->xmax, n->S, n->B});
    } else {
      collect_leaves(n->left.get(),  out);
      collect_leaves(n->right.get(), out);
    }
  }
};

// write out the deduction guide here
template<size_t D, typename R, typename... Args>
TreeBinning(dim_t<D>, int, std::function<R(double, double, Args...)>) -> TreeBinning<D, double, double, Args...>;

#endif
