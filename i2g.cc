#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>
#include <random>
#include <assert.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

//#define int int64_t
#define int int32_t
#include "../goki_check_cc/lieonn.hh"
typedef myfloat num_t;
#include "../goki_check_cc/goki.hh"

using std::cout;
using std::cerr;
using std::endl;
using std::atoi;
using std::string;
using std::to_string;
using std::vector;
using std::sort;
using std::binary_search;
using std::make_pair;

#include <stdlib.h>

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  assert(1 < argc);
  const auto m(argv[1][0]);
  const auto sz(std::atoi(argv[2]));
  if(m == '+') {
    std::vector<SimpleMatrix<num_t> > in0;
    if(! loadp2or3<num_t>(in0, argv[3])) return - 1;
    auto half(num_t(int(1)) / num_t(int(2)));
    auto in(rgb2d<num_t>(in0));
    std::vector<std::pair<int, int> > p;
    p.reserve(in.rows() * in.cols());
    for(int i = 0; i < in.rows(); i ++)
      for(int j = 0; j < in.cols(); j ++)
        if(in(i, j) < half) {
          p.emplace_back(make_pair(i, j));
          break;
        }
    while(true) {
      auto& lp(p[p.size() - 1]);
      if(0 <= lp.first - 1 && 0 <= lp.second - 1 &&
         in(lp.first - 1, lp.second - 1) < half)
        p.emplace_back(make_pair(lp.first - 1, lp.second - 1));
      else if(0 <= lp.first - 1 &&
              in(lp.first - 1, lp.second) < half)
        p.emplace_back(make_pair(lp.first - 1, lp.second));
      else if(0 <= lp.first - 1 && lp.second + 1 < in.cols() &&
              in(lp.first - 1, lp.second + 1) < half)
        p.emplace_back(make_pair(lp.first - 1, lp.second + 1));
      else if(                            lp.second + 1 < in.cols() &&
              in(lp.first, lp.second + 1) < half)
        p.emplace_back(make_pair(lp.first, lp.second + 1));
      else if(lp.first + 1 < in.rows() && lp.second + 1 < in.cols() &&
              in(lp.first + 1, lp.second + 1) < half)
        p.emplace_back(make_pair(lp.first + 1, lp.second + 1));
      else if(lp.first + 1 < in.rows()
              in(lp.first + 1, lp.second) < half)
        p.emplace_back(make_pair(lp.first + 1, lp.second));
      else if(lp.first + 1 < in.rows() && 0 < lp.second - 1 &&
              in(lp.first + 1, lp.second - 1) < half)
        p.emplace_back(make_pair(lp.first + 1, lp.second - 1));
      else if(                            0 < lp.second - 1 &&
              in(lp.first, lp.second - 1) < half)
      else break;
      auto q(p);
      q.erase(q.end() - 1);
      std::sort(q.begin(), q.end());
      if(binary_search(q.begin(), q.end(), p[p.size() - 1])) break;
      SimpleVector<num_t> py(p.size());
      SimpleVector<num_t> px(p.size());
      for(int i = 0; i < p.size(); i ++) {
        py[i] = p[i].first;
        px[i] = p[i].second;
      }
      SimpleMatrix<num_t> out(2, sz);
      out.row(0) = (dft<num_t>(- sz) * dft<num_t>(max(sz, py.size())).subMatrix(0, 0, sz, py.size()) * py.template cast<complex<num_t> >()).template real<num_t>();
      out.row(1) = (dft<num_t>(- sz) * dft<num_t>(max(sz, px.size())).subMatrix(0, 0, sz, px.size()) * px.template cast<complex<num_t> >()).template real<num_t>();
      savep2or3<num_t>(argv[4], out, true);
    }
  } else if(m == '-') {
    
  }
  return 0;
}

