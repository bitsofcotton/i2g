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
#include "../goki_check_cc/p0.hh"
#include "../goki_check_cc/p1.hh"
#include "../goki_check_cc/decompose.hh"
#include "../goki_check_cc/catg.hh"
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
    auto in1(rgb2d<num_t>(in0));
    for(int i = 0; i < in1.rows(); i ++)
      for(int j = 0; j < in1.cols(); j ++)
        in1(i, j) = in1(i, j) < half ? num_t(int(1)) : num_t(int(0));
    auto in(in1);
    const int shy[8] = {- 1, - 1, - 1, 0, 1, 1,   1,   0};
    const int shx[8] = {- 1,   0,   1, 1, 1, 0, - 1, - 1};
    for(int k = 0; k < 8; k ++)
      for(int i = 0; i < in1.rows(); i ++)
        for(int j = 0; j < in1.cols(); j ++)
          in(max(0, min(in.rows() - 1, shy[k] + i)),
             max(0, min(in.cols() - 1, shx[k] + j))) += in1(i, j);
    std::vector<std::pair<int, int> > p;
    p.reserve(in.rows() * in.cols());
    for(int i = 0; i < in.rows(); i ++)
      for(int j = 0; j < in.cols(); j ++)
        if(in(i, j) == num_t(int(1))) {
          p.emplace_back(make_pair(i, j));
          break;
        }
    auto mask(in);
    mask.O();
    mask(p[p.size() - 1].first, p[p.size() - 1].second) = num_t(int(1));
    int idx(7);
    while(true) {
      auto& lp(p[p.size() - 1]);
      int di;
      for(di = 0; di < 8; di ++) {
        auto q(make_pair(lp.first + shy[(idx + di) % 8], lp.second + shx[(idx + di) % 8]));
        if(0 <= q.first && q.first < in.rows() &&
           0 <= q.second && q.second < in.cols() &&
           in(q.first, q.second) == num_t(int(1)) &&
           mask(q.first, q.second) < half) {
          p.emplace_back(q);
          idx = (idx + di + 6) % 8;
          break;
        }
      }
      if(di == 8) break;
      mask(p[p.size() - 1].first, p[p.size() - 1].second) = num_t(int(1));
    }
    SimpleVector<num_t> py(p.size());
    SimpleVector<num_t> px(p.size());
    for(int i = 0; i < p.size(); i ++) {
      py[i] = p[i].first;
      px[i] = p[i].second;
    }
    SimpleMatrix<num_t> out(2, sz);
    out.row(0) = (dft<num_t>(- sz) * dft<num_t>(max(sz, py.size())).subMatrix(0, 0, sz, py.size()) * py.template cast<complex<num_t> >()).template real<num_t>();
    out.row(1) = (dft<num_t>(- sz) * dft<num_t>(max(sz, px.size())).subMatrix(0, 0, sz, px.size()) * px.template cast<complex<num_t> >()).template real<num_t>();
    vector<SimpleMatrix<num_t> > out0;
    out0.emplace_back(std::move(out));
    savep2or3<num_t>(argv[4], normalize<num_t>(out0), true);
  } else if(m == '-') {
    std::vector<SimpleMatrix<num_t> > in0;
    if(! loadp2or3<num_t>(in0, argv[3])) return - 1;
    auto in(rgb2d<num_t>(in0));
    SimpleMatrix<num_t> out(sz, sz);
    out.O(num_t(int(1)));
    auto my(in(0, 0));
    auto mx(in(1, 0));
    auto My(in(0, 0));
    auto Mx(in(1, 0));
    for(int i = 0; i < in.cols(); i ++) {
      my = min(my, in(0, i));
      mx = min(mx, in(1, i));
      My = max(My, in(0, i));
      Mx = max(Mx, in(1, i));
    }
    if(my == My) My += num_t(int(1));
    if(mx == Mx) Mx += num_t(int(1));
    for(int i = 0; i < in.cols(); i ++) {
      SimpleVector<num_t> st(2);
      SimpleVector<num_t> ed(2);
      st[0] = (in(0, i) - my) / (My - my) * num_t(sz);
      st[1] = (in(1, i) - mx) / (Mx - mx) * num_t(sz);
      ed[0] = (in(0, (i + 1) % in.cols()) - my) / (My - my) * num_t(sz);
      ed[1] = (in(1, (i + 1) % in.cols()) - mx) / (Mx - mx) * num_t(sz);
      drawMatchLine<num_t>(out, st, ed, num_t(int(0)));
    }
    std::vector<SimpleMatrix<num_t> > out0;
    out0.emplace_back(std::move(out));
    savep2or3<num_t>(argv[4], out0, true);
  }
  return 0;
}

