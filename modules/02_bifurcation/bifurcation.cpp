#define BOOST_PYTHON_STATIC_LIB
#include <boost/python.hpp>
#include "math_model.h"

namespace py = boost::python;

void RKutt(std::vector<double> &a, std::vector<double>&b,double x0, double y0, double z0, size_t t, double h, double time0, double w_q, size_t &iter) {
  size_t size = static_cast<int>(t / h + 2);
  std::vector<double> x(size);
  std::vector<double> y(size);
  std::vector<double> z(size);
  double yt = 0;
  double xt = 0;
  double zt = 0;
  double t0 = time0;
  size_t i = 0;
  x[i] = x0;
  y[i] = y0;
  z[i] = z0;
  
  for (t0; t0 < t; t0 += h) {
    double k1_y = fy(x[i], y[i], t0, w_q);
    double k2_y = fy(x[i] + k1_y / 2.0, y[i] + k1_y / 2.0, t0, w_q);
    double k3_y = fy(x[i] + k2_y / 2.0, y[i] + k2_y / 2.0, t0, w_q);
    double k4_y = fy(x[i] + k3_y, y[i] + k3_y, t0, w_q);
    yt = y[i] + h * (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
    y[i + 1] = yt;

    double k1_x = fx(x[i], y[i], t0, w_q);
    double k2_x = fx(x[i] + k1_x / 2.0, y[i] + k1_x / 2.0, t0, w_q);
    double k3_x = fx(x[i] + k2_x / 2.0, y[i] + k2_x / 2.0, t0, w_q);
    double k4_x = fx(x[i] + k3_x, y[i] + k3_x, t0, w_q);
    xt = x[i] + h * (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
    x[i + 1] = xt;

    double k1_z = fz(x[i], z[i], t0, w_q);
    double k2_z = fz(x[i] + k1_z / 2.0, z[i] + k1_z / 2.0, t0, w_q);
    double k3_z = fz(x[i] + k2_z / 2.0, z[i] + k2_z / 2.0, t0, w_q);
    double k4_z = fz(x[i] + k3_z, z[i] + k3_z, t0, w_q);
    zt = z[i] + h * (k1_z + 2 * k2_z + 2 * k3_z + k4_z) / 6;
    z[i + 1] = zt;

    if ((t0 > t_porog) && ((z[i - 2] - z[i - 1]) * (z[i - 1] - z[i]) < 0) && ((z[i - 2] - z[i - 1]) < 0)) {//
      a[iter] = z[i - 1];
      b[iter] = w_q;
      ++iter;
    }
    ++i;
  }
}

void bifurcation(py::list localM, py::list w_s) {
  size_t iter = 0;

  std::vector<double> lmax(9000000);
  std::vector<double> w_SS(9000000);
  int ili = static_cast<int>(i_py * 100000);
  
  int Iend = static_cast<int>(i_end * 100000);
  double start = omp_get_wtime();
#pragma omp parallel for shared(iter)
  for (int i = ili; i < Iend; i += 1) {
    double w_q = static_cast<double>(i) / 100000;
    RKutt(lmax, w_SS, x0, y0_, z0_, t, h, t0, w_q, iter);
  }
  double end = omp_get_wtime();
  printf("Parallel work time %f seconds\n", end - start);

  lmax.erase(lmax.begin() + iter, lmax.end());
  w_SS.erase(w_SS.begin() + iter, w_SS.end());
  
  std::vector<double>(lmax).swap(lmax);
  std::vector<double>(w_SS).swap(w_SS);

  lmax.erase(std::remove(lmax.begin(), lmax.end(), 0.0), lmax.end());
  w_SS.erase(std::remove(w_SS.begin(), w_SS.end(), 0.0), w_SS.end());

  for (size_t k = 0; k < lmax.size(); ++k)
    localM.append(lmax[k]);

  for (size_t j = 0; j < w_SS.size(); ++j)
    w_s.append(w_SS[j]);
}

BOOST_PYTHON_MODULE(bifurcation) {
  py::def("bifurcation", bifurcation);
}
