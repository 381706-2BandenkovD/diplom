#define BOOST_PYTHON_STATIC_LIB
#include <boost/python.hpp>
#include "math_model.h"
namespace py = boost::python;


void RK_ly(py::list x_list, py::list y_list, py::list z_list, py::list time_list, double w_q) {

    size_t size = static_cast<int>(t / h + 2);
    std::vector<double> x(size);
    std::vector<double> y(size);
    std::vector<double> z(size);
    std::vector<double> time_vec(size);

    double xt = 0;
    double yt = 0;
    double zt = 0;
    double l1 = 0;
    double l2 = 0;
    double l3 = 0;

    std::vector<double> dfx{1, 0, 0};
    std::vector<double> dfy{0, 1, 0};
    std::vector<double> dfz{0, 0, 1};
    std::vector<double> ort1(3);
    std::vector<double> ort2(3);
    std::vector<double> ort3(3);

    double t0 = 0;
    size_t i = 0;
    x[i] = x0;
    y[i] = y0_;
    z[i] = z0_;
    time_vec[i] = 0;

    for (t0; t0 < t; t0 += h) {
        double k1_x = fx(x[i], y[i], t0, w_q);
        double k2_x = fx(x[i] + k1_x / 2.0, y[i] + k1_x / 2.0, t0, w_q);
        double k3_x = fx(x[i] + k2_x / 2.0, y[i] + k2_x / 2.0, t0, w_q);
        double k4_x = fx(x[i] + k3_x, y[i] + k3_x, t0, w_q);
        xt = x[i] + h * (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
        x[i + 1] = xt;

        double k1_y = fy(x[i], y[i], t0, w_q);
        double k2_y = fy(x[i] + k1_y / 2.0, y[i] + k1_y / 2.0, t0, w_q);
        double k3_y = fy(x[i] + k2_y / 2.0, y[i] + k2_y / 2.0, t0, w_q);
        double k4_y = fy(x[i] + k3_y, y[i] + k3_y, t0, w_q);
        yt = y[i] + h * (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
        y[i + 1] = yt;

        double k1_z = fz(x[i], z[i], t0, w_q);
        double k2_z = fz(x[i] + k1_z / 2.0, z[i] + k1_z / 2.0, t0, w_q);
        double k3_z = fz(x[i] + k2_z / 2.0, z[i] + k2_z / 2.0, t0, w_q);
        double k4_z = fz(x[i] + k3_z, z[i] + k3_z, t0, w_q);
        zt = z[i] + h * (k1_z + 2 * k2_z + 2 * k3_z + k4_z) / 6;
        z[i + 1] = zt;
        time_vec[i + 1] = t0;

        std::vector<double> k11_x = Jacobian(dfx, x[i], y[i], z[i], t0, w_q);
        std::vector<double> tmp = addit_vec(dfx, k11_x, 1);
        std::vector<double> k22_x = Jacobian(tmp, x[i] + k1_x / 2.0 * h, y[i] + k1_x / 2.0 * h,
                                             z[i] + k1_x / 2.0 * h, t0, w_q);
        tmp = addit_vec(dfx, k22_x, 1);
        std::vector<double> k33_x = Jacobian(tmp, x[i] + k2_x / 2.0 * h, y[i] + k2_x / 2.0 * h,
                                             z[i] + k2_x / 2.0 * h, t0, w_q);
        tmp = addit_vec(dfx, k33_x, 0);
        std::vector<double> k44_x =
            Jacobian(tmp, x[i] + k3_x * h, y[i] + k3_x * h, z[i] + k3_x * h, t0, w_q);

        std::vector<double> k11_y = Jacobian(dfy, x[i], y[i], z[i], t0, w_q);
        std::vector<double> tmpy = addit_vec(dfy, k11_y, 1);
        std::vector<double> k22_y = Jacobian(tmpy, x[i] + k1_y / 2.0 * h, y[i] + k1_y / 2.0 * h,
                                             z[i] + k1_y / 2.0 * h, t0, w_q);
        tmpy = addit_vec(dfy, k22_y, 1);
        std::vector<double> k33_y = Jacobian(tmpy, x[i] + k2_y / 2.0 * h, y[i] + k2_y / 2.0 * h,
                                             z[i] + k2_y / 2.0 * h, t0, w_q);
        tmpy = addit_vec(dfy, k33_y, 0);
        std::vector<double> k44_y =
            Jacobian(tmpy, x[i] + k3_y * h, y[i] + k3_y * h, z[i] + k3_y * h, t0, w_q);

        std::vector<double> k11_z = Jacobian(dfz, x[i], y[i], z[i], t0, w_q);
        std::vector<double> tmpz = addit_vec(dfz, k11_z, 1);
        std::vector<double> k22_z = Jacobian(tmpz, x[i] + k1_z / 2.0 * h, y[i] + k1_z / 2.0 * h,
                                             z[i] + k1_z / 2.0 * h, t0, w_q);
        tmpz = addit_vec(dfz, k22_z, 1);
        std::vector<double> k33_z = Jacobian(tmpz, x[i] + k2_z / 2.0 * h, y[i] + k2_z / 2.0 * h,
                                             z[i] + k2_z / 2.0 * h, t0, w_q);
        tmpz = addit_vec(dfz, k33_z, 0);
        std::vector<double> k44_z =
            Jacobian(tmpz, x[i] + k3_z * h, y[i] + k3_z * h, z[i] + k3_z * h, t0, w_q);

        dfx = {dfx[0] + (k11_x[0] + 2.0 * k22_x[0] + 2.0 * k33_x[0] + k44_x[0]) / 6.0,
               dfx[1] + (k11_x[1] + 2.0 * k22_x[1] + 2.0 * k33_x[1] + k44_x[1]) / 6.0,
               dfx[2] + (k11_x[2] + 2.0 * k22_x[2] + 2.0 * k33_x[2] + k44_x[2]) / 6.0};
        dfy = {dfy[0] + (k11_y[0] + 2.0 * k22_y[0] + 2.0 * k33_y[0] + k44_y[0]) / 6.0,
               dfy[1] + (k11_y[1] + 2.0 * k22_y[1] + 2.0 * k33_y[1] + k44_y[1]) / 6.0,
               dfy[2] + (k11_y[2] + 2.0 * k22_y[2] + 2.0 * k33_y[2] + k44_y[2]) / 6.0};
        dfz = {dfz[0] + (k11_z[0] + 2.0 * k22_z[0] + 2.0 * k33_z[0] + k44_z[0]) / 6.0,
               dfz[1] + (k11_z[1] + 2.0 * k22_z[1] + 2.0 * k33_z[1] + k44_z[1]) / 6.0,
               dfz[2] + (k11_z[2] + 2.0 * k22_z[2] + 2.0 * k33_z[2] + k44_z[2]) / 6.0};

        ort1 = dfx;
        double n1 = norm(ort1);
        l1 += log(n1);
        dfx = {ort1[0] / n1, ort1[1] / n1, ort1[2] / n1};

        double dot_yx = dot(dfy, dfx);
        ort2 = minus_vec(dfy, dfx, dot_yx);
        double n2 = norm(ort2);
        l2 += log(n2);
        dfy = {ort2[0] / n2, ort2[1] / n2, ort2[1] / n2};

        double dot_zy = dot(dfz, dfy);
        ort3 = minus_vec(dfz, dfy, dot_zy);
        double dot_zx = dot(dfz, dfx);
        ort3 = minus_vec(ort3, dfx, dot_zx);
        double n3 = norm(ort3);
        l3 += log(n3);
        dfz = {ort3[0] / n3, ort3[1] / n3, ort3[2] / n3};
        ++i;
    }

    double lk1, lk2, lk3;
    double TT = static_cast<double>(t);
    lk1 = (l1 / TT) / log(2);
    lk2 = (l2 / TT) / log(2);
    lk3 = (l3 / TT) / log(2);
    std::cout<<" lya1 = "<<lk1<<" lya2 = "<<lk2<<" lya3 = "<<lk3;

    for (size_t k = 0; k < x.size(); ++k)
      x_list.append(x[k]);
    for (size_t k = 0; k < x.size(); ++k)
      y_list.append(y[k]);
    for (size_t k = 0; k < x.size(); ++k)
      z_list.append(z[k]);
    for (size_t k = 0; k < x.size(); ++k)
      time_list.append(time_vec[k]);

}


void RK_wiki(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double> &time,double x0, double y0, 
              double z0, size_t t, double h, double time0, double w_q) {
  double yt = 0;
  double xt = 0;
  double zt = 0;
  double t0 = time0;
  size_t i = 0;
  time[i] = t0;
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

    time[i + 1] = t0;
    ++i;
  }
}

void retuRelease(py::list funX, py::list funY, py::list funZ, py::list funT, double w_q) {
  size_t size = static_cast<int>(t / h + 1);
  std::vector<double> x(size);
  std::vector<double> y(size);
  std::vector<double> z(size);
  std::vector<double> time(size);

  RK_wiki(x, y, z, time, x0, y0_, z0_, t, h, t0, w_q);

  for (size_t i = 0; i < x.size(); ++i)
    funX.append(x[i]);

  for (size_t i = 0; i < y.size(); ++i)
    funY.append(y[i]);

  for (size_t i = 0; i < z.size(); ++i)
    funZ.append(z[i]);

  for (size_t i = 0; i < time.size(); ++i)
    funT.append(time[i]);
}

BOOST_PYTHON_MODULE(time_series) {
  py::def("retuRelease", retuRelease);
  py::def("RK_ly", RK_ly);
}
