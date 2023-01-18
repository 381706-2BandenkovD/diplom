#define BOOST_PYTHON_STATIC_LIB
#include <boost/python.hpp>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iostream>
#include <algorithm>

namespace py = boost::python;

const double a_p = 0.001, a_z = 0.0001, alfa_0 = 0.23, b_p = 0.01, b_z = 0.01, g_p = 0.001, q_p = 6,
q_z = 5.5, Z0 = 0, Z1 = 1, P0 = 0, P1 = 1, k_z = 0.15, k_p = 0.05, Q0 = 1.8, Qs = 3.2, t_porog = 125000,
a_r = 0.01, b_r = 0.01, R0 = 2, R1 = 1, q_r = 5.6, k_r = 0.1;

const size_t t = 200000;
const double t0 = 0, h = 0.01, x0 = 0.2, y0_ = 0.2, z0_ = 0.2;
//double i_py = 0.0, i_end = 0.012;
double i_py = 0.0045, i_end = 0.00584;

double fx(double Z, double P, double t, double w_q) {
    return (-(a_z + g_p * P) * Z + b_z * (Z0 - (Z0 - Z1) / (1 + exp(-(Q0 + alfa_0 * Z + Qs * sin(w_q * t) - q_z) / k_z))));
}

double fy(double Z, double P, double t, double w_q) {
    return (-a_p * P + b_p * (P0 - (P0 - P1) / (1 + exp(-(Q0 + alfa_0 * Z + Qs * sin(w_q * t) - q_p) / k_p))));
}

double fz(double Z, double R, double t, double w_q) {
    return (-a_r * R + b_r * (R0 - (R0 - R1) / (1 + exp(-(Q0 + alfa_0 * Z + Qs * sin(w_q * t) - q_r) / k_r))));
}

//-------------------------------LYAPUNOV-LOCAL-MAX------------------------------------------------------
double dot(const std::vector<double>& vec1, const std::vector<double>& vec2) {
  int size = vec1.size();
  if (size != vec2.size()) {
    std::cout << "V1 = " << vec1.size();
    std::cout << "\nV2 = " << vec2.size();
    throw 1;
  }
  double sum = 0;

  for (int i = 0; i < size; ++i)
    sum += vec1[i] * vec2[i];

  return sum;
}

std::vector<double> addit_vec(std::vector<double>& vec1, std::vector<double>& vec2, int fl) {
  int size = vec1.size();
  if (size != vec2.size()) {
    std::cout << "V1 = " << vec1.size();
    std::cout << "\nV2 = " << vec2.size();
    throw 1;
  }

  std::vector<double> tmp = vec2;
  if (fl == 1) {
    for (size_t i = 0; i < size; ++i)
      tmp[i] *= 0.5;
  }

  auto res = std::vector<double>(size);
  std::transform(vec1.begin(), vec1.end(), tmp.begin(), res.begin(), std::plus<double>());
  return res;
}

std::vector<double> minus_vec(std::vector<double>& vec1, std::vector<double>& vec2, double fl) {
  int size = vec1.size();
  if (size != vec2.size())
    throw 1;

  std::vector<double> tmp = vec2;
  for (size_t i = 0; i < size; ++i)
    tmp[i] *= fl;

  auto res = std::vector<double>(size);
  std::transform(vec1.begin(), vec1.end(), tmp.begin(), res.begin(), std::minus<double>());
  return res;
}

double norm(std::vector<double>& vec) {
  return sqrt(pow(vec[0], 2) + pow(vec[1], 2)+ pow(vec[2], 2));
}

std::vector<double> Jacobian(const std::vector<double>& df, double Z, double P, double R, double t, double w_q) {
  double fx_coff = (exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_z) / k_z) + 1);
  double fy_coff = (exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_p) / k_p) + 1);
  double fz_coff = (exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_r) / k_r) + 1);

  double fx_Z = -P * g_p - a_z - alfa_0 * b_z * (Z0 - Z1) * exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_z) / k_z) / (k_z * fx_coff * fx_coff);
  double fx_P = -Z * g_p;
  double fx_R = 0;
  double fy_Z = -alfa_0 * b_p * (P0 - P1) * exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_p) / k_p) / (k_p * fy_coff * fy_coff);
  double fy_P = -a_p;
  double fy_R = 0;
  double fz_Z = -alfa_0 * b_r * (R0 - R1) * exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_r) / k_r) / (k_r * fz_coff * fz_coff);
  double fz_P =  0;
  double fz_R =  -a_r;


  std::vector<double> d_fx{ fx_Z , fx_P, fx_R };
  std::vector<double> d_fy{ fy_Z, fy_P, fy_R };
  std::vector<double> d_fz{ fz_Z, fz_P, fz_R };

  double d1 = dot(d_fx, df) * h;
  double d2 = dot(d_fy, df) * h;
  double d3 = dot(d_fz, df) * h;
  std::vector<double> res = { d1, d2, d3 };
  return res;
}

void RK_lyapunov(std::vector<double> &lmax, std::vector<double>&w_ss, std::vector<double>& lya1, std::vector<double>& lya2,
                 std::vector<double>& lya3, double w_q, size_t &iter, size_t &j) {

  size_t size = static_cast<int>(t / h + 2);
  std::vector<double> x(size);
  std::vector<double> y(size);
  std::vector<double> z(size);
  
  double xt = 0;
  double yt = 0;
  double zt = 0;
  double l1 = 0;
  double l2 = 0;
  double l3 = 0;

  std::vector<double> dfx{ 1, 0, 0 };
  std::vector<double> dfy{ 0, 1, 0 };
  std::vector<double> dfz{ 0, 0, 1 };
  std::vector<double> ort1(3);
  std::vector<double> ort2(3);
  std::vector<double> ort3(3);

  double t0 = 0;
  size_t i = 0;
  x[i] = x0;
  y[i] = y0_;
  z[i] = z0_;
  
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

    std::vector<double> k11_x = Jacobian(dfx, x[i], y[i], z[i], t0, w_q);
    std::vector<double> tmp = addit_vec(dfx, k11_x, 1);
    std::vector<double> k22_x = Jacobian(tmp, x[i] + k1_x / 2.0 * h, y[i] + k1_x / 2.0 * h, z[i] + k1_x / 2.0 * h, t0, w_q);
    tmp = addit_vec(dfx, k22_x, 1);
    std::vector<double> k33_x = Jacobian(tmp, x[i] + k2_x / 2.0 * h, y[i] + k2_x / 2.0 * h, z[i] + k2_x / 2.0 * h, t0, w_q);
    tmp = addit_vec(dfx, k33_x, 0);
    std::vector<double> k44_x = Jacobian(tmp, x[i] + k3_x * h, y[i] + k3_x * h, z[i] + k3_x * h, t0, w_q);

    std::vector<double> k11_y = Jacobian(dfy, x[i], y[i], z[i], t0, w_q);
    std::vector<double> tmpy = addit_vec(dfy, k11_y, 1);
    std::vector<double> k22_y = Jacobian(tmpy, x[i] + k1_y / 2.0 * h, y[i] + k1_y / 2.0 * h, z[i] + k1_y / 2.0 * h, t0, w_q);
    tmpy = addit_vec(dfy, k22_y, 1);
    std::vector<double> k33_y = Jacobian(tmpy, x[i] + k2_y / 2.0 * h, y[i] + k2_y / 2.0 * h, z[i] + k2_y / 2.0 * h, t0, w_q);
    tmpy = addit_vec(dfy, k33_y, 0);
    std::vector<double> k44_y = Jacobian(tmpy, x[i] + k3_y * h, y[i] + k3_y * h, z[i] + k3_y * h, t0, w_q);

    std::vector<double> k11_z = Jacobian(dfz, x[i], y[i], z[i], t0, w_q);
    std::vector<double> tmpz = addit_vec(dfz, k11_z, 1);
    std::vector<double> k22_z = Jacobian(tmpz, x[i] + k1_z / 2.0 * h, y[i] + k1_z / 2.0 * h, z[i] + k1_z / 2.0 * h, t0, w_q);
    tmpz = addit_vec(dfz, k22_z, 1);
    std::vector<double> k33_z = Jacobian(tmpz, x[i] + k2_z / 2.0 * h, y[i] + k2_z / 2.0 * h, z[i] + k2_z / 2.0 * h, t0, w_q);
    tmpz = addit_vec(dfz, k33_z, 0);
    std::vector<double> k44_z = Jacobian(tmpz, x[i] + k3_z * h, y[i] + k3_z * h, z[i] + k3_z * h, t0, w_q);

    dfx = { dfx[0] + (k11_x[0] + 2.0 * k22_x[0] + 2.0 * k33_x[0] + k44_x[0]) / 6.0, dfx[1] + (k11_x[1] + 2.0 * k22_x[1] + 2.0 * k33_x[1] + k44_x[1]) / 6.0, 
            dfx[2] + (k11_x[2] + 2.0 * k22_x[2] + 2.0 * k33_x[2] + k44_x[2]) / 6.0 };
    dfy = { dfy[0] + (k11_y[0] + 2.0 * k22_y[0] + 2.0 * k33_y[0] + k44_y[0]) / 6.0, dfy[1] + (k11_y[1] + 2.0 * k22_y[1] + 2.0 * k33_y[1] + k44_y[1]) / 6.0,
            dfy[2] + (k11_y[2] + 2.0 * k22_y[2] + 2.0 * k33_y[2] + k44_y[2]) / 6.0 };
    dfz = { dfz[0] + (k11_z[0] + 2.0 * k22_z[0] + 2.0 * k33_z[0] + k44_z[0]) / 6.0, dfz[1] + (k11_z[1] + 2.0 * k22_z[1] + 2.0 * k33_z[1] + k44_z[1]) / 6.0,
            dfz[2] + (k11_z[2] + 2.0 * k22_z[2] + 2.0 * k33_z[2] + k44_z[2]) / 6.0 };

    ort1 = dfx;
    double n1 = norm(ort1);
    l1 += log(n1);
    dfx = { ort1[0] / n1, ort1[1] / n1, ort1[2] / n1 };
    
    double dot_yx = dot(dfy, dfx);
    ort2 = minus_vec(dfy, dfx, dot_yx);
    double n2 = norm(ort2);
    l2 += log(n2);
    dfy = { ort2[0] / n2, ort2[1] / n2, ort2[1] / n2 };

    double dot_zy = dot(dfz, dfy);
    double dot_zx = dot(dfz, dfx);
    ort3 = minus_vec(dfz, dfy, dot_zy);
    ort3 = minus_vec(ort3, dfx, dot_zx);
    double n3 = norm(ort3);
    l3 += log(n3);
    dfz = { ort3[0] / n3, ort3[1] / n3, ort3[2] / n3 };

    if ((t0 > t_porog) && ((z[i - 2] - z[i - 1]) * (z[i - 1] - z[i]) < 0) && ((z[i - 2] - z[i - 1]) < 0)) {//
      lmax[iter] = z[i - 1];
      w_ss[iter] = w_q;
      ++iter;
    }
    ++i;
  }

  double lk1, lk2, lk3;
  double TT = static_cast<double>(t);
  lk1 = (l1 / TT) / log(2);
  lk2 = (l2 / TT) / log(2);
  lk3 = (l3 / TT) / log(2);
  lya1[j] = lk1;
  lya2[j] = lk2;
  lya3[j] = lk3;
}

void lyapunov_solution(py::list localM, py::list w_s, py::list ly1, py::list ly2, py::list ly3, py::list w) {
  size_t iter = 0;
  size_t j = 0;
  std::vector<double> lmax(9000000);
  std::vector<double> w_SS(9000000);
  int ili = static_cast<int>(i_py * 1000000);
  int Iend = static_cast<int>(i_end * 1000000);
  std::vector<double> lyapun1(Iend - ili);
  std::vector<double> lyapun2(Iend - ili);
   std::vector<double> lyapun3(Iend - ili);
  std::vector<double> wwwq(Iend - ili);

  double start = omp_get_wtime();
#pragma omp parallel for shared(iter, j)
  for (int i = ili; i < Iend; i += 1) {
    double w_q = static_cast<double>(i) / 1000000;
    RK_lyapunov(lmax, w_SS, lyapun1, lyapun2, lyapun3, w_q, iter, j);
    wwwq[j] = w_q;
    ++j;
  }
  double end = omp_get_wtime();
  printf("Parallel lyapunov work time %f seconds\n", end - start);

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

  for (size_t l_i = 0; l_i < lyapun1.size(); ++l_i)
    ly1.append(lyapun1[l_i]);

  for (size_t l_j = 0; l_j < lyapun2.size(); ++l_j)
    ly2.append(lyapun2[l_j]);

  for (size_t l_k = 0; l_k < lyapun3.size(); ++l_k)
    ly3.append(lyapun3[l_k]);

  for (size_t j = 0; j < wwwq.size(); ++j)
    w.append(wwwq[j]);
}

//-------------------------------LYAPUNOV-REALISATION------------------------------------------------------
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
//----------------------BIFURCATION----------------------------------------------------------------
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

//-----------------------REALISATION---------------------------------------------------------------

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

void retuRelQQ(py::list funX, py::list funY, py::list funZ, py::list funT, double w_q, double qq) {
  size_t size = static_cast<int>(t / h + 1);
  std::vector<double> x(size + 1);
  std::vector<double> y(size + 1);
  std::vector<double> z(size + 1);
  std::vector<double> time(size + 1);

  //q_r = qq;

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

BOOST_PYTHON_MODULE(lyapunov) {
  py::def("bifurcation", bifurcation);
  py::def("retuRelQQ", retuRelQQ);
  py::def("retuRelease", retuRelease);
  py::def("lyapunov_solution", lyapunov_solution);
  py::def("RK_ly", RK_ly);
}
