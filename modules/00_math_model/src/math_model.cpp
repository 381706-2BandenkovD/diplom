#include "math_model.h"

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
