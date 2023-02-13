#ifndef MATH_MODEL_H_
#define MATH_MODEL_H_
#include <cmath>
#include <vector>
#include <omp.h>
#include <iostream>
#include <algorithm>

const double a_p = 0.001, a_z = 0.0001, alfa_0 = 0.23, b_p = 0.01, b_z = 0.01, g_p = 0.001, q_p = 6,
q_z = 5.5, Z0 = 0, Z1 = 1, P0 = 0, P1 = 1, k_z = 0.15, k_p = 0.05, Q0 = 1.8, Qs = 3.2, t_porog = 125000,
a_r = 0.01, b_r = 0.01, R0 = 2, R1 = 1, q_r = 5.6, k_r = 0.1;

const size_t t = 200000;
const double t0 = 0, h = 0.1, x0 = 0.2, y0_ = 0.2, z0_ = 0.2;
//double i_py = 0.0, i_end = 0.012;
extern double i_py, i_end;

double fx(double Z, double P, double t, double w_q);

double fy(double Z, double P, double t, double w_q);

double fz(double Z, double R, double t, double w_q);

double dot(const std::vector<double>& vec1, const std::vector<double>& vec2);

std::vector<double> addit_vec(std::vector<double>& vec1, std::vector<double>& vec2, int fl);

std::vector<double> minus_vec(std::vector<double>& vec1, std::vector<double>& vec2, double fl);

double norm(std::vector<double>& vec);

std::vector<double> Jacobian(const std::vector<double>& df, double Z, double P, double R, double t, double w_q);

#endif  // MATH_MODEL_H_