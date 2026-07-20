#include "kmat/stats.hpp"

#include <cmath>

namespace kmat {

namespace {

// Regularized incomplete beta I_x(a,b) via continued fraction (Numerical Recipes style).
double betacf(double a, double b, double x) {
  constexpr int max_iter = 200;
  constexpr double eps = 3.0e-7;
  constexpr double fpmin = 1.0e-30;

  const double qab = a + b;
  const double qap = a + 1.0;
  const double qam = a - 1.0;
  double c = 1.0;
  double d = 1.0 - qab * x / qap;
  if (std::abs(d) < fpmin) {
    d = fpmin;
  }
  d = 1.0 / d;
  double h = d;

  for (int m = 1; m <= max_iter; ++m) {
    const int m2 = 2 * m;
    double aa = m * (b - static_cast<double>(m)) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (std::abs(d) < fpmin) {
      d = fpmin;
    }
    c = 1.0 + aa / c;
    if (std::abs(c) < fpmin) {
      c = fpmin;
    }
    d = 1.0 / d;
    h *= d * c;

    aa = -(a + static_cast<double>(m)) * (qab + static_cast<double>(m)) * x /
         ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    if (std::abs(d) < fpmin) {
      d = fpmin;
    }
    c = 1.0 + aa / c;
    if (std::abs(c) < fpmin) {
      c = fpmin;
    }
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::abs(del - 1.0) < eps) {
      break;
    }
  }
  return h;
}

double betai(double a, double b, double x) {
  if (x <= 0.0) {
    return 0.0;
  }
  if (x >= 1.0) {
    return 1.0;
  }

  const double ln_beta =
      std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
  const double front = std::exp(std::log(x) * a + std::log(1.0 - x) * b - ln_beta) / a;

  if (x < (a + 1.0) / (a + b + 2.0)) {
    return front * betacf(a, b, x);
  }
  return 1.0 - front * betacf(b, a, 1.0 - x);
}

double student_t_cdf(double t, int df) {
  if (df <= 0) {
    return 0.5;
  }
  const double x = static_cast<double>(df) / (static_cast<double>(df) + t * t);
  const double a = static_cast<double>(df) / 2.0;
  const double ib = betai(a, 0.5, x);
  if (t >= 0.0) {
    return 1.0 - 0.5 * ib;
  }
  return 0.5 * ib;
}

}  // namespace

double student_t_pvalue_two_sided(double t_stat, int df) {
  if (df <= 0) {
    return 1.0;
  }
  const double t_abs = std::abs(t_stat);
  return 2.0 * (1.0 - student_t_cdf(t_abs, df));
}

}  // namespace kmat
