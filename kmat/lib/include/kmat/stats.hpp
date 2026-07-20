#pragma once

namespace kmat {

/// Two-sided p-value for Student's t with `df` degrees of freedom.
double student_t_pvalue_two_sided(double t_stat, int df);

}  // namespace kmat
