#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "kmat/stats.hpp"

TEST_CASE("student_t_pvalue_two_sided sane at small df", "[stats]") {
  using Catch::Matchers::WithinAbs;

  const double p_center = kmat::student_t_pvalue_two_sided(0.0, 2);
  REQUIRE_THAT(p_center, WithinAbs(1.0, 1e-6));

  const double p_tail = kmat::student_t_pvalue_two_sided(10.0, 2);
  REQUIRE(p_tail < 0.01);
  REQUIRE(p_tail > 0.0);
}
