#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "accusphgeom/constructions/gca_constlat_intersection.hpp"
#include "test_helpers.hpp"

namespace {

using accusphgeom::constructions::gca_constlat_intersection;
using accusphgeom::numeric::Vec3;
using test_sigexp_helpers::parse_i64;
using test_sigexp_helpers::parse_scalar_sigexp;
using test_sigexp_helpers::parse_vec3_sigexp;
using test_sigexp_helpers::split_csv_line;

struct Row {
  std::int64_t case_id{};
  Vec3<double> a0{};
  Vec3<double> a1{};
  double z0{};
  double baseline_x{};
  double baseline_y{};
};

constexpr double kMachineEpsilon = 2.22e-16;
constexpr double kPoleTol = 3.0 * kMachineEpsilon;
constexpr double kEquatorTol = 100.0 * kMachineEpsilon;
constexpr double kPi = 3.141592653589793238462643383279502884;

double latitude_deg(const Vec3<double>& p) {
  return std::asin(p[2]) * 180.0 / kPi;
}

bool in_lat_bin(double latitude_deg, double lat_min_deg, double lat_max_deg) {
  return latitude_deg >= lat_min_deg && latitude_deg <= lat_max_deg;
}

double case_tolerance(const Row& row) {
  const double lat0 = latitude_deg(row.a0);
  const double lat1 = latitude_deg(row.a1);
  if (in_lat_bin(lat0, 0.0, 1.0) && in_lat_bin(lat1, 0.0, 1.0)) {
    return kEquatorTol;  // 100u near the equator.
  }
  return kPoleTol;  // 3u away from the equator.
}

std::vector<Row> read_rows(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open CSV: " + path);
  }

  std::string line;
  if (!std::getline(in, line)) {
    throw std::runtime_error("CSV is empty: " + path);
  }

  std::vector<Row> rows;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    const auto fields = split_csv_line(line);
    if (fields.size() != 19) {
      throw std::runtime_error("unexpected field count in CSV row");
    }

    Row row;
    row.case_id = parse_i64(fields[0]);
    row.a0 = parse_vec3_sigexp(fields, 1);
    row.a1 = parse_vec3_sigexp(fields, 7);
    row.z0 = parse_scalar_sigexp(fields, 13);
    row.baseline_x = parse_scalar_sigexp(fields, 15);
    row.baseline_y = parse_scalar_sigexp(fields, 17);
    rows.push_back(row);
  }
  return rows;
}

void print_failure_debug(const Row& row,
                         const Vec3<double>& result,
                         double err_xy,
                         double tol) {
  std::cerr << std::setprecision(17);
  std::cerr << "[FAIL] case_id=" << row.case_id << " err_xy=" << err_xy << '\n';
  std::cerr << "  a0=(" << row.a0[0] << ", " << row.a0[1] << ", " << row.a0[2]
            << ")\n";
  std::cerr << "  a1=(" << row.a1[0] << ", " << row.a1[1] << ", " << row.a1[2]
            << ")\n";
  std::cerr << "  z0=" << row.z0 << '\n';
  std::cerr << "  baseline=(" << row.baseline_x << ", " << row.baseline_y << ")\n";
  std::cerr << "  result=(" << result[0] << ", " << result[1] << ", "
            << result[2] << ")\n";
  std::cerr << "  tolerance=" << tol << '\n';
}

}  // namespace

int main() {
  const std::string csv_path =
      std::string(ACCUSPHGEOM_SOURCE_DIR) +
      "/tests/input/gca_constlat_cases_with_baseline.csv";
  const std::vector<Row> rows = read_rows(csv_path);

  if (rows.size() != 200) {
    std::cerr << "[FAIL] expected 200 rows, got " << rows.size() << '\n';
    return EXIT_FAILURE;
  }

  int failures = 0;
  for (const Row& row : rows) {
    try {
      const double tol = case_tolerance(row);
      const Vec3<double> result = gca_constlat_intersection(row.a0, row.a1, row.z0);
      const double dx = result[0] - row.baseline_x;
      const double dy = result[1] - row.baseline_y;
      const double err_xy = std::sqrt(dx * dx + dy * dy);
      if (err_xy >= tol) {
        ++failures;
        print_failure_debug(row, result, err_xy, tol);
      }
    } catch (const std::exception& ex) {
      ++failures;
      std::cerr << std::setprecision(17);
      std::cerr << "[FAIL] case_id=" << row.case_id
                << " threw exception: " << ex.what() << '\n';
      std::cerr << "  a0=(" << row.a0[0] << ", " << row.a0[1] << ", " << row.a0[2]
                << ")\n";
      std::cerr << "  a1=(" << row.a1[0] << ", " << row.a1[1] << ", " << row.a1[2]
                << ")\n";
      std::cerr << "  z0=" << row.z0 << '\n';
      std::cerr << "  baseline=(" << row.baseline_x << ", " << row.baseline_y
                << ")\n";
      try {
        std::cerr << "  tolerance=" << case_tolerance(row) << '\n';
      } catch (const std::exception&) {
      }
    }
  }

  if (failures != 0) {
    std::cerr << "[FAIL] " << failures << " / " << rows.size()
              << " baseline cases failed.\n";
    return EXIT_FAILURE;
  }

  std::cout << "[PASS] " << rows.size()
            << " gca_constlat_intersection baseline cases passed.\n";
  return EXIT_SUCCESS;
}
