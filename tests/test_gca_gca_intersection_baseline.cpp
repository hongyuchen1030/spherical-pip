#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "accusphgeom/constructions/gca_gca_intersection.hpp"

namespace {

using accusphgeom::constructions::gca_gca_intersection;
using accusphgeom::numeric::Vec3;

struct Row {
  std::int64_t pair_id{};
  Vec3<double> a0{};
  Vec3<double> a1{};
  Vec3<double> b0{};
  Vec3<double> b1{};
  Vec3<double> baseline{};
};

inline double fromSigExp(std::int64_t sig, std::int64_t exp) {
  return std::ldexp(static_cast<double>(sig), static_cast<int>(exp));
}

std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, ',')) {
    fields.push_back(field);
  }
  return fields;
}

std::int64_t parse_i64(const std::string& text) {
  std::size_t pos = 0;
  const std::int64_t value = std::stoll(text, &pos);
  if (pos != text.size()) {
    throw std::runtime_error("invalid integer field: " + text);
  }
  return value;
}

Vec3<double> parse_vec3_sigexp(const std::vector<std::string>& fields, std::size_t start) {
  return {
      fromSigExp(parse_i64(fields.at(start)), parse_i64(fields.at(start + 1))),
      fromSigExp(parse_i64(fields.at(start + 2)), parse_i64(fields.at(start + 3))),
      fromSigExp(parse_i64(fields.at(start + 4)), parse_i64(fields.at(start + 5))),
  };
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
    if (fields.size() != 32) {
      throw std::runtime_error("unexpected field count in CSV row");
    }

    Row row;
    row.pair_id = parse_i64(fields[0]);
    row.a0 = parse_vec3_sigexp(fields, 2);
    row.a1 = parse_vec3_sigexp(fields, 8);
    row.b0 = parse_vec3_sigexp(fields, 14);
    row.b1 = parse_vec3_sigexp(fields, 20);
    row.baseline = parse_vec3_sigexp(fields, 26);
    rows.push_back(row);
  }
  return rows;
}

Vec3<double> subtract(const Vec3<double>& lhs, const Vec3<double>& rhs) {
  return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
}

double norm(const Vec3<double>& v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

}  // namespace

int main() {
  const std::string csv_path =
      std::string(ACCUSPHGEOM_SOURCE_DIR) +
      "/tests/input/gca_gca_pairs_seed20251104_N100_with_baseline.csv";
  const std::vector<Row> rows = read_rows(csv_path);
  const double tol = 1e-8;

  if (rows.size() != 31) {
    std::cerr << "[FAIL] expected 31 rows, got " << rows.size() << '\n';
    return EXIT_FAILURE;
  }

  int failures = 0;
  for (const Row& row : rows) {
    try {
      const Vec3<double> result =
          gca_gca_intersection(row.a0, row.a1, row.b0, row.b1);
      const double err = norm(subtract(result, row.baseline));
      if (!(err < tol)) {
        ++failures;
        std::cerr << "[FAIL] pair_id=" << row.pair_id
                  << " error=" << err << '\n';
      }
    } catch (const std::exception& ex) {
      ++failures;
      std::cerr << "[FAIL] pair_id=" << row.pair_id
                << " threw exception: " << ex.what() << '\n';
    }
  }

  if (failures != 0) {
    std::cerr << "[FAIL] " << failures << " / " << rows.size()
              << " baseline cases failed.\n";
    return EXIT_FAILURE;
  }

  std::cout << "[PASS] " << rows.size()
            << " gca_gca_intersection baseline cases passed.\n";
  return EXIT_SUCCESS;
}
