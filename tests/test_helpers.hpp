#pragma once

#include <cmath>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "accusphgeom/constructions/gca_gca_intersection.hpp"

namespace test_sigexp_helpers {

using accusphgeom::numeric::Vec3;

inline double fromSigExp(std::int64_t sig, std::int64_t exp) {
  return std::ldexp(static_cast<double>(sig), static_cast<int>(exp));
}

inline std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, ',')) {
    fields.push_back(field);
  }
  return fields;
}

inline std::int64_t parse_i64(const std::string& text) {
  std::size_t pos = 0;
  const std::int64_t value = std::stoll(text, &pos);
  if (pos != text.size()) {
    throw std::runtime_error("invalid integer field: " + text);
  }
  return value;
}

inline double parse_scalar_sigexp(const std::vector<std::string>& fields,
                                  std::size_t start) {
  return fromSigExp(parse_i64(fields.at(start)), parse_i64(fields.at(start + 1)));
}

inline Vec3<double> parse_vec3_sigexp(const std::vector<std::string>& fields,
                                      std::size_t start) {
  return {
      parse_scalar_sigexp(fields, start),
      parse_scalar_sigexp(fields, start + 2),
      parse_scalar_sigexp(fields, start + 4),
  };
}

inline Vec3<double> subtract(const Vec3<double>& lhs, const Vec3<double>& rhs) {
  return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
}

inline double norm(const Vec3<double>& v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

}  // namespace test_sigexp_helpers
