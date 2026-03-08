#include "MultiPrecision_psm.h"

int main() {
  GEO::expansion_nt a(1.0), b(2.0);
  GEO::expansion_nt c = a*b - a;
  return (c.sign() == GEO::POSITIVE) ? 0 : 1;
}