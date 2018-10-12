#include <iostream>
#include "Clothoid.hh"
#include "PolyLine.hh"

int main(int argc, const char * argv[]) {
  G2lib::ClothoidList curve;
  G2lib::PolyLine poly;
  int result = 14;

  curve.push_back_G1(0.0, 0.0, 0.0, 0.0, 3.0, G2lib::m_pi);
  std::cout << "Building polyline...";
  poly.build(curve, 0.01);
  std::cout << " done!" << std::endl;
  std::cout << "Checking result: ";
  if (poly.numSegment() ==result)
    std::cout << "CORRECT" << std::endl;
  else
    std::cout << "WRONG" << std::endl;
  
  return 0;
}