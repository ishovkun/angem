#pragma once

namespace angem {

enum VTK_ID : int
{
  VertexID = 1,
  PolyVertexID = 2,
  LineID = 3,
  PolyLineID = 4,
  TriangleID = 5,
  TriangleStripID = 6,
  GeneralPolygonID = 7,
  PixelID = 8,
  QuadrangleID  = 9,
  TetrahedronID = 10,
  VoxelID = 11,
  HexahedronID = 12,
  WedgeID = 13,
  PyramidID = 14,
  PentagonalPrismID = 15,
  HexagonalPrismID = 16,
  QuadraticEdgeID = 21,
  QuadraticTriangleID = 22,
  QuadraticQuadID = 23,
  QuadraticTetraID = 24,
  QudraticHexahedronID = 25,
  QuadraticWedgeID = 26,
  QuadraticPyramidID = 27,
  BiquadraticQuadID = 28,
  TriquadraticHexahedronID = 29,
  QuadraticLinearQuadID = 30,
  QuadraticLinearWedgeID = 31,
  BiquadraticWedgeID = 32,
  BiquadraticQuadraticHexahedronID = 33,
  BiquadraticTriangleID = 34,
  CubicLineID = 35,
  QuadraticPolygonID = 36,
  GeneralPolyhedronID = 42,
  InvalidElementID = -1
};

}  // end namespace angem
