#pragma once
//#include "ReferenceShapeTraits.hpp"
#include "ReferenceShape.hpp"
#include "ReferenceShapeFace.hpp"
#include "ReferenceShapeFaceTriangle.hpp"
#include "ReferenceShapeFaceQuadrilateral.hpp"
#include "ReferenceShapeVolume.hpp"
#include "ReferenceShapeVolumeTetrahedron.hpp"
#include "ReferenceShapeVolumePyramid.hpp"
#include "ReferenceShapeVolumeWedge.hpp"
#include "ReferenceShapeVolumeHexahedron.hpp"
#include "ReferenceShapeNode.hpp"
#include "ReferenceShapeEdge.hpp"

template <DimensionType::enum_t _dimension,
	  typename ReferenceShapeTraits<_dimension>::type_t::enum_t _cellType >
struct CellToNodes
{
  
public: using refshape_t = ReferenceShape<_dimension,_cellType>;  
public: std::array<int_t,refshape_t::NumEntities[DimensionType::Node]> m_nodeIds;
public: int_t m_cod;
  
};

using EdgeToNodes 		= CellToNodes<DimensionType::Edge, EdgeType::Edge >;
using TriangleToNodes 		= CellToNodes<DimensionType::Face, FaceType::Triangle >;
using QuadrilateralToNodes 	= CellToNodes<DimensionType::Face, FaceType::Quadrilateral >;
using TetrahedronToNodes 	= CellToNodes<DimensionType::Volume, VolumeType::Tetrahedron>;
using PyramidToNodes 		= CellToNodes<DimensionType::Volume, VolumeType::Pyramid>;
using WedgeToNodes 		= CellToNodes<DimensionType::Volume, VolumeType::Wedge>;
using HexahedronToNodes 	= CellToNodes<DimensionType::Volume, VolumeType::Hexahedron>;
