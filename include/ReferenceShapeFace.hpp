#pragma once
#include "ReferenceShapeEdge.hpp"

template <FaceType::enum_t _faceType> using ReferenceShapeFace = ReferenceShape<DimensionType::Face,_faceType>;
