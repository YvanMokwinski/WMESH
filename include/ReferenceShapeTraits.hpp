#pragma once

#include "DimensionType.hpp"
#include "NodeType.hpp"
#include "EdgeType.hpp"
#include "FaceType.hpp"
#include "VolumeType.hpp"


template <DimensionType::enum_t _dimensionType> struct ReferenceShapeTraits;
  
template <> struct ReferenceShapeTraits<DimensionType::Volume>
{
public: using type_t = VolumeType;
};
  
template <> struct ReferenceShapeTraits<DimensionType::Face>
{
public: using type_t = FaceType;
};
  
template <> struct ReferenceShapeTraits<DimensionType::Edge>
{
public: using type_t = EdgeType;
};
  
template <> struct ReferenceShapeTraits<DimensionType::Node>
{
public: using type_t = NodeType;
};
