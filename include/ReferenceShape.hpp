#pragma once
#include "ReferenceShapeTraits.hpp"

template <DimensionType::enum_t _dimensionType, typename ReferenceShapeTraits<_dimensionType>::type_t::enum_t _cellType >
struct ReferenceShape;
