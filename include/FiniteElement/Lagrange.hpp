#pragma once

#include "DimensionType.hpp"
#include "ReferenceShapeTraits.hpp"

namespace FiniteElement
{
  
  template <unsigned int _degree,
	    DimensionType::enum_t _dimensionType,
	    typename ReferenceShapeTraits<_dimensionType>::type_t::enum_t _cellType >
  struct Lagrange;

};
