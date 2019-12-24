#pragma once

#include "LagrangeFace.hpp"

namespace FiniteElement
{
  template <unsigned int _degree,VolumeType::enum_t _volumeType> using LagrangeVolume = Lagrange<_degree,
												 DimensionType::Volume,
												 _volumeType>;
  
};  
