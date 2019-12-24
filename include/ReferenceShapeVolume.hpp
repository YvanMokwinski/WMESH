#pragma once
#include "ReferenceShapeFace.hpp"

template <VolumeType::enum_t _volumeType> using ReferenceShapeVolume = ReferenceShape<DimensionType::Volume,_volumeType>;
