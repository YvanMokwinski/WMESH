#pragma once

#include "MeshEntity.hpp"
#include "CRTP_MeshTopology.hpp"
template <DimensionType::enum_t _dimension> class MeshTopology;
#if 0

template <DimensionType::enum_t _dimension> struct CRTP_MeshTopology_traits<MeshTopology<_dimension> >
{
public: static constexpr bool IsMixed = true;
public: static constexpr DimensionType::enum_t Dimension = _dimension;
public: using meshentity_t 	= uint64_t;
public: using cell_t 		= uint64_t;
public: using cellidx_t 	= uint64_t;
public: using entitykind_t 	= typename ReferenceShapeTraits<_dimension>::type_t;
};
#endif
