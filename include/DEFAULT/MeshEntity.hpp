#pragma once

#include "CRTP_MeshEntity.hpp"
#include "GenericEncoding.hpp"


template <DimensionType::enum_t _dimension> struct SpecificEncoding;

template <> struct SpecificEncoding<DimensionType::Volume>
{
public: static constexpr unsigned int NumDigitsForType = 2;
public: static constexpr unsigned int NumDigitsForRefinementLevel = 4;
public: static constexpr unsigned int NumDigits = NumDigitsForType + NumDigitsForRefinementLevel;  
public: using celltype_t = VolumeType::enum_t;
};

template <> struct SpecificEncoding<DimensionType::Face>
{
public: static constexpr unsigned int NumDigits = 1;
public: static constexpr unsigned int NumDigitsForRefinementLevel = 4;
public: using celltype_t = FaceType::enum_t;
};





template <DimensionType::enum_t _dimension> struct MeshEntity;
template <DimensionType::enum_t _dimension> struct CRTP_MeshEntityTraits< MeshEntity<_dimension> >
{
public: using cellidx_t = int_t;
public: using celltype_t = typename SpecificEncoding<_dimension>::celltype_t;
public: static constexpr unsigned int NumDigits = SpecificEncoding<_dimension>::NumDigits;
};


template <DimensionType::enum_t _dimension> struct MeshEntity : public CRTP_MeshEntity<MeshEntity<_dimension> >
{
private: using this_t = MeshEntity<_dimension>;
private: using traits_t = CRTP_MeshEntityTraits<this_t>;
private: using celltype_t = typename traits_t::celltype_t;  
private: using cellidx_t = typename traits_t::cellidx_t; 
private: using Encoder = GenericEncoding<int_t,traits_t::NumDigits>;
public: int_t m_encoding;
  
public: inline celltype_t GetType() const noexcept
  {
    return (celltype_t)Encoder::Low(m_encoding);
  };
  
public: inline cellidx_t GetIndex() const noexcept
  {
    return Encoder::Up(m_encoding);
  };
  
public: inline MeshEntity(const celltype_t cellType_,
			  const cellidx_t  cellIndex_)
  : m_encoding(Encoder::Encod(cellIndex_,(int_t)cellType_))
  {
    std::cout << "up "  << Encoder::s_upLimit << std::endl;
  };
  
};
