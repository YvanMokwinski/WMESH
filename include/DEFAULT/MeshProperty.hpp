#pragma once

#include "CRTP_MeshProperty.hpp"


  template <typename _data_t,
	    unsigned int _numComponents,
	    DimensionType::enum_t _dimensionEntity>
  class MeshProperty;


template <typename _data_t,
	  unsigned int _numComponents,
	  DimensionType::enum_t _dimensionEntity>
struct Traits_CRTP_MeshProperty<MeshProperty<_data_t,_numComponents, _dimensionEntity> >
{
public: static constexpr DimensionType::enum_t DimensionEntity = _dimensionEntity;
public: using meshtopology_t 	= MeshTopology3D;
public: using data_t 		= _data_t;
};

  template <typename _data_t,
	    unsigned int _numComponents,
	    DimensionType::enum_t _dimensionEntity>
  class MeshProperty : public CRTP_MeshProperty<MeshProperty<_data_t,_numComponents, _dimensionEntity> >
  {
  private: using this_t = MeshProperty<_data_t,_numComponents,_dimensionEntity>;
  private: using traits_t = Traits_CRTP_MeshProperty<this_t>;  
  public: static constexpr DimensionType::enum_t DimensionEntity = traits_t::DimensionEntity;
  public: using meshtopology_t 	= typename traits_t::meshtopology_t;
  public: using data_t 		= typename traits_t::data_t;
  public: using entity_t	= typename traits_t::meshtopology_t::template entity_t<DimensionEntity>;
    
  public: MeshProperty(const meshtopology_t* meshtopology_)
    : m_meshtopology(meshtopology_),
      m_data(new _data_t[m_meshtopology->template GetNumEntities<_dimensionEntity>()])
    {      
    };
    
  public: inline const meshtopology_t* GetMeshTopology() const noexcept 
    {
      return this->m_meshtopology;
    };
    
  public: inline unsigned int GetNumComponents() const noexcept 
    {
      return _numComponents;
    };
    
  public: template <typename _array_data_t> inline void GetData(const entity_t entity_,
									_array_data_t& array_data_) const noexcept
    {
      const auto entityIndex = this->m_meshtopology->template GetEntityIndex<_dimensionEntity>(entity_);
      for (unsigned int componentIndex=0;componentIndex<_numComponents;++componentIndex)
	{
	  array_data_[componentIndex] = this->m_data[entityIndex * _numComponents + componentIndex];
	}
    };
    
  private: const meshtopology_t* m_meshtopology;
  private: _data_t* m_data;    

  };
