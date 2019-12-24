#pragma once

#include "CRTP_MeshTopology.hpp"

template <class _derivedClass> struct Traits_CRTP_MeshProperty
{

};

template <typename _derivedClass> class CRTP_MeshProperty
{
  //
  // CRTP asImp() methods
  //
private: inline const _derivedClass & asImp() const { return static_cast<const _derivedClass&>(*this); }    
private: inline _derivedClass & asImp()  { return static_cast<_derivedClass&>(*this); }    
  
private: using traits_t = Traits_CRTP_MeshProperty<_derivedClass>;  
public: static constexpr DimensionType::enum_t DimensionEntity = traits_t::DimensionEntity;
public: using meshtopology_t 	= typename traits_t::meshtopology_t;
public: using data_t 		= typename traits_t::data_t;
public: using entity_t 		= typename traits_t::meshtopology_t::template entity_t<DimensionEntity>;
  
  //
  // Delete copy constructor.
  //
public: CRTP_MeshProperty(const CRTP_MeshProperty<_derivedClass>&) = delete;
  
  //
  // Protected constructor.
  //
protected: inline CRTP_MeshProperty(){};  

public: inline const meshtopology_t* GetMeshTopology() const noexcept 
  {
    return asImp().GetMeshTopology();
  };
  
public: inline unsigned int GetNumComponents() const noexcept 
  {
    return asImp().GetNumComponents();
  };

public: template <typename _array_data_t> inline void GetData(const entity_t entity_,
							      _array_data_t& array_data_) const noexcept
  {
    asImp().GetData(entity_,
		    array_data_);
  };
  
};


