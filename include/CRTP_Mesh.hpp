#pragma once
#include "CRTP_MeshTopology.hpp"
//!
//! @brief Traits for the class CRTP_Mesh.
//!
template <class _derivedClass> struct Traits_CRTP_Mesh
{
};

template <typename _derivedClass> class CRTP_Mesh
{
  //
  // CRTP asImp() methods
  //
private: inline const _derivedClass & asImp() const { return static_cast<const _derivedClass&>(*this); }    
private: inline _derivedClass & asImp()  { return static_cast<_derivedClass&>(*this); }    
  
  //
  // Typedefs
  //
private: using traits_t = Traits_CRTP_Mesh<_derivedClass>;
  
public:  using meshtopology_t = typename traits_t::meshtopology_t;
public:  using meshgeometry_t = typename traits_t::meshgeometry_t;
  
  //
  // Delete copy constructor.
  //
public: CRTP_Mesh(const CRTP_MeshTopology<_derivedClass>&) = delete;
  
  //
  // Protected constructor.
  //
protected: inline CRTP_Mesh(){};  
  
public: const meshtopology_t * GetTopology() const noexcept
  {
    return asImp().GetTopology();
  };

public: const meshgeometry_t * GetGeometry() const noexcept
  {
    return asImp().GetGeometry();
  };

};


