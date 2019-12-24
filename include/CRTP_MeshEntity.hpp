#pragma once

#include "DimensionType.hpp"

template <typename _derivedClass> struct CRTP_MeshEntityTraits;


template <typename _derivedClass> struct CRTP_MeshEntity
{
private: inline const _derivedClass & asImp() const { return static_cast<const _derivedClass&>(*this); }    
private: inline _derivedClass & asImp()  { return static_cast<_derivedClass&>(*this); }    

private: using traits_t = CRTP_MeshEntityTraits<_derivedClass>;

  //
  // Delete copy constructor.
  //
public: CRTP_MeshEntity(const CRTP_MeshEntity<_derivedClass>&) = delete;  

  
  //
  // Delete copy constructor.
  //
protected: inline CRTP_MeshEntity(){};

public: using celltype_t = typename traits_t::celltype_t;
public: using cellidx_t  = typename traits_t::cellidx_t;
  
public: inline celltype_t GetType() const noexcept
  {
    return asImp().GetType();
  };
  
public: inline cellidx_t GetIndex() const noexcept
  {
    return asImp().GetIndex();
  };
  
};
