#pragma once

template <class _derivedClass> struct Traits_CRTP_FiniteElementSpace
{
};

template <typename _derivedClass> class CRTP_FiniteElementSpace
{
  //
  // CRTP asImp() methods
  //
private: inline const _derivedClass & asImp() const { return static_cast<const _derivedClass&>(*this); }    
private: inline _derivedClass & asImp()  { return static_cast<_derivedClass&>(*this); }    
  
  //
  // Typedefs
  //
private: using traits_t = Traits_CRTP_FiniteElementSpace<_derivedClass>;
public: using meshtopology_t = typename traits_t::meshtopology_t;
public: using meshgeometry_t = typename traits_t::meshgeometry_t;
public: using fe_t = typename traits_t::fe_t;
  
  //
  // Delete copy constructor.
  //
public: CRTP_FiniteElementSpace(const CRTP_FiniteElementSpaceTopology<_derivedClass>&) = delete;
  
  //
  // Protected constructor.
  //
protected: inline CRTP_FiniteElementSpace(){};  

  //
  // The finite element.
  //  
public: const meshtopology_t * GetTopology() const noexcept
  {
    return asImp().GetTopology();
  };
  
public: const meshgeometry_t * GetGeometry() const noexcept
  {
    return asImp().GetGeometry();
  };

};


