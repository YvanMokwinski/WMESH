#pragma once

template <class _derivedClass> struct Traits_CRTP_Space
{
};

template <typename _derivedClass> class CRTP_Space
{
  //
  // CRTP asImp() methods
  //
private: inline const _derivedClass & asImp() const { return static_cast<const _derivedClass&>(*this); }    
private: inline _derivedClass & asImp()  { return static_cast<_derivedClass&>(*this); }    
  
  //
  // Typedefs
  //
private: using traits_t = Traits_CRTP_Space<_derivedClass>;
  
public: using fe_t 		= typename traits_t::fe_t;
public: using spacetopology_t 	= typename traits_t::spacetopology_t;
public: using spacegeometry_t 	= typename traits_t::spacegeometry_t;
  
  //
  // Delete copy constructor.
  //
public: CRTP_Space(const CRTP_SpaceTopology<_derivedClass>&) = delete;
  
  //
  // Protected constructor.
  //
protected: inline CRTP_Space(){};  

  //
  // The finite element.
  //  
public: const spacetopology_t * GetTopology() const noexcept
  {
    return asImp().GetTopology();
  };
  
public: const spacegeometry_t * GetGeometry() const noexcept
  {
    return asImp().GetGeometry();
  };

};


