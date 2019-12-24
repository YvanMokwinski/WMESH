#pragma once

template <typename _impl_t> class CRTP
{
  //!
  //! @brief Constant static cast of the derived class.
  //!  
protected: inline const _impl_t& 	asImp() const;  
  //!
  //! @brief Static cast of the derived class.
  //!  
protected: inline _impl_t& 	asImp();

};


//
//
//
template <typename _impl_t> 
inline const _impl_t & CRTP<_impl_t>::asImp() const
{ return static_cast<const _impl_t&>(*this); };

//
//
//
template <typename _impl_t> 
inline _impl_t & CRTP<_impl_t>::asImp()
{ return static_cast<_impl_t&>(*this); };
