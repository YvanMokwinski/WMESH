#pragma once

#include "MeshEntity.hpp"

namespace AHF
{
  
  //!
  //! @brief The mesh node.
  //!
  template<> struct MeshEntity<DimensionType::Node>
  {
  private: using this_t = MeshEntity<DimensionType::Node>;
  private: int_t m_index;

  public: inline MeshEntity() noexcept
    {
    };
    
  public: inline MeshEntity(const int_t index_) noexcept
    : m_index(index_)
    {
    };
    
  public: inline int_t GetIndex() const noexcept
    {
      return this->m_index;
    };
    
  public: inline NodeType::enum_t GetKind() const noexcept
    {
      return NodeType::Node;
    };
    
  public: inline this_t& operator++() noexcept { ++this->m_index; return *this; };
  public: inline bool	 operator!=(const this_t& that_) const noexcept { return this->m_index != that_.m_index; };    
  };

};
