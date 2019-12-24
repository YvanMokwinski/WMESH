#pragma once

#include "MeshEntity.hpp"

namespace AHF
{
  //!
  //! @brief AHF representation of a 3D-entity.
  //!
  //! @remark
  //! The type of the cell is encoded in the mesh entity (2 bits)
  //! is enough to enumerate the 4 classical volumes.
  //!
  template<> struct MeshEntity<DimensionType::Volume>
  {
  private: int_t m_id;
  private: using this_t = MeshEntity<DimensionType::Volume>;

  public: inline MeshEntity() noexcept
    {
    };
    
  public: inline MeshEntity(const int_t index_,
			    const VolumeType::enum_t volumeKind_) noexcept
    : m_id(GenericEncoding<int_t,2>::Encod(index_,(int_t)volumeKind_))
    {
    };
    
  public: inline int_t GetIndex() const noexcept
    {
      return GenericEncoding<int_t,2>::Up(this->m_id);
    };
    
  public: inline VolumeType::enum_t GetKind() const noexcept
    {
      return (VolumeType::enum_t)GenericEncoding<int_t,2>::Low(this->m_id);
    };
    
  public: inline this_t& 	operator++() noexcept
    {
      ++this->m_id;
      return *this;
    };
    
  public: inline bool		operator!=(const this_t& that_) const noexcept
    {
      return this->m_id != that_.m_id;
    };
    
  };
  

  using MeshNode = MeshEntity<DimensionType::Node>;

  //
  // CECI EST UNE CONNERIE.
  //
  
  template <DimensionType::enum_t _dimension> struct ahfentity_iterator : public std::iterator<std::forward_iterator_tag, int_t>
  {
  MeshEntity<_dimension> m_val;
public:
  ahfentity_iterator(const MeshEntity<_dimension> v) : m_val(v) {};
  ~ahfentity_iterator() {};
  inline ahfentity_iterator& 		operator++() noexcept { ++this->m_val; return *this; };
  inline MeshEntity<_dimension> 	operator* () const  noexcept  { return this->m_val; };
  inline bool  			operator!=(const ahfentity_iterator& rhs) const noexcept { return m_val != rhs.m_val; };
};
  
  template <DimensionType::enum_t _dimension> struct RangeMeshEntities
  {
  protected: MeshEntity<_dimension> m_start;
  protected: MeshEntity<_dimension> m_bound;
    //  protected: _int_t m_start;
    //  protected: _int_t m_bound;
  public: using iterator = ahfentity_iterator<_dimension>;

  public: inline iterator begin() 
    {
      return iterator(m_start);
    };
    
  public: inline iterator end() 
    {
      return iterator(m_bound);     
    };
    
  public: inline RangeMeshEntities(const MeshEntity<_dimension> start_,
				   const MeshEntity<_dimension> bound_)
    : m_start(start_), m_bound(bound_)
    {
    };
    
  };

};

#if 0
template <DimensionType::enum_t _dimension>
std::ostream& operator<< (std::ostream &out_,
			  const AHF::MeshEntity<_dimension>& meshEntity_)
{
  out_ << *((const int_t*)&meshEntity_);
  return out_;
};
#endif
