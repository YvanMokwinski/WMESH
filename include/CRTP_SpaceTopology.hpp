#pragma once

#include "MeshEntity.hpp"
#include "entity_iterator.hpp"


template <class _derivedClass> struct Traits_CRTP_GraphHandle
{
public: static constexpr bool IsMixed = false;
public: static constexpr DimensionType::enum_t Dimension = DimensionType::Node;
};


template <class _derivedClass> class CRTP_GraphHandle
{  
private: const _int_t  m_numElements;
  
public: CRTP_GraphHandle(const CRTP_GraphHandle<_derivedClass>&) = delete;

  //!
  //! @brief Get the number of entities for a specific dimension.
  //! @return The number of entities for a specific dimension.
  //!
public: template <DimensionType::enum_t _dimension> inline int_t GetNumEntities() const noexcept
  {
    return asImp().GetNumEntities<_dimension>();
  };

  //!
  //! @brief Get the number of entities for a specific dimension.
  //! @return The number of entities for a specific dimension.
  //!
public: template <DimensionType::enum_t _dimension> inline int_t GetNumEntities() const noexcept
  {
    return asImp().GetNumEntities<_dimension>();
  };

  
};

template <typename int_t> class GraphHandle : public CRTP_GraphHandle<GraphHandle>
{
  //!
  //! @brief Constant static cast of the derived class.
  //!  
private: inline const _derivedClass & asImp() const { return static_cast<const _derivedClass&>(*this); }    
  //!
  //! @brief Static cast of the derived class.
  //!  
private: inline _derivedClass & asImp()  { return static_cast<_derivedClass&>(*this); }    
  //!
  //! @brief Associated traits of the class.
  //!  
private: using traits_t = Traits_CRTP_GraphHandle<_derivedClass>;

  
private: const int_t  m_numElements;
private: const int_t* m_begin;
private: const int_t* m_indices;
  
public: constexpr SparsityPatternHandle(const int_t*b_,
					const int_t*i_)
  m_begin(b_),m_indices(i_)
  {
  };

  
};

template <class _derivedClass> struct Traits_CRTP_SpaceTopology
{
public: static constexpr bool IsMixed = false;
public: static constexpr DimensionType::enum_t Dimension = DimensionType::Node;
};

template <typename _derivedClass> class CRTP_SpaceTopology
{
  //!
  //! @brief Constant static cast of the derived class.
  //!  
private: inline const _derivedClass & asImp() const { return static_cast<const _derivedClass&>(*this); }    
  //!
  //! @brief Static cast of the derived class.
  //!  
private: inline _derivedClass & asImp()  { return static_cast<_derivedClass&>(*this); }    
  
  //
  // Typedefs
  //
public: using cellidx_t = int_t;

  //!
  //! @brief Associated traits of the class.
  //!  
private: using traits_t = Traits_CRTP_SpaceTopology<_derivedClass>;

  //
  // Type of the entity kind (tet,hex etc..)
  //
public: using entitykind_t 	= typename traits_t::entitykind_t;

  //!
  //! @brief Type of the mesh entity.
  //!  
public: template <DimensionType::enum_t _dimension>
using entity_t = typename traits_t::template entity_t<_dimension>;
  
  //!
  //! @brief Type of the range of entities.
  //!  
public: template <DimensionType::enum_t _dimension> using range_t = typename traits_t::template range_t<_dimension>;
  
  //!
  //! @brief The dimension as a constant expression.
  //!  
public: static constexpr DimensionType::enum_t Dimension = traits_t::Dimension;
  
public: using cell_t 	= entity_t<Dimension>;
public: using node_t 	= entity_t<DimensionType::Node>;
  
  //!
  //! @brief Deleted copy constructor.
  //!  
public: CRTP_SpaceTopology(const CRTP_SpaceTopology<_derivedClass>&) = delete;

  //!
  //! @brief Protected default constructor.
  //!  
protected: inline CRTP_SpaceTopology(){};  
  
  //!
  //! @brief Get the dimension of the topology.
  //! @return The dimension of the topology.
  //!
public: inline DimensionType::enum_t GetDimension() const noexcept
  {
    return asImp().GetDimension();
  };

  //!
  //! @brief Get the number of entities for a specific dimension.
  //! @return The number of entities for a specific dimension.
  //!
public: template <DimensionType::enum_t _dimension> inline int_t GetNumEntities() const noexcept
  {
    return asImp().GetNumEntities<_dimension>();
  };
  
  //!
  //! @brief Get the number of entities for a specific dimension and a specific type.
  //! @return The number of entities for a specific dimension and a specific type.
  //!
public: template <DimensionType::enum_t _dimension>
inline int_t GetNumEntities(const typename ReferenceShapeTraits<_dimension>::type_t::enum_t cellType_) const noexcept
  {
    return asImp().GetNumEntities<_dimension>(cellType_);
  };

  
  //!
  //! @brief Get the index of an entity of a specific dimension.
  //! @param entity_ the mesh entity.
  //! @return The index of the mesh entity.
  //!
public: template <DimensionType::enum_t _dimension>
inline int_t GetEntityIndex(const entity_t<_dimension>& entity_) const noexcept
  {
    return asImp().GetEntityIndex<_dimension>(entity_);
  };
  
  //!
  //! @brief Get the nodes of a specific entity of a cell.
  //! @param target_ the dimension of the target entity. 
  //! @param cell_ the cell 
  //! @param localEntityIndex_ the local entity index.
  //! @param entityToNodes_ the nodes.
  //! @return the number of nodes of the entity.
  //!
public: template <DimensionType::enum_t _target>
inline unsigned int GetEntityToDofs(const cell_t& 	cell_,
				    const unsigned int localEntityIndex_,
				    node_t 		entityToDofs_[]) const noexcept
  {
    return asImp().GetEntityToDofs<_target>(cell_,
					    localEntityIndex_,
					    entityToDofs_);
  };

  //!
  //! @brief Get the set of entities of a given dimension.
  //! @param entityKind_ The kind of the entity.
  //! @return The range.
  //!
public: template <DimensionType::enum_t _dimension>
inline range_t<_dimension> GetEntities(const typename ReferenceShapeTraits<_dimension>::type_t::enum_t entityKind_) const noexcept
  {
    return asImp().GetEntities<_dimension>(entityKind_);
  };
  
public: template <DimensionType::enum_t _source,
		  DimensionType::enum_t _target>
inline unsigned int GetEntityToEntities(const entity_t<_source>& sourceEntity_,
					entity_t<_target> targetEntities_[]) const noexcept
  {
    return asImp().GetEntityToEntities<_source,_target>(sourceEntity_,
							targetEntities_);
  };
  
};







#include "MeshTopology.hpp"
#if 0
#include "MeshProperty.hpp"
#endif



template <> struct CRTP_MeshTopology_traits<MeshTopology<DimensionType::Volume> >
{
public: static constexpr bool IsMixed = true;
public: static constexpr DimensionType::enum_t Dimension = DimensionType::Volume;
public: using cell_t 		= int_t;
public: using cellidx_t 	= int_t;
public: using entitykind_t 	= VolumeType;
public: template <DimensionType::enum_t _dimension> using range_t = RangeMeshEntities<int_t>;

public: template <DimensionType::enum_t _dimension> using entity_t = int_t;
  
  
#if 0
public:   	template <DimensionType::enum_t _dimension> using range_t = typename traits_t::range_t;
public: template <DimensionType::enum_t _dimension>
inline range_t<_dimension> GetEntities(const typename ReferenceShapeTraits<_dimension>::type_t entityKind_) const noexcept
  {
    return asImp().GetEntities<_dimension>(entityKind_);
  };
#endif
  
};



  
  

template<> class SpaceTopology<DimensionType::Volume>
  : public CRTP_SpaceTopology<SpaceTopology<DimensionType::Volume> >
{
  
protected: using this_t 	= SpaceTopology<DimensionType::Volume>;
protected: using traits_t 	= CRTP_SpaceTopology< this_t >::traits_t;  

protected: using cell_t 	= typename traits_t::cell_t;
protected: using cellidx_t 	= typename traits_t::cellidx_t;
public:    using entitykind_t 	= typename traits_t::entitykind_t;

public: template <DimensionType::enum_t _dimension> using entity_t = traits_t::entity_t<_dimension>;
  
  
public: static constexpr DimensionType::enum_t Dimension = traits_t::Dimension;
  
public:  
  
public:

public: template <DimensionType::enum_t _source,
		  DimensionType::enum_t _target>
inline unsigned int GetEntityToEntities(const entity_t<_source>& sourceEntity_,
					entity_t<_target> targetEntities_[]) const noexcept;
  
public: template <DimensionType::enum_t _dimension> inline int_t GetEntityIndex(const entity_t<_dimension>& entity_) const noexcept
  {
    return entity_;
  };
  
  template <VolumeType::enum_t _volumeType> struct RangeEntities
  {    
  protected: const this_t& m_meshTopology;
  public: using iterator = entity_iterator<cell_t>;
  public: inline iterator begin() 
    {
      return iterator(0);
    };
  public: inline iterator end() 
    {
      return iterator(m_meshTopology.m_numCellsPerKind[_volumeType]);     
    };
    
  public: inline RangeEntities(const this_t& meshTopology_)
    : m_meshTopology(meshTopology_)
    {
    };    
  };
  
public: template <DimensionType::enum_t _dimension>
inline range_t<_dimension> GetEntities(const typename ReferenceShapeTraits<_dimension>::type_t::enum_t entityKind_) const noexcept;

  
public: template <VolumeType::enum_t _volumeType> inline RangeEntities<_volumeType> GetRangeEntities()
  {
    return RangeEntities<_volumeType>(*this);
  };
  
public: inline DimensionType::enum_t GetDimension() const noexcept
  {
    return DimensionType::Volume;
  };


public: inline cellidx_t GetCellIndex(const cell_t cell_) const noexcept
  {
    return cell_;
  };
  
public: inline VolumeType::enum_t GetCellType(const cell_t cell_) const noexcept
  {
    return VolumeType::Tetrahedron;
  };

public: SpaceTopology(Input::Medit&inputMedit) noexcept
  {
    
    //
    // Initialize the number of cells.
    //
    for(const auto volumeKind : VolumeType::All)
      {
	this->m_numCellsPerKind[volumeKind] = inputMedit.GetNumVolumes(volumeKind);
	
	std::cout << "### "  << volumeKind << " " << this->m_numCellsPerKind[volumeKind] << std::endl;
      }

    //
    // The number of nodes
    //
    this->m_numEntities[0] = inputMedit.GetNumNodes();
    
    //
    // Calculate the total number of cells.
    //
    this->m_numEntities[3] = this->m_numCellsPerKind[0] + this->m_numCellsPerKind[1] + this->m_numCellsPerKind[2] + this->m_numCellsPerKind[3];

    //
    // Allocate Tetrahedron
    //
    if (this->m_numCellsPerKind[VolumeType::Tetrahedron]>0)
      {
	this->m_tetrahedraToNodes = new TetrahedronToNodes[this->m_numCellsPerKind[VolumeType::Tetrahedron]];
	this->m_cellsToNodes[VolumeType::Tetrahedron] = (int_t*)this->m_tetrahedraToNodes;
	
	inputMedit.ReadTopology(this->m_tetrahedraToNodes);	
#if 0
	for (int_t i = 0;i< this->m_numCellsPerKind[VolumeType::Tetrahedron];++i)
	  {
	    std::cout << this->m_tetrahedraToNodes[i].m_node0
		      << " " << this->m_tetrahedraToNodes[i].m_node1
		      << " " <<  this->m_tetrahedraToNodes[i].m_node2
		      << " " << this->m_tetrahedraToNodes[i].m_node3 << std::endl;
	  }
#endif
      }
    
    //
    // Allocate Hexahedron
    //
    if (this->m_numCellsPerKind[VolumeType::Hexahedron]>0)
      {
	this->m_hexahedraToNodes = new HexahedronToNodes[this->m_numCellsPerKind[VolumeType::Hexahedron]];
	this->m_cellsToNodes[VolumeType::Hexahedron] = (int_t*)this->m_hexahedraToNodes;
	inputMedit.ReadTopology(this->m_hexahedraToNodes);      
      }
    
    //
    // Allocate Pyramids
    //
    if (this->m_numCellsPerKind[VolumeType::Pyramid]>0)
      {
	this->m_pyramidsToNodes = new PyramidToNodes[this->m_numCellsPerKind[VolumeType::Pyramid]];
	this->m_cellsToNodes[VolumeType::Pyramid] = (int_t*)this->m_pyramidsToNodes;
	
	inputMedit.ReadTopology(this->m_pyramidsToNodes);
      }
    
    //
    // Allocate Wedges
    //
    if (this->m_numCellsPerKind[VolumeType::Wedge]>0)
      {
	this->m_wedgesToNodes = new WedgeToNodes[this->m_numCellsPerKind[VolumeType::Wedge]];
	this->m_cellsToNodes[VolumeType::Wedge] = (int_t*)this->m_wedgesToNodes;

	inputMedit.ReadTopology(this->m_wedgesToNodes);	
      }


    int j =0;
    for(const auto i : GetRangeEntities<VolumeType::Tetrahedron>())
      {
	j = (j<i)?i:j;

      }
    std::cout << "aaa " << j << std::endl;


    j =0;
    for(const auto i : GetRangeEntities<VolumeType::Hexahedron>())
      {
	j = (j<i)?i:j;
      }
    std::cout << "aaa " << j << std::endl;    
    
  };

public: inline ~SpaceTopology()
  {

    if(nullptr != this->m_wedgesToNodes)
      {
	delete[] this->m_wedgesToNodes;
	this->m_wedgesToNodes = nullptr;
      }

    if(nullptr != this->m_pyramidsToNodes)
      {
	delete[] this->m_pyramidsToNodes;
	this->m_pyramidsToNodes = nullptr;
      }

    if(nullptr != this->m_hexahedraToNodes)
      {
	delete[] this->m_hexahedraToNodes;
	this->m_hexahedraToNodes = nullptr;
      }
	
    if(nullptr != this->m_tetrahedraToNodes)
      {
	delete[] this->m_tetrahedraToNodes;
	this->m_tetrahedraToNodes = nullptr;
      }

  };

public: template <DimensionType::enum_t _dimension> inline int_t GetNumEntities() const noexcept
  {
    return this->m_numEntities[_dimension];
  };
  
public: template <DimensionType::enum_t _dimension>
inline int_t GetNumEntities(const typename ReferenceShapeTraits<_dimension>::type_t::enum_t cellType_) const noexcept
  {
    return this->m_numCellsPerKind[cellType_];
  };
  
private:
  
  int_t			m_numEntities[4]{0};  
  int_t   		m_numCellsPerKind[4]{0};
  int_t * 		m_cellsToNodes[4]{nullptr};
  int_t * 		m_cellCodes[4]{nullptr};
  HexahedronToNodes* 	m_hexahedraToNodes 	= nullptr;
  TetrahedronToNodes* 	m_tetrahedraToNodes 	= nullptr;
  WedgeToNodes* 	m_wedgesToNodes 	= nullptr;
  PyramidToNodes* 	m_pyramidsToNodes 	= nullptr;
  
};

using SpaceTopology3D = SpaceTopology<DimensionType::Volume>;

#if 0
public: template <DimensionType::enum_t _source,
		  DimensionType::enum_t _target>
inline unsigned int GetEntityToEntities(const entity_t<_source>& sourceEntity_,
					entity_t<_target> targetEntities_[]) const noexcept
  {
    return 0;
  };
#endif


template <>
inline unsigned int SpaceTopology3D::GetEntityToEntities<DimensionType::Volume,DimensionType::Node>(const SpaceTopology3D::entity_t<DimensionType::Volume>& sourceEntity_,
														       SpaceTopology3D::entity_t<DimensionType::Node> targetEntities_[]) const noexcept
{
  //
  // Get the cell
  //
  entitykind_t::enum_t cellKind_;
  if (sourceEntity_ < m_numCellsPerKind[0])
    {
      cellKind_ = entitykind_t::Tetrahedron;
    }
  else if (sourceEntity_ < m_numCellsPerKind[0] + m_numCellsPerKind[1])
    {
      cellKind_ = entitykind_t::Pyramid;
    }
  else if (sourceEntity_ < m_numCellsPerKind[0] + m_numCellsPerKind[1] + m_numCellsPerKind[2])
    {
      cellKind_ = entitykind_t::Wedge;

    }
  else
    {
      cellKind_ = entitykind_t::Hexahedron;
    }
  const unsigned int numNodesInCell = entitykind_t::GetNumNodes(cellKind_);
  const auto cellIndex = sourceEntity_;
  switch(cellKind_)
    {
    case entitykind_t::Tetrahedron:
      {
	for (unsigned int localNodeIndex = 0;localNodeIndex < numNodesInCell;++localNodeIndex)
	  {
	    targetEntities_[localNodeIndex] = m_tetrahedraToNodes[cellIndex].m_nodeIds[localNodeIndex];
	  }
	break;
      }
    case entitykind_t::Pyramid:
      {
	for (unsigned int localNodeIndex = 0;localNodeIndex < numNodesInCell;++localNodeIndex)
	  {
	    targetEntities_[localNodeIndex] = m_pyramidsToNodes[cellIndex].m_nodeIds[localNodeIndex];
	  }
	break;
      }
    case entitykind_t::Wedge:
      {
	for (unsigned int localNodeIndex = 0;localNodeIndex < numNodesInCell;++localNodeIndex)
	  {
	    targetEntities_[localNodeIndex] = m_wedgesToNodes[cellIndex].m_nodeIds[localNodeIndex];
	  }
	break;
      }
    case entitykind_t::Hexahedron:
      {
	for (unsigned int localNodeIndex = 0;localNodeIndex < numNodesInCell;++localNodeIndex)
	  {
	    targetEntities_[localNodeIndex] = m_hexahedraToNodes[cellIndex].m_nodeIds[localNodeIndex];
	  }
	break;
      }
    }
  return numNodesInCell;

};

template <>
inline RangeMeshEntities<int_t> SpaceTopology3D::GetEntities<DimensionType::Volume>(const VolumeType::enum_t entityKind_) const noexcept
{
  return RangeMeshEntities<int_t>(0,m_numCellsPerKind[entityKind_]);
};

template <>
inline RangeMeshEntities<int_t> SpaceTopology3D::GetEntities<DimensionType::Node>(const NodeType::enum_t entityKind_) const noexcept
{
  return RangeMeshEntities<int_t>(0,m_numEntities[DimensionType::Node]);
};
