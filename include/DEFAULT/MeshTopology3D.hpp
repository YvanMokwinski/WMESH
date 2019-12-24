#pragma once

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



  
  

template<> class MeshTopology<DimensionType::Volume>
  : public CRTP_MeshTopology<MeshTopology<DimensionType::Volume> >
{
  
protected: using this_t 	= MeshTopology<DimensionType::Volume>;
protected: using traits_t 	= CRTP_MeshTopology< this_t >::traits_t;  

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

public: MeshTopology(Input::Medit&inputMedit) noexcept
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

public: inline ~MeshTopology()
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

using MeshTopology3D = MeshTopology<DimensionType::Volume>;

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
inline unsigned int MeshTopology3D::GetEntityToEntities<DimensionType::Volume,DimensionType::Node>(const MeshTopology3D::entity_t<DimensionType::Volume>& sourceEntity_,
														       MeshTopology3D::entity_t<DimensionType::Node> targetEntities_[]) const noexcept
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
inline RangeMeshEntities<int_t> MeshTopology3D::GetEntities<DimensionType::Volume>(const VolumeType::enum_t entityKind_) const noexcept
{
  return RangeMeshEntities<int_t>(0,m_numCellsPerKind[entityKind_]);
};

template <>
inline RangeMeshEntities<int_t> MeshTopology3D::GetEntities<DimensionType::Node>(const NodeType::enum_t entityKind_) const noexcept
{
  return RangeMeshEntities<int_t>(0,m_numEntities[DimensionType::Node]);
};
