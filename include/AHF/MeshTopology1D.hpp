#pragma once

#include "CRTP_MeshTopology.hpp"
#include "AHF/MeshEntity.hpp"
#include "AHF/MeshTopology.hpp"

namespace AHF
{
  using MeshTopology1D = MeshTopology<DimensionType::Edge>;
};

//!
//! Implementation of the AHF 1D topology.
//!
template<> struct CRTP_MeshTopology_traits<AHF::MeshTopology1D>
{
  //
  // Specify the topological dimension of the topology.
  //
public: static constexpr DimensionType::enum_t Dimension = DimensionType::Edge;

  //
  // Specify the type of the entity kind.
  //
public: using entitykind_t 	= EdgeType;
  
public: template <DimensionType::enum_t _dimension> using entity_t = AHF::MeshEntity<_dimension>;
public: using celltype_t 	= EdgeType;
  //public: template <DimensionType::enum_t _dimension> using range_t = AHF::RangeMeshEntities<_dimension>;
  
};

namespace AHF
{
  
  template <> class MeshTopology<DimensionType::Edge>
    : public CRTP_MeshTopology< MeshTopology<DimensionType::Edge> >
  {    
  protected: using this_t 	= MeshTopology<DimensionType::Edge>;
  protected: using traits_t 	= CRTP_MeshTopology< this_t >::traits_t;    
  protected: using celltype_t 	= typename traits_t::celltype_t;
  protected: using cellkind_t 	= typename celltype_t::enum_t;
  public: template <DimensionType::enum_t _dimension> using entity_t = traits_t::entity_t<_dimension>;
    
  public: inline DimensionType::enum_t GetDimension() const noexcept
    {
      return traits_t::Dimension;
    };

  public: template <DimensionType::enum_t _dimension> inline int_t GetEntityIndex(const entity_t<_dimension>& entity_) const noexcept
    {
      return entity_.GetIndex();
    };
    
  public: inline int_t GetNumNodes() const noexcept
    {
      return this->m_numNodes;
    };
    
  public: inline int_t GetNumCells() const noexcept
    {
      return this->m_numEdges;
    };
    
  public: inline int_t GetNumCells(const cellkind_t cellKind_) const noexcept
    {
      return this->m_numEdges;
    };

  public: template <cellkind_t _cellKind> inline int_t GetNumCells() const noexcept
    {
      return this->m_numEdges;
    };


  public: template<cellkind_t _cellKind,typename _array_int_t> void SetCellToNodes(const int_t 		cellIndex_,
										   _array_int_t&	cellToNodes_) noexcept
    {
      const unsigned int numNodesInCell = celltype_t::GetNumNodes(_cellKind);
      int_t at = this->m_beginCellsToNodes[_cellKind] + cellIndex_ * ( numNodesInCell + 1 );
      for (unsigned int i =0;i<numNodesInCell;++i)
	{
	  this->m_cellsToNodes[at++] = MeshNode(cellToNodes_[i]);
	}
    };
    
  public: template <typename _array_int_t>
  inline unsigned int GetCellToNodes(const cellkind_t 	cellKind_,
				     const int_t  	cellIndex_,
				     _array_int_t&	cellToNodes_) const noexcept
    {
      const unsigned int numNodesInCell = celltype_t::GetNumNodes(cellKind_);
      int_t at = this->m_beginCellsToNodes[cellKind_] + cellIndex_ * ( numNodesInCell + 1 );
      for (unsigned int i =0;i<numNodesInCell;++i)
	{
	  cellToNodes_[i] = this->m_cellsToNodes[at++];
	}
      return numNodesInCell;
    };
    
  public: template <DimensionType::enum_t _dimension>
  inline int_t GetNumEntities(const typename ReferenceShapeTraits<_dimension>::type_t::enum_t cellType_) const noexcept;
    
  public: template <DimensionType::enum_t _dimension> inline int_t GetNumEntities() const noexcept;
#if 0
  public: template <DimensionType::enum_t _dimension>
  inline range_t<_dimension> GetEntities(const typename ReferenceShapeTraits<_dimension>::type_t::enum_t entityKind_) const noexcept;
#endif
    
  public: template <DimensionType::enum_t _source,
		    DimensionType::enum_t _target>
  inline unsigned int GetEntityToEntities	(const entity_t<_source>& sourceEntity_,
						 entity_t<_target> targetEntities_[]) const noexcept;
    
  public: template <DimensionType::enum_t _target>
  inline unsigned int GetEntityToNodes		(const entity_t<Dimension>& 	cell_,
						 const unsigned int 		localEntityIndex_,
						 entity_t<DimensionType::Node> 	entityToNodes_[]) const noexcept;

  public: MeshTopology(const int_t numNodes_,
		       const int_t numCells_[]) noexcept
  : m_numNodes(numNodes_)      
    {
      
      for(const auto cellKind : celltype_t::All)
	{	 
	  this->m_numCells[cellKind] = numCells_[cellKind];
	  this->m_numCells[celltype_t::NumKinds] += numCells_[cellKind];
	}

      
      this->m_beginCellsToNodes[0] = 0;
      for(const auto cellKind : celltype_t::All)
	{
	  const int_t numCellsOfKind = this->m_numCells[cellKind];
	  int_t size = 0;
	  if (numCellsOfKind > 0)
	    {
	      size = numCellsOfKind * celltype_t::GetNumNodes(cellKind);
	    }
	  this->m_beginCellsToNodes[cellKind+1] = this->m_beginCellsToNodes[cellKind] + size;
	}

      {
	int_t size = 0;
	for(const auto cellKind : celltype_t::All)
	  {
	    size += 
	      ( 1 + celltype_t::GetNumNodes(cellKind) ) * this->m_numCells[cellKind];
	  }
	std::cout << "size " << size << std::endl;
	this->m_sizeCellsToNodes = size;
	this->m_cellsToNodes = new MeshNode[size];
      }

    };
    
  public: MeshTopology(Input::Medit&inputMedit) noexcept
    {
      this->m_numNodes = inputMedit.GetNumNodes();
      //
      // Initialize the number of cells.
      //
      for(const auto cellKind : celltype_t::All)
	{
	  this->m_numCells[cellKind] = inputMedit.GetNumEdges(cellKind);
	}

      this->m_numCells[celltype_t::NumKinds] = 0;
      for(const auto cellKind : celltype_t::All)
	{
	  this->m_numCells[celltype_t::NumKinds] += this->m_numCells[cellKind];
	}
      
      {
	int_t size = 0;
	for(const auto cellKind : celltype_t::All)
	  {
	    size += 
	      ( 1 + celltype_t::GetNumNodes(cellKind) ) * this->m_numCells[cellKind];
	  }
	std::cout << "size " << size << std::endl;
	this->m_sizeCellsToNodes = size;
	this->m_cellsToNodes = new MeshNode[size];
      }

      this->m_beginCellsToNodes[0] = 0;
      for(const auto cellKind : celltype_t::All)
	{
	  const int_t numCellsOfKind = this->m_numCells[cellKind];
	  int_t size = 0;
	  if (numCellsOfKind > 0)
	    {
	      size = numCellsOfKind * celltype_t::GetNumNodes(cellKind);
	      inputMedit.ReadTopology(cellKind,
				      (int_t*)&this->m_cellsToNodes[this->m_beginCellsToNodes[cellKind]],
				      celltype_t::GetNumNodes(cellKind)+1,
				      (int_t*)&this->m_cellsToNodes[this->m_beginCellsToNodes[cellKind] + celltype_t::GetNumNodes(cellKind)],
				      celltype_t::GetNumNodes(cellKind)+1);
	      
	    }
	  this->m_beginCellsToNodes[cellKind+1] = this->m_beginCellsToNodes[cellKind] + size;
	}
      std::cout << "compute" << std::endl;
      this->ComputeNodesToCells();
      std::cout << "compute done" << std::endl;
    };
    
  public: inline ~MeshTopology()
    {    
      if(nullptr != this->m_cellsToNodes)
	{
	  delete[] this->m_cellsToNodes;
	  this->m_cellsToNodes = nullptr;
	}
    };

    
  protected: void ComputeNodesToCells() noexcept
    {
      std::cout << "compute ..." << std::endl;
      this->m_beginNodesToCells = new int_t[this->m_numNodes+1];
      std::cout << "compute ..." << this->m_sizeCellsToNodes - m_numCells[celltype_t::NumKinds] * 1  << std::endl;
      this->m_nodesToCells 	= new MeshCell[ this->m_sizeCellsToNodes - m_numCells[celltype_t::NumKinds] * 1 ];
      std::cout << "compute ..." << std::endl;

      this->m_beginNodesToCells[0] = 0;
      MeshNode cellToNodes[8];
      for(const auto cellKind : celltype_t::All)
	{
	  
	  for(int_t cellIndex=0;cellIndex<this->m_numCells[cellKind];++cellIndex)
	    {
	      const auto numNodesInCell = this->GetCellToNodes(cellKind,
							       cellIndex,
							       cellToNodes);
	      
	      for (unsigned int localNodeIndex = 0;localNodeIndex < numNodesInCell;++localNodeIndex)
		{
		  this->m_beginNodesToCells[cellToNodes[localNodeIndex].GetIndex()+1] += 1;
		}	      
	    }
	  
	}
      
      for (int_t i=2;i<=m_numNodes;++i)
	{
	  this->m_beginNodesToCells[i] += this->m_beginNodesToCells[i-1];
	}
      
      for(const auto cellKind : celltype_t::All)
	{
	  
	  for(int_t cellIndex=0;cellIndex<this->m_numCells[cellKind];++cellIndex)
	    {
	      const auto numNodesInCell = this->GetCellToNodes(cellKind,
							       cellIndex,
							       cellToNodes);
	      
	      for (unsigned int localNodeIndex = 0;localNodeIndex < numNodesInCell;++localNodeIndex)
		{
		  this->m_nodesToCells[this->m_beginNodesToCells[cellToNodes[localNodeIndex].GetIndex()]] = MeshCell(cellIndex, cellKind);
		  ++this->m_beginNodesToCells[cellToNodes[localNodeIndex].GetIndex()];
		}	      
	    }	  
	}

      
      for (int_t i=m_numNodes;i>0;--i)
	{
	  this->m_beginNodesToCells[i] = this->m_beginNodesToCells[i-1];
	}      
      this->m_beginNodesToCells[0] = 0;
      
    };
    
  private:
    int_t 	m_numNodes;
    int_t 	m_numCells[celltype_t::NumKinds+1];    
    int_t 	m_beginCellsToNodes[celltype_t::NumKinds+1];
    int_t	m_sizeCellsToNodes;
    MeshNode*	m_cellsToNodes;
    int_t* 	m_beginNodesToCells;
    MeshCell* 	m_nodesToCells;
  };
    
  template <>
  inline unsigned int MeshTopology1D::GetEntityToEntities<DimensionType::Edge,DimensionType::Node>
  (const entity_t<DimensionType::Edge>& sourceEntity_,
   entity_t<DimensionType::Node> targetEntities_[]) const noexcept
  {
    const auto cellKind_ = sourceEntity_.GetKind();
    const auto cellIndex_ = sourceEntity_.GetIndex();
    const unsigned int numNodesInCell = celltype_t::GetNumNodes(cellKind_);
    int_t at = this->m_beginCellsToNodes[cellKind_] + cellIndex_ * ( numNodesInCell + 1 );
    for (unsigned int i =0;i<numNodesInCell;++i)
      {
	targetEntities_[i] = this->m_cellsToNodes[at++];
      }
    return numNodesInCell;
  };

  template <>
  inline unsigned int MeshTopology1D::GetEntityToNodes<DimensionType::Edge>(const entity_t<DimensionType::Volume>& 	cell_,
									    const unsigned int 				localEntityIndex_,
									    entity_t<DimensionType::Node> 		entityToNodes_[]) const noexcept
  {
    const auto cellKind 	= cell_.GetKind();
    const auto cellIndex 	= cell_.GetIndex();
    
    const unsigned int numNodesInCell = celltype_t::GetNumNodes(cellKind);    
    const int_t at = this->m_beginCellsToNodes[cellKind] + cellIndex * ( numNodesInCell + 1 );

    
    const auto edgeToNodes = celltype_t::GetEdgeToNodes(cellKind,
							localEntityIndex_);    
    entityToNodes_[0] = this->m_cellsToNodes[at + edgeToNodes[0]];
    entityToNodes_[1] = this->m_cellsToNodes[at + edgeToNodes[1]];
    return 2;
    
#if 0
    for (const auto localNodeIndex : VolumeType::EdgeToNodes[cellKind][localEntityIndex_])
      {
	entityToNodes_[localNodeIndex] = this->m_cellsToNodes[at + localNodeIndex];
      }
#endif

  };
  
  template <> inline int_t MeshTopology1D::GetNumEntities<DimensionType::Node>() const noexcept
  {
    return this->GetNumNodes();
  };

  template <> inline int_t MeshTopology1D::GetNumEntities<DimensionType::Edge>() const noexcept
  {
    return this->GetNumCells();
  };
  
  
  template <>
  inline int_t MeshTopology1D::GetNumEntities<DimensionType::Edge>(const typename VolumeType::enum_t cellType_) const noexcept
  {
    return this->m_numCells[cellType_];
  };
  
  template <>
  inline RangeMeshEntities<DimensionType::Node> MeshTopology1D::GetEntities<DimensionType::Node>(const NodeType::enum_t entityKind_) const noexcept
  {
    return RangeMeshEntities<DimensionType::Node>(MeshNode(0),MeshNode(m_numNodes));
  };
  
  template <>
  inline RangeMeshEntities<DimensionType::Edge> MeshTopology1D::GetEntities<DimensionType::Edge>(const VolumeType::enum_t entityKind_) const noexcept
  {
    return RangeMeshEntities<DimensionType::Edge>(MeshCell(0,entityKind_),
						  MeshCell(m_numCells[entityKind_],entityKind_));
  };
  
};
