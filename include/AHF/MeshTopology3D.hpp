#pragma once

#include "CRTP_MeshTopology.hpp"
#include "AHF/MeshTopology.hpp"
#include "Input/Medit.hpp"

namespace AHF
{
  using MeshTopology3D = MeshTopology<DimensionType::Volume>;
};

//
// Specify traits
//
template<> struct CRTP_MeshTopology_traits<AHF::MeshTopology3D>
{
  //
  // Specify the topological dimension of the topology.
  //
public: static constexpr DimensionType::enum_t Dimension = DimensionType::Volume;
  //
  // Specify the type of the entity kind.
  //
public: using entitykind_t 	= VolumeType;  
public: template <DimensionType::enum_t _dimension> using entity_t = AHF::MeshEntity<_dimension>;
public: using celltype_t 	= VolumeType;
public: template <DimensionType::enum_t _dimension> using range_t = AHF::RangeMeshEntities<_dimension>;

};


//
//
//
namespace AHF
{

  //
  // Specialization of the mesh topology as a 3D mesh topology.
  //
  template <> class MeshTopology<DimensionType::Volume>
    : public CRTP_MeshTopology< MeshTopology<DimensionType::Volume> >
  {

    //
    // Definition of protected alias.
    //
  protected: using this_t 	= MeshTopology<DimensionType::Volume>;
  protected: using traits_t 	= CRTP_MeshTopology< this_t >::traits_t;    
  protected: using celltype_t 	= typename traits_t::celltype_t;
  protected: using cellkind_t 	= typename celltype_t::enum_t;
  protected: using cell_t =  typename traits_t::entity_t<DimensionType::Volume>;
  public: template <DimensionType::enum_t _dimension> using entity_t = traits_t::entity_t<_dimension>;
    

    //!
    //!
    //!
  public: MeshTopology(const int_t numNodes_, const int_t numCells_[]) noexcept;
    
    //!
    //!
    //!
  public: inline ~MeshTopology();
    
    //!
    //!
    //!
  protected: void ComputeNodesToCells() noexcept;

    //!
    //!
    //!
  public: inline DimensionType::enum_t 	GetDimension() 	const noexcept;
    
    //!
    //!
    //!
  public: inline int_t 	GetNumNodes() 	const noexcept;
    
    //!
    //!
    //!
  public: inline int_t 	GetNumCells() 	const noexcept;
    
    //!
    //!
    //!
  public: inline int_t 	GetNumCells(const cellkind_t cellKind_) const noexcept;

    //!
    //!
    //!
  public: template <cellkind_t _cellKind>
  inline int_t 		GetNumCells		() const noexcept;

    //!
    //!
    //!
  public: template <DimensionType::enum_t _dimension>
  inline int_t 		GetEntityIndex		(const entity_t<_dimension>& entity_) const noexcept;
    
    //!
    //!
    //!
  public: template <DimensionType::enum_t _dimension>
  inline range_t<_dimension> GetEntities(const typename ReferenceShapeTraits<_dimension>::type_t::enum_t entityKind_) const noexcept;

    //!
    //!
    //!
  public: template <DimensionType::enum_t _dimension>
  inline int_t 		GetNumEntities		(const typename ReferenceShapeTraits<_dimension>::type_t::enum_t cellType_) const noexcept;
    
    //!
    //!
    //!
  public: template <DimensionType::enum_t _dimension>
  inline int_t 		GetNumEntities		() const noexcept;

    //!
    //!
    //!
  public: template <DimensionType::enum_t _source,
		    DimensionType::enum_t _target>
  inline unsigned int 	GetEntityToEntities	(const entity_t<_source>& sourceEntity_,
						 entity_t<_target> targetEntities_[]) const noexcept;
    
  public: template <DimensionType::enum_t _target>
  inline unsigned int 	GetEntityToNodes	(const entity_t<Dimension>& 	cell_,
						 const unsigned int 		localEntityIndex_,
						 entity_t<DimensionType::Node> 	entityToNodes_[]) const noexcept;
    
  public: template<cellkind_t _cellKind,typename _array_int_t>
  inline void 		SetCellToNodes		(const int_t	cellIndex_,
						 _array_int_t&	cellToNodes_) noexcept;
    
  public: template <typename _array_int_t>
  inline unsigned int 	GetCellToNodes		(const cellkind_t 	cellKind_,
						 const int_t  	cellIndex_,
						 _array_int_t&	cellToNodes_) const noexcept;
    //
    // Array of size celltype_t::NumKinds
    //
  public: inline static MeshTopology3D * Build(Input::Medit&inputMedit);

    
  private:     	int_t 		m_numNodes;
  private:     	int_t 		m_numCells[celltype_t::NumKinds+1];    
  private:     	int_t 		m_beginCellsToNodes[celltype_t::NumKinds+1];
  private:     	int_t		m_sizeCellsToNodes;
  public:      	int_t*		m_cellsToNodes;
    
  private:     	int_t* 		m_beginNodesToCells;
  private: 	cell_t* 	m_nodesToCells;
  };

  //
  // Array of size celltype_t::NumKinds
  //
  inline MeshTopology3D * MeshTopology3D::Build(Input::Medit&inputMedit)
  {
    //
    // Initialize the number of cells.
    //
    int_t nu[4];
    for(const auto cellKind : celltype_t::All)
      {
	nu[cellKind] = inputMedit.GetNumVolumes(cellKind);
      }

    MeshTopology * meshTopology = new MeshTopology(inputMedit.GetNumNodes(),nu);
    meshTopology->m_beginCellsToNodes[0] = 0;
    for(const auto cellKind : celltype_t::All)
      {
	const int_t numCellsOfKind = meshTopology->m_numCells[cellKind];
	int_t size = 0;
	if (numCellsOfKind > 0)
	  {
	    size = numCellsOfKind * (1 + celltype_t::GetNumNodes(cellKind) );
	    inputMedit.ReadTopology(cellKind,
				    (int_t*)&meshTopology->m_cellsToNodes[meshTopology->m_beginCellsToNodes[cellKind]],
				    celltype_t::GetNumNodes(cellKind)+1,
				    (int_t*)&meshTopology->m_cellsToNodes[meshTopology->m_beginCellsToNodes[cellKind] + celltype_t::GetNumNodes(cellKind)],
				    celltype_t::GetNumNodes(cellKind)+1);
	  }
	meshTopology->m_beginCellsToNodes[cellKind+1] = meshTopology->m_beginCellsToNodes[cellKind] + size;
	  
      }

      
    std::cout << "compute nodes to cells" << std::endl;
    meshTopology->ComputeNodesToCells();
    std::cout << "compute nodes to cells done" << std::endl;
    return meshTopology;
  };
  
  
  template<MeshTopology3D::cellkind_t _cellKind,typename _array_int_t>
  void MeshTopology3D::SetCellToNodes(const int_t	cellIndex_,
				      _array_int_t&	cellToNodes_) noexcept
  {
    const unsigned int numNodesInCell = celltype_t::GetNumNodes(_cellKind);
    int_t at = this->m_beginCellsToNodes[_cellKind] + cellIndex_ * ( numNodesInCell + 1 );
    for (unsigned int i =0;i<numNodesInCell;++i)
      {
	this->m_cellsToNodes[at++] = MeshNode(cellToNodes_[i]).GetIndex();
      }
  };

    
  template <typename _array_int_t>
  inline unsigned int MeshTopology3D::GetCellToNodes(const cellkind_t 	cellKind_,
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


    

  
  inline DimensionType::enum_t MeshTopology3D::GetDimension() const noexcept
  {
    return traits_t::Dimension;
  };

  //
  // Get the entity index.
  //
  template <DimensionType::enum_t _dimension>
  inline int_t MeshTopology3D::GetEntityIndex(const entity_t<_dimension>& entity_) const noexcept
  {
    return entity_.GetIndex();
  };


  inline int_t MeshTopology3D::GetNumNodes() const noexcept
  {
    return this->m_numNodes;
  };
  
  inline int_t MeshTopology3D::GetNumCells() const noexcept
  {
    return this->m_numCells[celltype_t::NumKinds];
  };
  
  inline int_t MeshTopology3D::GetNumCells(const MeshTopology3D::cellkind_t cellKind_) const noexcept
  {
    return this->m_numCells[cellKind_];
  };
  
  template <MeshTopology3D::cellkind_t _cellKind> inline int_t MeshTopology3D::GetNumCells() const noexcept
  {
    return this->m_numCells[_cellKind];
  };

  MeshTopology3D::MeshTopology(const int_t numNodes_,
			       const int_t numCells_[]) noexcept
  : m_numNodes(numNodes_) 
  {
      
    //
    // 
    //
    this->m_numCells[celltype_t::NumKinds] = 0;
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
	    size = numCellsOfKind * (1 + celltype_t::GetNumNodes(cellKind) );
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
      //	std::cout << "size " << size << std::endl;
      this->m_sizeCellsToNodes = size;
      this->m_cellsToNodes = new int_t[size];
    }

  };

  inline MeshTopology3D::~MeshTopology()
  {    
    if(nullptr != this->m_cellsToNodes)
      {
	delete[] this->m_cellsToNodes;
	this->m_cellsToNodes = nullptr;
      }
  };





  void MeshTopology3D::ComputeNodesToCells() noexcept
  {
    //      std::cout << "compute ..." << std::endl;
    this->m_beginNodesToCells = new int_t[this->m_numNodes+1];
    for (int i=0;i<=m_numNodes;++i)
      {
	this->m_beginNodesToCells[i] = 0;
      }
    //      std::cout << "compute ..." << this->m_sizeCellsToNodes - m_numCells[celltype_t::NumKinds] * 1  << std::endl;
    this->m_nodesToCells 	= new cell_t[ this->m_sizeCellsToNodes - m_numCells[celltype_t::NumKinds] * 1 ];
    //      std::cout << "compute ..." << std::endl;
    this->m_beginNodesToCells[0] = 0;
    MeshNode cellToNodes[8];
    for(const auto cellKind : celltype_t::All)
      {
	//	  std::cout << "cellKind " << cellKind << std::endl;
	for(int_t cellIndex=0;cellIndex<this->m_numCells[cellKind];++cellIndex)
	  {
	    const auto numNodesInCell = this->GetCellToNodes(cellKind,
							     cellIndex,
							     cellToNodes);
	    // std::cout << "cellIndex " << cellIndex << std::endl;
	    for (unsigned int localNodeIndex = 0;localNodeIndex < numNodesInCell;++localNodeIndex)
	      {
		//		  std::cout << " " << cellToNodes[localNodeIndex].GetIndex()+1;
		this->m_beginNodesToCells[cellToNodes[localNodeIndex].GetIndex()+1] += 1;
	      }
	    //	      std::cout << std::endl;
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
		const auto index = cellToNodes[localNodeIndex].GetIndex();
		this->m_nodesToCells[this->m_beginNodesToCells[index]] = cell_t(cellIndex, cellKind);
		++this->m_beginNodesToCells[index];
	      }	      
	  }	  
      }
      
    for (int_t i=m_numNodes;i>0;--i)
      {
	this->m_beginNodesToCells[i] = this->m_beginNodesToCells[i-1];
      }      
    this->m_beginNodesToCells[0] = 0;
      
  };
  
  
  template <>
  inline unsigned int MeshTopology3D::GetEntityToEntities<DimensionType::Volume,DimensionType::Node>
  (const entity_t<DimensionType::Volume>& sourceEntity_,
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
  inline unsigned int MeshTopology3D::GetEntityToNodes<DimensionType::Edge>(const entity_t<DimensionType::Volume>& 	cell_,
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

  };
  
  
  template <>
  inline unsigned int MeshTopology3D::GetEntityToNodes<DimensionType::Face>(const entity_t<DimensionType::Volume>& 	cell_,
									    const unsigned int 			localEntityIndex_,
									    entity_t<DimensionType::Node> 		entityToNodes_[]) const noexcept
  {
    const auto cellKind 		= cell_.GetKind();
    const auto cellIndex 		= cell_.GetIndex();    
    const unsigned int numNodesInCell 	= celltype_t::GetNumNodes(cellKind);    
    const int_t at 			= this->m_beginCellsToNodes[cellKind] + cellIndex * ( numNodesInCell + 1 );
    
    const auto faceToNodes 		= celltype_t::GetFaceToNodes(cellKind,
								     localEntityIndex_);    
    
    const auto faceType 		= celltype_t::GetFaceType(cellKind,
								  localEntityIndex_);
    
    const auto numNodesInFace 		= FaceType::GetNumNodes(faceType);    
    for (unsigned int localNodeIndex = 0;localNodeIndex < numNodesInFace;++localNodeIndex)
      {
	entityToNodes_[localNodeIndex] = this->m_cellsToNodes[at + faceToNodes[localNodeIndex]];
      }
        
    return numNodesInFace;
  };

  
  template <> inline int_t MeshTopology3D::GetNumEntities<DimensionType::Node>() const noexcept
  {
    return this->GetNumNodes();
  };

  template <> inline int_t MeshTopology3D::GetNumEntities<DimensionType::Volume>() const noexcept
  {
    return this->GetNumCells();
  };
  
  
  template <>
  inline int_t MeshTopology3D::GetNumEntities<DimensionType::Volume>(const typename VolumeType::enum_t cellType_) const noexcept
  {
    return this->m_numCells[cellType_];
  };
  
  template <>
  inline RangeMeshEntities<DimensionType::Node> MeshTopology3D::GetEntities<DimensionType::Node>(const NodeType::enum_t entityKind_) const noexcept
  {
    return RangeMeshEntities<DimensionType::Node>(MeshNode(0),MeshNode(m_numNodes));
  };
  
  template <>
  inline RangeMeshEntities<DimensionType::Volume> MeshTopology3D::GetEntities<DimensionType::Volume>(const VolumeType::enum_t entityKind_) const noexcept
  {
    return RangeMeshEntities<DimensionType::Volume>(cell_t(0, entityKind_),
						    cell_t(m_numCells[entityKind_],entityKind_));
  };
  
};
