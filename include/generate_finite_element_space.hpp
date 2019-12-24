#pragma once

#include "Config.hpp"

inline void get_cell_to_nodes(unsigned int 		numCellNodes_,
			      const int_t*__restrict__ 	cellsToNodes_,
			      int_t 			cellsToNodesLd_,
			      int_t 			cellIndex_,
			      int_t*__restrict__ 	cnc_)
{
  for (unsigned int localNodeIndex=0;localNodeIndex<numCellNodes_;++localNodeIndex)
    {
      cnc_[localNodeIndex] = cellsToNodes_[cellIndex_*cellsToNodesLd_+localNodeIndex];      
    }
}

inline void HalfFaceDecomposition(const int_t halfFaceIndex_,
				  const int_t numDofsPerCell_,
				  const int_t numDofsPerInteriorFace_,
				  int_t* cellIndex_,
				  int_t* localFaceIndex_)
{
  cellIndex_[0] 	= halfFaceIndex_ / numDofsPerCell_;
  localFaceIndex_[0] 	= ( halfFaceIndex_ % numDofsPerCell_ ) / numDofsPerInteriorFace_;
#ifndef NDEBUG
  if (localFaceIndex_[0]<0 || localFaceIndex_[0] >=4)
    {
      std::cerr << "error"<< std::endl;
      exit(1);
    }
#endif
};

inline void HalfEdgeDecomposition(const int_t halfEdgeIndex_,
				  const int_t numDofsPerCell_,
				  const int_t numDofsPerInteriorEdge_,
				  int_t* cellIndex_,
				  int_t* localEdgeIndex_)
{
  cellIndex_[0] = halfEdgeIndex_ / numDofsPerCell_;
  localEdgeIndex_[0] = ( halfEdgeIndex_ % numDofsPerCell_ ) / numDofsPerInteriorEdge_;
#ifndef NDEBUG
  if (localEdgeIndex_[0]<0 || localEdgeIndex_[0] >=6)
    {
      std::cerr << "error"<< std::endl;
      exit(1);
    }
#endif
};


template <VolumeType::enum_t _volumeType> inline void get_edge_to_nodes(const int_t 	cnc_[],
									const int_t 	localEdgeIndex_,
									int_t 		cncedge_[])		    
{
  cncedge_[0] = cnc_[ReferenceShapeVolume<_volumeType>::EdgesToNodes[localEdgeIndex_][0]];
  cncedge_[1] = cnc_[ReferenceShapeVolume<_volumeType>::EdgesToNodes[localEdgeIndex_][1]];
};



template <VolumeType::enum_t _volumeType,FaceType::enum_t _faceType>
struct Extract
{
public: static inline void GetFaceToNodes(const int_t 	cnc_[],
					  unsigned int 	localFaceIndex_,
					  int_t 	cncface_[]);
};

template <VolumeType::enum_t _volumeType> struct Extract<_volumeType,FaceType::Triangle>
{
public: static inline void GetFaceToNodes(const int_t 	cnc_[],
					  unsigned int 	localFaceIndex_,
					  int_t 	cncface_[])
  {
    cncface_[0] = cnc_[ReferenceShapeVolume<_volumeType>::TrianglesToNodes[localFaceIndex_][0]];
    cncface_[1] = cnc_[ReferenceShapeVolume<_volumeType>::TrianglesToNodes[localFaceIndex_][1]];
    cncface_[2] = cnc_[ReferenceShapeVolume<_volumeType>::TrianglesToNodes[localFaceIndex_][2]];
  };
};

template <VolumeType::enum_t _volumeType> struct Extract<_volumeType,FaceType::Quadrilateral>
{
public: static inline void GetFaceToNodes(const int_t 	cnc_[],
					  unsigned int	localFaceIndex_,
					  int_t 	cncface_[]) noexcept
  {
    cncface_[0] = cnc_[ReferenceShapeVolume<_volumeType>::QuadrilateralsToNodes[localFaceIndex_][0]];
    cncface_[1] = cnc_[ReferenceShapeVolume<_volumeType>::QuadrilateralsToNodes[localFaceIndex_][1]];
    cncface_[2] = cnc_[ReferenceShapeVolume<_volumeType>::QuadrilateralsToNodes[localFaceIndex_][2]];
    cncface_[3] = cnc_[ReferenceShapeVolume<_volumeType>::QuadrilateralsToNodes[localFaceIndex_][3]];
  };
};

//
// Definition of a MeshHalfFacet
//
template <typename _derived> struct Traits_CRTP_MeshHalfFacet;

template <typename _derived> struct CRTP_MeshHalfFacet
{
private: inline const _derived & asImp() const { return static_cast<const _derived&>(*this); }    
private: inline _derived & asImp()  { return static_cast<_derived&>(*this); }    

private: using traits_t = Traits_CRTP_MeshHalfFacet<_derived>;
public: CRTP_MeshHalfFacet(const CRTP_MeshHalfFacet<_derived>&) = delete;  

protected: inline CRTP_MeshHalfFacet(){};
  
public: using cell_t 	= typename traits_t::cell_t;
  
public: inline cell_t		GetMeshCell() 		const noexcept
  {
    return asImp().GetMeshCell();
  };
  
public: inline unsigned int	GetLocalFacetIndex() 	const noexcept
  {
    return asImp().GetLocalFacetIndex();
  };  
};




#if 0
template <unsigned int _degree> struct FiniteElement::Lagrange<_degree,DimensionType::Volume,VolumeType::Pyramid>
{
private: static constexpr unsigned int s_numVerticesInCell 	 = ReferenceShapeVolumePyramid::NbNodes;
private: static constexpr unsigned int s_numEdgesInCell 	 = ReferenceShapeVolumePyramid::NbEdges;
private: static constexpr unsigned int s_numFacesInCell 	 = ReferenceShapeVolumePyramid::NbFaces;
  //
  // 6*n = (_degree+1)*(_degree+1)*(_degree+1) + ( 2 * (_degree-1) ) * 8
  //
  
public: static constexpr unsigned int NumDofsPerCell 	 	= ((_degree+1)*(_degree+1)*(_degree)) / 2;
public: static constexpr unsigned int NumDofsPerInteriorEdge 	= (_degree>0) ? _degree-1 : 0;
  
public: static constexpr unsigned int NumDofsPerInteriorFace[]
  { (_degree>2) ? ((_degree-1)*(_degree-2))/2 : 0,
      (_degree>0) ? (_degree-1)*(_degree-1) : 0 };
  
public: static constexpr unsigned int NumDofsPerInteriorVolume = NumDofsPerCell
    - s_numVerticesInCell
    - s_numEdgesInCell * NumDofsPerInteriorEdge
    - 4 * NumDofsPerInteriorFace[0]
    - 1 * NumDofsPerInteriorFace[1];
  
};
#endif


template <typename fe_t>
void GenerateVolume(const int_t numCells_,
		    const int_t*__restrict__ cellsToNodes_,
		    const int_t cellsToNodesLd_,
		    int_t* __restrict__ numDofs_,
		    int_t* __restrict__ cellsToDofs_,
		    const int_t cellsToDofsLD_)
{
  int_t dofIndex = numDofs_[0];
  static constexpr unsigned int s_numDofsPerInteriorVolume = fe_t::NumDofsPerEntities[DimensionType::Volume][VolumeType::Tetrahedron];  
  if (s_numDofsPerInteriorVolume>0)
    {
      for (int_t cellIndex=0;cellIndex<numCells_;++cellIndex)
	{
	  auto p = &cellsToDofs_[cellsToDofsLD_*cellIndex];
	  for (unsigned int k = 0;k<s_numDofsPerInteriorVolume;++k)
	    {
	      p[k] = ++dofIndex;
	    }
	}
    }
  
#ifndef NDEBUG
  std::cout << "volume table " << std::endl;    
  //
  // mark face
  //
  if (s_numDofsPerInteriorVolume>0)
    {
      for (int_t cellIndex=0;cellIndex<numCells_;++cellIndex)
	{
	  auto p = &cellsToDofs_[cellsToDofsLD_*cellIndex];
	  for (unsigned int k = 0;k<s_numDofsPerInteriorVolume;++k)
	    {
	      std::cout << " " << p[k];
	    }
	  std::cout << std::endl;
	}
    }
#endif
  numDofs_[0] = dofIndex;
};

#include "Hasher.hpp"



template <typename fe_t,FaceType::enum_t _faceType, typename int_t>
unsigned long long int GenerateFace(const int_t 		numCells_,
				    const int_t*__restrict__ 	cellsToNodes_,
				    const int_t 		cellsToNodesLd_,
				    int_t*			numDofs_,
				    int_t* __restrict__ 	cellsToDofs_,
				    const int_t 		cellsToDofsLd_,
				    Hasher<int_t>&		hasher_)
{
  int_t dofIndex = numDofs_[0];
  using CellShape = ReferenceShapeVolume<fe_t::VolumeKind>;
  static constexpr unsigned int _degree 		 = fe_t::Degree;
  static constexpr unsigned int s_numVerticesInCell 	 = CellShape::NbNodes;
  static constexpr unsigned int s_numEdgesInCell 	 = CellShape::NbEdges;
  static constexpr unsigned int s_numFacesInCell 	 = CellShape::NbFaces;

  static constexpr unsigned int s_numDofsPerCell 	 	= fe_t::NumDofs;
  static constexpr unsigned int s_numDofsPerInteriorEdge 	= fe_t::NumDofsPerEntities[DimensionType::Edge][EdgeType::Edge];
  static constexpr unsigned int s_numDofsPerInteriorFace 	= fe_t::NumDofsPerEntities[DimensionType::Face][_faceType];
  
  std::cout << "s_numDofsPerCell           " << s_numDofsPerCell           << std::endl;
  std::cout << "s_numDofsPerInteriorFace   " << s_numDofsPerInteriorFace   << std::endl;

  unsigned long long int numFaces = 0;
  
  //
  // Set dofs over the faces.
  //
  if (s_numDofsPerInteriorFace > 0)
    {
      static constexpr const  int_t s = (_degree-1)*(_degree-1);

      int_t* permloc = new int_t[s > s_numDofsPerInteriorFace ? s : s_numDofsPerInteriorFace];
      {      
	int_t faceToNodes[4];	  
	int_t cellToNodes[s_numVerticesInCell];
	//	static constexpr const unsigned int s_localShift = s_numVerticesInCell + s_numEdgesInCell * s_numDofsPerInteriorEdge;
	for (int_t cellIndex=0;cellIndex<numCells_;++cellIndex)
	  {
	    //
	    // Copy connectivity.
	    //
	    get_cell_to_nodes(s_numVerticesInCell,cellsToNodes_,
						cellsToNodesLd_,					    
						cellIndex,
						cellToNodes);
	    
	    int_t at = cellsToDofsLd_ * cellIndex; // + s_localShift;
	    for (unsigned int localFaceIndex=0;localFaceIndex<s_numFacesInCell;++localFaceIndex)
	      {	      
		//
		// Extract the face.
		//
		Extract<fe_t::VolumeKind,_faceType>::GetFaceToNodes(cellToNodes,
							       localFaceIndex,
							       faceToNodes);
		
		//
		// Compute hash value.
		//
		const int_t hashValue = FaceCompare<_faceType>::HashCode(numCells_,
									 faceToNodes);
		
		//
		// Assign the last half edge index with the same hash value.
		//
		//		cellsToDofs_[at] = link[hashValue];
		cellsToDofs_[at] = hasher_[hashValue];
		//
		// Ujpdate the last half edge index with respect to the hash value.
		//
		hasher_[hashValue] = -at;
		at += s_numDofsPerInteriorFace;
	      }	
	  }
      }
      
      int_t faceToNodes[4];	  
      int_t cellToNodes[s_numVerticesInCell];
      int_t interiorFaceDofIndices[s_numDofsPerInteriorFace];    
      for (int_t cellIndex=numCells_-1;cellIndex>=0;--cellIndex)
	{	
	  
	  get_cell_to_nodes(s_numVerticesInCell,cellsToNodes_,
					      cellsToNodesLd_,					    
					      cellIndex,
					      cellToNodes);
	  
	  for (int localFaceIndex = s_numFacesInCell-1;localFaceIndex>=0;--localFaceIndex)
	    {
	      //s_numVerticesInCell + s_numEdgesInCell * s_numDofsPerInteriorEdge
	      const int_t initialHalfFaceIndex = -(cellsToDofsLd_*cellIndex + s_numDofsPerInteriorFace*localFaceIndex + 0);
	      int_t nextHalfFaceIndex = cellsToDofs_[-initialHalfFaceIndex];
	      if (nextHalfFaceIndex < 0)
		{
		  ++numFaces;
		  Extract<fe_t::VolumeKind,_faceType>::GetFaceToNodes(cellToNodes,
								 localFaceIndex,
								 faceToNodes);
		  
		  //
		  // Assign the interior dofs to the current new face
		  //
		  FaceCompare<_faceType>::template OrientationPermutation<_degree>((char)1,
										   permloc);
		  
		  int_t oldDofIndex = dofIndex;
		  
		  {
		    const int_t at = cellsToDofsLd_ * cellIndex
		      // + s_numVerticesInCell
		      // + s_numEdgesInCell * s_numDofsPerInteriorEdge
		      + s_numDofsPerInteriorFace * localFaceIndex;
		    for (unsigned int i=0;i<s_numDofsPerInteriorFace;++i)
		      {
			cellsToDofs_[at+i] = dofIndex + 1 + permloc[i];
		      }
		    dofIndex += s_numDofsPerInteriorFace;
		  }
		  
		  int_t lastDifferentHalfFaceIndex = initialHalfFaceIndex;
		  while (nextHalfFaceIndex < 0 && nextHalfFaceIndex != Hasher<int_t>::s_default_hash)
		    {		    
		      int_t testedCellIndex;
		      int_t testedLocalFaceIndex;
		      // - s_numVerticesInCell - s_numEdgesInCell * s_numDofsPerInteriorEdge
		      HalfFaceDecomposition(-nextHalfFaceIndex,
					    s_numDofsPerCell,
					    s_numDofsPerInteriorFace,
					    &testedCellIndex,
					    &testedLocalFaceIndex);
		      
		      //
		      // We need to extract the next half face before we overwrite it.
		      //
		      const int_t nextHalfFaceIndexBackup = nextHalfFaceIndex;
		      nextHalfFaceIndex = cellsToDofs_[-nextHalfFaceIndex];
		      
		      int_t testedCellToNodes[s_numVerticesInCell];
		      get_cell_to_nodes(s_numVerticesInCell,cellsToNodes_,
							  cellsToNodesLd_,					    
							  testedCellIndex,
							  testedCellToNodes);
		    
		      int_t testedFaceToNodes[4];
		      Extract<fe_t::VolumeKind,_faceType>::GetFaceToNodes(testedCellToNodes,
								     testedLocalFaceIndex,
								     testedFaceToNodes);
		      
		      if (FaceCompare<_faceType>::AreSame(faceToNodes,
							  testedFaceToNodes))
			{
			  FaceCompare<_faceType>::template OrientationPermutation<_degree>
			    (FaceCompare<_faceType>::Orientation(faceToNodes,
								 testedFaceToNodes),
			     permloc);

			  //
			  // We need to rebuild the linked list.
			  //
			  if (lastDifferentHalfFaceIndex != initialHalfFaceIndex)
			    {				
			      cellsToDofs_[-lastDifferentHalfFaceIndex] = nextHalfFaceIndex;
			    }

			  // s_numVerticesInCell+ s_numEdgesInCell * s_numDofsPerInteriorEdge
			  const int_t at = cellsToDofsLd_*testedCellIndex + s_numDofsPerInteriorFace*testedLocalFaceIndex;
			  
			  // we copy 
			  for (unsigned int i=0;i<s_numDofsPerInteriorFace;++i)
			    {
			      cellsToDofs_[at + i] =  oldDofIndex + 1 + permloc[i];
			    }
			  
			  break;
			
			}
		      else
			{
			  //
			  // this is NOT the same face.
			  //
			  lastDifferentHalfFaceIndex = nextHalfFaceIndexBackup;
			}
		    
		    }		
		}
	    }
	}

#ifndef NDEBUG
      std::cout << "face table " << std::endl;
  
      //
      // mark face
      //
      for (int_t cellIndex=0;cellIndex<numCells_;++cellIndex)
	{
	  for (int_t localFaceIndex=0;localFaceIndex<s_numFacesInCell;++localFaceIndex)
	    {
	      for (unsigned int k = 0;k<s_numDofsPerInteriorFace;++k)
		{
		  //
		  // only the first dof on edge.
		  //
		  // + s_numVerticesInCell + s_numEdgesInCell * s_numDofsPerInteriorEdge + 
		  std::cout << " " << cellsToDofs_[cellsToDofsLd_ * cellIndex + s_numDofsPerInteriorFace * localFaceIndex+k];
		}
	    }      
	  std::cout << std::endl;
	}

#endif
      std::cout << "numFaces " << numFaces << std::endl;
      std::cout << "dofIndex " << dofIndex << std::endl;
    }
  
  numDofs_[0] = dofIndex;
  return numFaces;  
};


template <typename fe_t, typename int_t>
void GenerateEdgeHashBasedLinkedList(const int_t numCells_,
				     const int_t*__restrict__ cellsToNodes_,
				     const int_t cellsToNodesLd_,
				     int_t* __restrict__ cellsToDofs_,
				     const int_t cellsToDofsLd_,
				     Hasher<int_t>&hasher_)
{
  using Cell = ReferenceShapeVolume<fe_t::VolumeKind>;
  static constexpr unsigned int s_numVerticesInCell 	 	= Cell::NbNodes;
  static constexpr unsigned int s_numEdgesInCell 	 	= Cell::NbEdges;
  static constexpr unsigned int s_numDofsPerCell 	 	= fe_t::NumDofs;
  static constexpr unsigned int s_numDofsPerInteriorEdge 	= fe_t::NumDofsPerEntities[DimensionType::Edge][EdgeType::Edge];
  
  std::cout << "s_numDofsPerCell           " << s_numDofsPerCell           << std::endl;
  std::cout << "s_numDofsPerInteriorEdge   " << s_numDofsPerInteriorEdge   << std::endl;
  if (s_numDofsPerInteriorEdge>0)
    {
      int_t edgeToNodes[2];	  
      int_t cellToNodes[s_numVerticesInCell];
      for (int_t cellIndex=0;cellIndex<numCells_;++cellIndex)
	{
	  //
	  // Copy connectivity.
	  //
	  get_cell_to_nodes(s_numVerticesInCell,cellsToNodes_,
					      cellsToNodesLd_,
					      cellIndex,
					      cellToNodes);
	    
	  int_t at = cellsToDofsLd_ * cellIndex;// + s_numVerticesInCell;
	  for (unsigned int localEdgeIndex=0;localEdgeIndex<s_numEdgesInCell;++localEdgeIndex)
	    {	      
	      //
	      // Extract the edgeToNodes.
	      //
	      get_edge_to_nodes<fe_t::VolumeKind>(cellToNodes,
					  localEdgeIndex,
					  edgeToNodes);
		
	      //
	      // Compute hash value.
	      //
	      const int_t hashValue = EdgeCompare::HashCode(numCells_,
							    edgeToNodes);
		
	      //
	      // Assign the last half edge index with the same hash value.
	      //
	      cellsToDofs_[at] = hasher_[hashValue];
		
	      //
	      // Update the last half edge index with respect to the hash value.
	      //
	      hasher_[hashValue] = -at;
	      at += s_numDofsPerInteriorEdge;
	    }	
	}
    }
  
};



template <typename fe_t, typename int_t>
unsigned long long int GenerateEdge(const int_t numCells_,
				    const int_t*__restrict__ cellsToNodes_,
				    const int_t cellsToNodesLd_,
				    int_t*numDofs_,
				    int_t* __restrict__ cellsToDofs_,
				    const int_t cellsToDofsLd_,
				    Hasher<int_t>&hasher_)
{
  int_t dofIndex = numDofs_[0];

  using CellShape = ReferenceShapeVolume<fe_t::VolumeKind>;
  static constexpr unsigned int s_numVerticesInCell 	 = CellShape::NbNodes;
  static constexpr unsigned int s_numEdgesInCell 	 = ReferenceShapeVolumeTetrahedron::NbEdges;
  

  static constexpr unsigned int s_numDofsPerCell 	 	= fe_t::NumDofs;
  static constexpr unsigned int s_numDofsPerInteriorEdge 	= fe_t::NumDofsPerEntities[DimensionType::Edge][EdgeType::Edge];
  
  std::cout << "s_numDofsPerCell           " << s_numDofsPerCell           << std::endl;
  std::cout << "s_numDofsPerInteriorEdge   " << s_numDofsPerInteriorEdge   << std::endl;
  
  unsigned long long int numEdges = 0;
  if (s_numDofsPerInteriorEdge>0)
    {      
      static constexpr const int_t default_hash = Hasher<int_t>::s_default_hash;
      GenerateEdgeHashBasedLinkedList<fe_t,int_t>(numCells_,
						  cellsToNodes_,
						  cellsToNodesLd_,
						  cellsToDofs_,
						  cellsToDofsLd_,
						  hasher_);
      int_t edgeToNodes[2];	  
      int_t cellToNodes[s_numVerticesInCell];
      int_t interiorEdgeDofIndices[s_numDofsPerInteriorEdge];    
      for (int_t cellIndex=numCells_-1;cellIndex>=0;--cellIndex)
	{	
	  get_cell_to_nodes(s_numVerticesInCell,cellsToNodes_,
					      cellsToNodesLd_,
					      cellIndex,
					      cellToNodes);

	  //
	  // reverse edge.
	  //
	  for (int localEdgeIndex = s_numEdgesInCell-1;localEdgeIndex>=0;--localEdgeIndex)
	    {
	      //s_numVerticesInCell
	      const int_t initialHalfEdgeIndex = -(cellsToDofsLd_*cellIndex + s_numDofsPerInteriorEdge*localEdgeIndex + 0);
	      int_t nextHalfEdgeIndex = cellsToDofs_[-initialHalfEdgeIndex];
	      if (nextHalfEdgeIndex < 0)
		{
		  ++numEdges;
		  get_edge_to_nodes<fe_t::VolumeKind>(cellToNodes,
					      localEdgeIndex,
					      edgeToNodes);
		  
		  const bool edgeHasPositiveOrientation = EdgeCompare::HasPositiveOrientation(edgeToNodes);
		  //
		  // that's the first visit of the edge
		  //
		  if (edgeHasPositiveOrientation)
		    {
		      for (unsigned int i=0;i<s_numDofsPerInteriorEdge;++i)
			{
			  interiorEdgeDofIndices[i] = ++dofIndex;
			}		    
		    }
		  else
		    {
		      for (unsigned int i=0;i<s_numDofsPerInteriorEdge;++i)
			{
			  interiorEdgeDofIndices[s_numDofsPerInteriorEdge - 1 - i] = ++dofIndex;
			}
		    }
		
		  //
		  // Assign the interior dofs to the current new edge
		  //
		  {
		    // s_numVerticesInCell
		    const int_t at = cellsToDofsLd_ * cellIndex  + s_numDofsPerInteriorEdge * localEdgeIndex;
		    for (unsigned int i=0;i<s_numDofsPerInteriorEdge;++i)
		      {
			cellsToDofs_[at+i] = interiorEdgeDofIndices[i];
		      }
		  }
		
		  int_t lastDifferentHalfEdgeIndex = initialHalfEdgeIndex;
		  while (nextHalfEdgeIndex < 0 && nextHalfEdgeIndex != default_hash)
		    {		    
		      int_t testedCellIndex;
		      int_t testedLocalEdgeIndex;
		      HalfEdgeDecomposition(-nextHalfEdgeIndex,// - s_numVerticesInCell,
					    s_numDofsPerCell,
					    s_numDofsPerInteriorEdge,
					    &testedCellIndex,
					    &testedLocalEdgeIndex);
		    
		      //
		      // We need to extract the next half edge before we overwrite it.
		      //
		      const int_t nextHalfEdgeIndexBackup = nextHalfEdgeIndex;
		      nextHalfEdgeIndex = cellsToDofs_[-nextHalfEdgeIndex];
		    
		      int_t testedCellToNodes[s_numVerticesInCell];
		      get_cell_to_nodes(s_numVerticesInCell,cellsToNodes_,
							  cellsToNodesLd_,
							  testedCellIndex,
							  testedCellToNodes);
		    
		      int_t testedEdgeToNodes[2];
		      get_edge_to_nodes<fe_t::VolumeKind>(testedCellToNodes,
						  testedLocalEdgeIndex,
						  testedEdgeToNodes);
		    
		      if (EdgeCompare::AreSame(edgeToNodes,
					       testedEdgeToNodes))
			{
			  //
			  // We need to rebuild the linked list.
			  //
			  if (lastDifferentHalfEdgeIndex != initialHalfEdgeIndex)
			    {				
			      cellsToDofs_[-lastDifferentHalfEdgeIndex] = nextHalfEdgeIndex;
			    }
			  
			  // s_numVerticesInCell
			  const int_t at = s_numDofsPerCell*testedCellIndex+s_numDofsPerInteriorEdge*testedLocalEdgeIndex;
			  
			  if (EdgeCompare::HasPositiveOrientation(testedEdgeToNodes)
			      == edgeHasPositiveOrientation)
			    {
			      // copy 
			      for (unsigned int i=0;i<s_numDofsPerInteriorEdge;++i)
				{
				  cellsToDofs_[at + i] = interiorEdgeDofIndices[i];
				}
			    }
			  else
			    {
			      // copy reverse
			      for (unsigned int i=0;i<s_numDofsPerInteriorEdge;++i)
				{
				  cellsToDofs_[at + i] = interiorEdgeDofIndices[s_numDofsPerInteriorEdge - 1 - i];
				}
			    }
			}
		      else
			{
			  //
			  // this is NOT the same edge.
			  //
			  lastDifferentHalfEdgeIndex = nextHalfEdgeIndexBackup;
			}
		    }		
		}
	    }
	}
    
#ifndef NDEBUG      
      std::cout << "edge table " << std::endl;
      
      //
      // mark edge
      //
      for (int_t cellIndex=0;cellIndex<numCells_;++cellIndex)
	{
	  for (unsigned int localEdgeIndex=0;localEdgeIndex<s_numEdgesInCell;++localEdgeIndex)
	    {
	      for (unsigned int k = 0;k<s_numDofsPerInteriorEdge;++k)
		{
		  //
		  // only the first dof on edge.
		  //
		  // s_numVerticesInCell + 
		  std::cout << " " << cellsToDofs_[cellsToDofsLd_ * cellIndex + s_numDofsPerInteriorEdge * localEdgeIndex+k];
		}
	    }      
	  std::cout << std::endl;
	}    
#endif
    }
  
  numDofs_[0] = dofIndex;
  return numEdges;
};








//
//
// mixed shapes
//
// order of the shapes:
//
// tets, wedges, pyramids, hexahedrons
//
// Do edges,
//
// [-max -n0 -n1 ... -n2 ... -n3] [-max ... -n4 ... -n5  ... -n6  ... -n7]
//
// treat edges hexahedron
// assign dofs
// foreach edge in inplace linked lists
//
//     if (linked > lowerbound)
//     {
//        then it belongs to hexahedron
//        we continue to assign dofs
//     }
//    else
//     {
//        it points to the first edge with same hash code in another shape
//        we marked this edge as a shared edge
//        
//     }
//
// Do tets
//  Do triangles and leave the boundary ones, we will treat them after.
//
// Do wedges
//  Do triangles and leave the boundary ones, we will treat them after, they can be shared by other types.
//  Do quadrilaterals and leave the boundary ones, we will treat them after, they can be shared by other types.
//
// Do pyramids
//  Do triangles and leave the boundary ones, we will treat them after, they can be shared by other types.
//  Do quadrilaterals and leave the boundary ones, we will treat them after, they can be shared by other types.
//
// Do hexahedron
//  Do triangles and leave the boundary ones, we will treat them after.
//  Do quadrilaterals and leave the boundary ones, we will treat them after.
//
// Do remaining triangles
//
// Do remaining quadrilaterals
// 
//

template <unsigned int _degree,
	  VolumeType::enum_t _volumeType,
	  typename int_t>
int_t FiniteElementSpaceSize()
{
  return FiniteElement::Lagrange<_degree, DimensionType::Volume, _volumeType>::NumDofs;
}


template <unsigned int 		_degree,
	  VolumeType::enum_t 	_volumeType,
	  typename int_t>
void  GenerateFiniteElementSpace(const int_t 			numCells_,
				 const int_t*__restrict__ 	cellsToNodes_,
				 const int_t 			cellsToNodesLd_,
				 int_t*				numDofs_,
				 int_t*__restrict__ 		cellsToDofs_,
				 const int_t 			cellsToDofsLd_)
{
  //
  // The cell type.
  //
  using CellShape = ReferenceShapeVolume<_volumeType>;

  static constexpr unsigned int s_numNodesInCell = CellShape::NbNodes;
  static constexpr unsigned int s_numEdgesInCell = CellShape::NbEdges;
  static constexpr unsigned int s_numFacesInCell = CellShape::NbFaces;
  static constexpr unsigned int s_numFaceTInCell = CellShape::NbTriangleFaces;
  static constexpr unsigned int s_numFaceQInCell = CellShape::NbQuadrilateralFaces;
  
  //
  // The finite element type.
  //
  using fe_t = FiniteElement::Lagrange<_degree, DimensionType::Volume, _volumeType>;

  //
  // Get the number of dofs on each support.
  //
  static constexpr unsigned int s_numDofsPerCell
    = fe_t::NumDofs;
  static constexpr unsigned int s_numDofsPerInteriorEdge
    = fe_t::NumDofsPerEntities[DimensionType::Edge][EdgeType::Edge];
  static constexpr unsigned int s_numDofsPerInteriorFaceTriangle
    = fe_t::NumDofsPerEntities[DimensionType::Face][FaceType::Triangle];
  static constexpr unsigned int s_numDofsPerInteriorFaceQuadrilateral
    = fe_t::NumDofsPerEntities[DimensionType::Face][FaceType::Quadrilateral];
  static constexpr unsigned int s_numDofsPerInteriorVolume
    = fe_t::NumDofsPerEntities[DimensionType::Volume][_volumeType];    

  //
  // Print information.
  //
  std::cout << "s_numDofsPerCell                        "
	    << s_numDofsPerCell
	    << std::endl;
  
  std::cout << "s_numDofsPerInteriorEdge                "
	    << s_numDofsPerInteriorEdge
	    << std::endl;
  
  std::cout << "s_numDofsPerInteriorFaceTriangle        "
	    << s_numDofsPerInteriorFaceTriangle
	    << std::endl;
  
  std::cout << "s_numDofsPerInteriorFaceQuadrilateral   "
	    << s_numDofsPerInteriorFaceQuadrilateral
	    << std::endl;
  
  std::cout << "s_numDofsPerInteriorVolume              "
	    << s_numDofsPerInteriorVolume
	    << std::endl;
  
  //
  // Set dofs over the vertices
  //
  int_t dofIndex = numDofs_[0];
  
  //
  // Copy
  //
  int_t shift = 0;
  bool hasP1 = true;
  if ( hasP1 )
    {
      std::cout << "COPY P1 space ..." << std::endl;
      //
      // Copy.
      //
      for (int_t cellIndex=0;cellIndex<numCells_;++cellIndex)
	{	  
	  auto cellToNodes = &cellsToNodes_[cellsToNodesLd_*cellIndex];
	  auto cellToDofs = &cellsToDofs_[cellsToDofsLd_*cellIndex];
	  for (unsigned int j=0;j<s_numNodesInCell;++j)
	    {
	      cellToDofs[j] = cellToNodes[j];
	    }  
	}

      //
      // Max
      //
      for (int_t cellIndex=0;cellIndex<numCells_;++cellIndex)
	{
	  auto cellToNodes = &cellsToNodes_[cellsToNodesLd_*cellIndex];
	  
	  for (unsigned int j=0;j<s_numNodesInCell;++j)
	    {
	      auto nodeIndex = cellToNodes[j];
	      if (nodeIndex > dofIndex)
		{
		  dofIndex = nodeIndex;
		}
	    }
	  
	}
      std::cout << "COPY P1 space done." << std::endl;
      shift += 4;
    }
  
  unsigned long long int numEdges = 0;
  unsigned long long int numFaces = 0;
  
  //
  // Edges.
  //
  if (s_numDofsPerInteriorEdge>0)
    {
      //
      // Set dofs over the edges.
      //
      Hasher<int_t> hasher(numCells_);      
      numEdges = GenerateEdge<fe_t>(numCells_,
				    cellsToNodes_,
				    cellsToNodesLd_,
				    &dofIndex,
				    cellsToDofs_ + shift,
				    cellsToDofsLd_,
				    hasher);
      shift += s_numDofsPerInteriorEdge * s_numEdgesInCell;
    }
  
  std::cout << "numEdges " << numEdges << std::endl;
  std::cout << "dofIndex " << dofIndex << std::endl;
  
  //
  // Generate faces.
  //
  {    
    Hasher<int_t> hasher(numCells_);
    if (s_numDofsPerInteriorFaceTriangle>0)
      {
	std::cout << "generate face " << dofIndex << std::endl;
	numFaces =  GenerateFace<fe_t, FaceType::Triangle, int_t>(numCells_,
										  cellsToNodes_,
										  cellsToNodesLd_,
										  &dofIndex,
										  cellsToDofs_ + shift,
										  cellsToDofsLd_,
										  hasher);
	std::cout << "generate face done " << dofIndex << std::endl;	
	shift += s_numDofsPerInteriorFaceTriangle * s_numFaceTInCell;
      }
    
    if (s_numDofsPerInteriorFaceQuadrilateral>0)
      {
	hasher.Reset();
	numFaces +=  GenerateFace<fe_t, FaceType::Quadrilateral, int_t>(numCells_,
										     cellsToNodes_,
										     cellsToNodesLd_,
										     &dofIndex,
										     cellsToDofs_ + shift,
										     cellsToDofsLd_,
										     hasher);
	shift += s_numDofsPerInteriorFaceTriangle * s_numFaceQInCell;	
      }
    
  }

  std::cout << "numFaces " << numFaces << std::endl;  
  std::cout << "dofIndex " << dofIndex << std::endl;
  if (s_numDofsPerInteriorVolume>0)
    {
      GenerateVolume<fe_t>(numCells_,
			   cellsToNodes_,
			   cellsToNodesLd_,
			   &dofIndex,
			   cellsToDofs_+shift,
			   cellsToDofsLd_);
    }

  std::cout << "NUMDOFS " << dofIndex << std::endl;
  
#if 0
  
  if (s_numDofsPerInteriorVolume>0)
    {
      //
      // Set dofs over the volumes.
      //
      for (int_t cellIndex=0;cellIndex<numCells_;++cellIndex)
	{
	  const int_t at = s_numDofsPerCell * cellIndex + s_numNodesInCell + s_numEdgesInCell * s_numDofsPerInteriorEdge + s_numDofsPerInteriorFace * s_numFacesInCell;
	  for (unsigned int k = 0;k<s_numDofsPerInteriorVolume;++k)
	    {
	      cellsToDofs[at + k] = ++dofIndex;
	    }     
	}
    }
  
#ifndef NDEBUG
  std::cout << "volume table " << std::endl;    
  //
  // mark face
  //
  if (s_numDofsPerInteriorVolume>0)
    {
      for (int_t cellIndex=0;cellIndex<numCells_;++cellIndex)
	{
	  const int_t at = s_numDofsPerCell * cellIndex + s_numNodesInCell + s_numEdgesInCell * s_numDofsPerInteriorEdge + s_numDofsPerInteriorFace * s_numFacesInCell;
	  for (unsigned int k = 0;k<s_numDofsPerInteriorVolume;++k)
	    {
	      //
	      // only the first dof on edge.
	      //
	      std::cout << " " << cellsToDofs[at +k];
	    }
	  std::cout << std::endl;
	}
    }
#endif
#endif  
  numDofs_[0] = ++dofIndex;
  //  return cellsToDofs;  
};


