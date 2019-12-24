#pragma once

#include "Hasher.hpp"

#include "FiniteElement/Lagrange.hpp"
#include "FiniteElement/LagrangeNode.hpp"
#include "FiniteElement/LagrangeEdge.hpp"
#include "FiniteElement/LagrangeFaceTriangle.hpp"
#include "FiniteElement/LagrangeFaceQuadrilateral.hpp"
#include "FiniteElement/LagrangeVolumeHexahedron.hpp"
#include "FiniteElement/LagrangeVolumeTetrahedron.hpp"
#include "FiniteElement/LagrangeVolumeWedge.hpp"



template <VolumeType::enum_t _volumeType,FaceType::enum_t _faceType,typename _int_t>
class calculator_cells_to_cells_t
{  

  using CellShape = ReferenceShapeVolume<_volumeType>;

  static constexpr unsigned int s_numVerticesInCell 	 = CellShape::NbNodes;
  static constexpr unsigned int s_numEdgesInCell 	 = CellShape::NbEdges;
  static constexpr unsigned int s_numFacesInCell 	 = CellShape::NbFaces;

private: static inline void MyHalfFaceDecomposition(const int_t 	halfFaceIndex_,
						    const int_t 	ld_,
						    int_t* 		cellIndex_,
						    int_t* 		localFaceIndex_)
  {
    cellIndex_[0] = halfFaceIndex_ / ld_;
    localFaceIndex_[0] = ( halfFaceIndex_ % ld_ );
  };

  template<unsigned int n_>
  static inline void GetCellToNodes(const _int_t*__restrict__ 	cellsToNodes_,
			  const _int_t 			cellsToNodesLd_,				
			  _int_t 			cellIndex,
			  _int_t * 			cellToNodes)
  {
    for (unsigned int i =0;i<n_;++i)
      {
	cellToNodes[i] = cellsToNodes_[cellsToNodesLd_*cellIndex+i];
      }
  }

private: static inline void InitHasher(const _int_t 			numCells_,
				const _int_t*__restrict__ 	cellsToNodes_,
				const _int_t 			cellsToNodesLd_,
				_int_t*__restrict__ 		cellsToCells_,
				const _int_t 			cellsToCellsLd_,
				Hasher<_int_t>& hasher)
  {

    _int_t faceToNodes[4];	  
    _int_t cellToNodes[8];
    for (_int_t cellIndex = 0;cellIndex < numCells_;++cellIndex)
      {
	GetCellToNodes<s_numVerticesInCell>(cellsToNodes_,
					    cellsToNodesLd_,					    
					    cellIndex,
					    cellToNodes);
	for (unsigned int localFaceIndex = 0;localFaceIndex < s_numFacesInCell;++localFaceIndex)
	  { 

	    Extract<_volumeType,_faceType>::GetFaceToNodes(cellToNodes,
							   localFaceIndex,
							   faceToNodes);

	    auto hash_value = FaceCompare<_faceType>::HashCode(numCells_,
							       faceToNodes);
	    
	    cellsToCells_[cellIndex * cellsToCellsLd_ + localFaceIndex] = hasher[hash_value];	    
	    hasher[hash_value] = -(cellIndex * cellsToCellsLd_ + localFaceIndex + 1);
	  }
      }
#if 0
    for (_int_t cellIndex = 0;cellIndex < numCells_;++cellIndex)
      {
	for (unsigned int localFaceIndex = 0;localFaceIndex < s_numFacesInCell;++localFaceIndex)
	  {
	    std::cout << " " << cellsToCells_[cellIndex * cellsToCellsLd_ + localFaceIndex];
	  }
	std::cout << std::endl;
      }
#endif
    
    
  };
  //
  // 
  //
public: static void calculate(const _int_t 			numCells_,
			      const _int_t*__restrict__ 	cellsToNodes_,
			      const _int_t 			cellsToNodesLd_,
			      _int_t*__restrict__ 		cellsToCells_,
			      const _int_t 			cellsToCellsLd_,
			      unsigned long long int *__restrict__ outNumFaces_,
			      unsigned long long int *__restrict__ outNumFacesInterior_,
			      unsigned long long int *__restrict__ outNumFacesBoundary_)
  {
    
Hasher<_int_t> hasher(numCells_);
 

    InitHasher(numCells_,
	       cellsToNodes_,
	       cellsToNodesLd_,
	       cellsToCells_,
	       cellsToCellsLd_,					  
	       hasher);

    //
    //
    //
    _int_t faceIndex = 0;    
    unsigned long long int numFaces = 0;
    unsigned long long int numInteriorFaces = 0;
    
    _int_t faceToNodes[4];	  
    _int_t cellToNodes[8];

    _int_t testedFaceToNodes[4];
    _int_t testedCellToNodes[8];
    
    //
    // 4 volumes, hence 2 bits.
    //

    //
    // Loop over the cells.
    // 

    for (_int_t cellIndex = numCells_ -1;cellIndex >= 0;--cellIndex)
      {
	//
	// Get the cell-to-nodes
	//
	GetCellToNodes<s_numVerticesInCell>(cellsToNodes_,
					    cellsToNodesLd_,					    
					    cellIndex,
					    cellToNodes);
	
	//
	// Loop over the faces.
	//
	for (int localFaceIndex = s_numFacesInCell - 1;localFaceIndex >= 0;--localFaceIndex)
	  { 
	    const _int_t initialHalfFaceIndex = -(cellsToCellsLd_ * cellIndex + localFaceIndex + 1);
	    _int_t nextHalfFaceIndex = cellsToCells_[-initialHalfFaceIndex-1];
	    if (nextHalfFaceIndex < 0)
	      {
		++numFaces;
		
		//
		// Extract the face-to-nodes.
		//
		Extract<_volumeType,_faceType>::GetFaceToNodes(cellToNodes,
							       localFaceIndex,
							       faceToNodes);

		//
		// Assign the face index.
		//
		cellsToCells_[cellsToCellsLd_ * cellIndex + localFaceIndex] = 0;
		
		//
		// Now traverse the candidates.
		//
		_int_t lastDifferentHalfFaceIndex = initialHalfFaceIndex;
		while (nextHalfFaceIndex < 0 && nextHalfFaceIndex != Hasher<_int_t>::s_default_hash)
		  {

		    //
		    // Get the cell and the local face index.
		    //
		    _int_t
		      testedCellIndex,
		      testedLocalFaceIndex;
		    MyHalfFaceDecomposition(-nextHalfFaceIndex-1,
					    cellsToCellsLd_,
					    &testedCellIndex,
					    &testedLocalFaceIndex);
		    //
		    // We need to extract the next half face before we overwrite it.
		    //
		    const auto nextHalfFaceIndexBackup = nextHalfFaceIndex;
		    nextHalfFaceIndex = cellsToCells_[-nextHalfFaceIndex-1];

		    //
		    // Get the cell-to-nodes of the cell it is testing.
		    //
		    GetCellToNodes<s_numVerticesInCell>(cellsToNodes_,
							cellsToNodesLd_,					    
							testedCellIndex,
							testedCellToNodes);

		    //
		    // Get the facet-to-nodes of the face it is testing.
		    //		  
		    Extract<_volumeType,_faceType>::GetFaceToNodes(testedCellToNodes,
								   testedLocalFaceIndex,
								   testedFaceToNodes);
		  

		    //
		    // Compare the face.
		    //
		    if (FaceCompare<_faceType>::AreSame(faceToNodes,
							testedFaceToNodes))
		      {

			//
			// SAME FACE: it might need to rebuild the linked list.
			//
			if (lastDifferentHalfFaceIndex != initialHalfFaceIndex)
			  {				
			    cellsToCells_[-lastDifferentHalfFaceIndex-1] = nextHalfFaceIndex;
			  }

			//
			// Assign and break the loop.
			//

			// encoder.Encod(cellIndex+1, _volumeType);
			
			cellsToCells_[cellsToCellsLd_ * testedCellIndex + testedLocalFaceIndex] =  cellIndex + 1;
			cellsToCells_[cellsToCellsLd_ * cellIndex + localFaceIndex] =  testedCellIndex + 1;
			++numInteriorFaces;
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

    outNumFaces_[0] = numFaces;
    outNumFacesInterior_[0] = numInteriorFaces;
    outNumFacesBoundary_[0] = numFaces - numInteriorFaces; 
  };

};
