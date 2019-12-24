#pragma once

#include "MeshTopology.hpp"

template <> class MeshTopology<DimensionType::Face>
{
  
public: MeshTopology(const char * filename_) noexcept
  {
    Input::Medit inputMedit(filename_);

    //
    // Initialize the number of cells.
    //
    for(const auto faceType : FaceType::All)
      {
	this->m_numCellsPerKind[faceType] = inputMedit.GetNumFaces(faceType);	
	std::cout << "### "  << faceType << " " << this->m_numCellsPerKind[faceType] << std::endl;
      }

    //
    // Calculate the total number of cells.
    //
    this->m_numCells = this->m_numCellsPerKind[0] + this->m_numCellsPerKind[1];

    //
    // Allocate Triangles
    //
    if (this->m_numCellsPerKind[FaceType::Triangle]>0)
      {
	this->m_trianglesToNodes = new TriangleToNodes[this->m_numCellsPerKind[FaceType::Triangle]];
	inputMedit.ReadTopology(this->m_trianglesToNodes);	
      }
    
    //
    // Allocate Quadrilateral
    //
    if (this->m_numCellsPerKind[FaceType::Quadrilateral]>0)
      {
	this->m_quadrilateralsToNodes = new QuadrilateralToNodes[this->m_numCellsPerKind[FaceType::Quadrilateral]];
	inputMedit.ReadTopology(this->m_quadrilateralsToNodes);	
      }    
  };

public: inline ~MeshTopology()
  {

    if(nullptr != this->m_trianglesToNodes)
      {
	delete[] this->m_trianglesToNodes;
	this->m_trianglesToNodes = nullptr;
      }

    if(nullptr != this->m_quadrilateralsToNodes)
      {
	delete[] this->m_quadrilateralsToNodes;
	this->m_quadrilateralsToNodes = nullptr;
      }
  };
  
private:
  
  int_t 		m_numCells;
  int_t   		m_numCellsPerKind[2]{0};
  int_t * 		m_cellsToNodes[2]{nullptr};
  int_t * 		m_cellCodes[2]{nullptr};
  
  TriangleToNodes*	m_trianglesToNodes 	= nullptr;
  QuadrilateralToNodes*	m_quadrilateralsToNodes = nullptr;
};


