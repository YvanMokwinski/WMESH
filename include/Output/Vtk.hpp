#pragma once

#include <string>
#include "ReferenceShape.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "CRTP_Mesh.hpp"
#include "CRTP_MeshGeometry.hpp"

namespace Output
{
  namespace Vtk
  {
    class Writer
    {
    protected:
      std::string m_filename;
      std::ofstream m_out;
    public : 
      Writer(const std::string& filename_)
	: m_filename(filename_+".vtk"),
	  m_out(m_filename.c_str()) 
      {
	this->m_out.precision(15);
	this->m_out.setf(std::ios::scientific);
      };
      
      virtual ~Writer()
      {
	this->m_out.close();
      };
      
      template <typename _type> Writer& operator<<(const _type &type_)
      {
	m_out << type_;
	return *this;
      };
      
    };

    typedef enum __eVtkElement{ ERROR=0,
				VERTEX=1,
				POLYVERTEX=2,
				LINE=3,
				POLYLINE=4,
				TRIANGLE=5,
				TRIANGLE_STRIP=6,
				POLYGON=7,
				PIXEL=8,
				QUAD=9,
				TETRA=10,
				VOXEL=11,
				HEXA=12,
				WEDGE=13,
				PYRAMID=14,
				QUADRATICEDGE=21,
				QUADRATICTRIANGLE=22,
				QUADRATICQUAD=23,
				QUADRATICTETRA=24,
				QUADRATICHEXA=25,
				ALL} eVtkElement;
  
    inline std::ostream& operator<< (std::ostream &out_,
				     const VolumeType::enum_t &volumeType_) noexcept
    {
      switch(volumeType_)
	{
	case VolumeType::Tetrahedron:
	  {
	    out_ << TETRA;
	    break;
	  }
	case VolumeType::Wedge:
	  {
	    out_ << WEDGE;
	    break;
	  }
	case VolumeType::Hexahedron:
	  {
	    out_ << HEXA;
	    break;
	  }
	  
	case VolumeType::Pyramid:
	  {
	    out_ << PYRAMID;
	    break;
	  }
	}
      return out_;
    };  
  
    std::ostream& operator<< (std::ostream &out_,
			      const FaceType::enum_t &faceType_)    
    {
      switch(faceType_)
	{
	case FaceType::Quadrilateral:
	  {
	    out_ << QUAD;
	    break;
	  }
	case FaceType::Triangle:
	  {
	    out_ << TRIANGLE;
	    break;
	  }
	}
      return out_;
    };  


    template <typename _derivedClass>  std::ostream& operator<< (std::ostream &out_,
								 const CRTP_Mesh<_derivedClass> &mesh_)
    {
      //      typedef typename CRTP_ReadOnlyMesh<_derivedClass>::PointType PointType;
      //      PointType p;      

    
      const auto topology 	= mesh_.GetTopology();	
      const auto geometry 	= mesh_.GetGeometry();
      const auto numPoints  	= geometry->GetNumPoints();
    
      out_ << "# vtk DataFile Version 3.0"
	   << std::endl
	   << "vtk output"
	   << std::endl
	   << "ASCII"
	   << std::endl
	   << "DATASET UNSTRUCTURED_GRID"
	   << std::endl
	   << "POINTS "
	   << numPoints
	   << " double"
	   << std::endl
	   << *geometry
	   << *topology;
      return out_;
    };

    template <long unsigned int _dimension,typename _real_t> std::ostream& operator<< (std::ostream &out_,
										  const Point<_dimension,_real_t>&point_)    
    {
      for (const auto x : point_)
	{
	  out_ << " " << x;
	}
      // out_ << " " << point_.m_cod;
      return out_;
    };  

    template <typename _derivedClass>  std::ostream& operator<< (std::ostream &out_,
								 const CRTP_MeshGeometry<_derivedClass> &meshGeometry_)
    {

      const auto meshtopology = meshGeometry_.GetMeshTopology();
      for (const auto node : meshtopology->template GetEntities<DimensionType::Node>(NodeType::Node) )
	{
	  out_ << meshGeometry_.GetPoint(node) << std::endl;	
	}

      return out_;
    };

  
    template <typename _derivedClass>  std::ostream& operator<< (std::ostream &out_,
								 const CRTP_MeshTopology<_derivedClass> &meshTopology_)
    {
      using this_t = CRTP_MeshTopology<_derivedClass>;
      using entitykind_t = typename this_t::entitykind_t;
      static constexpr DimensionType::enum_t meshDimension = this_t::Dimension;
      
      const auto numCells = meshTopology_.template GetNumEntities<meshDimension>();    
      unsigned long long int total_numCellsData = 0;

      for( const auto cellType : entitykind_t::All)
	{
	  const auto numCellsOfType = meshTopology_.template GetNumEntities<meshDimension>(cellType);
	  total_numCellsData += numCellsOfType * ( 1 + VolumeType::GetNumNodes(cellType) );
	  std::cout << 1 + VolumeType::GetNumNodes(cellType)  << std::endl;
	}
      
      out_ << "CELLS "
	   << numCells
	   << " "
	   << total_numCellsData
	   << std::endl;
    
      for( const auto cellType : entitykind_t::All)
	{
	  typename this_t::template entity_t<DimensionType::Node> cellToNodes[8];
	  for (const auto cell : meshTopology_.template GetEntities<meshDimension>(cellType) )
	    {
	      const auto numNodesInCell = meshTopology_.template GetEntityToEntities<meshDimension,DimensionType::Node>(cell,
															cellToNodes);
	      out_ << numNodesInCell;
	      for (unsigned int i=0;i<numNodesInCell;++i)
		{
		  out_ << " " << meshTopology_.template GetEntityIndex<DimensionType::Node>(cellToNodes[i]);
		}
	      out_ << std::endl;
	    }
	}
    
      out_ << "CELL_TYPES " << numCells << std::endl;
      for( const auto cellType : entitykind_t::All)
	{ 
	  const auto numCellsOfType = meshTopology_.template GetNumEntities<meshDimension>(cellType);
	  for (int_t i=0;i<numCellsOfType;++i)	
	    {
	      out_ << cellType << std::endl;
	    }
	}
      return out_;

    };
    
  };  
  
};
