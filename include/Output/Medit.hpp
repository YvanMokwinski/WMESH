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
  class Medit
  {
  protected:
    std::string m_filename;
    std::ofstream m_out;
  public : 
    Medit(const std::string& filename_)
      : m_filename(filename_+".mesh"),
	m_out(m_filename.c_str()) 
    {
      this->m_out.precision(15);
      this->m_out.setf(std::ios::scientific);
    };
      
    virtual ~Medit()
    {
      this->m_out.close();
    };
      
    template <typename _type> Medit& operator<<(const _type &type_)
    {
      m_out << type_;
      return *this;
    };
      
  };


  template <long unsigned int _dimension,typename _real_t> std::ostream& operator<< (std::ostream &out_,
										const Point<_dimension,_real_t>&point_)    
  {
    for (const auto x : point_)
      {
	out_ << " " << x;
      }
    out_ << " " << point_.m_cod;
    return out_;
  };  

  
  std::ostream& operator<< (std::ostream &out_,
			    const FaceType::enum_t &faceType_)    
  {
    switch(faceType_)
      {
      case FaceType::Quadrilateral:
	{
	  out_ << "Quadrilaterals";
	  break;
	}
      case FaceType::Triangle:
	{
	  out_ << "Triangles";
	  break;
	}
      }
    return out_;
  };  

  std::ostream& operator<< (std::ostream &out_,
			    const VolumeType::enum_t &volumeType_)    
  {	
    switch(volumeType_)
      {
      case VolumeType::Tetrahedron:
	{
	  out_ << "Tetrahedra";
	  break;
	}
      case VolumeType::Hexahedron:
	{
	  out_ << "Hexahedra";
	  break;
	}
      case VolumeType::Pyramid:
	{
	  out_ << "Pyramids";
	  break;
	}
      case VolumeType::Wedge:
	{
	  out_ << "Prisms";
	  break;
	}
      }
    return out_;
  };


  template <typename _derivedClass>  std::ostream& operator<< (std::ostream &out_,
							       const CRTP_Mesh<_derivedClass> &mesh_)
  {
    out_ << "MeshVersionFormatted"
	 << std::endl
	 << "1"
	 << std::endl
	 << *mesh_.GetGeometry()
	 << std::endl
	 << *mesh_.GetTopology()
	 << "End" << std::endl;
    return out_;

  };

  template <typename _derivedClass>  std::ostream& operator<< (std::ostream &out_,
							       const CRTP_MeshGeometry<_derivedClass> &meshGeometry_)
  {

    const auto meshtopology = meshGeometry_.GetMeshTopology();
    out_ << std::endl
	 << "Dimension"
	 << std::endl
	 << meshGeometry_.GetDimension()
	 << "Vertices"
	 << std::endl
	 << meshtopology->template GetNumEntities<DimensionType::Node>()
	 << std::endl;
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
    for( const auto cellType : entitykind_t::All)
      {
	const auto numCellsOfType = meshTopology_.template GetNumEntities<meshDimension>(cellType);
	if ( numCellsOfType > 0 )
	  {
	    out_ << cellType
		 << std::endl
		 << numCellsOfType
		 << std::endl;

	    typename this_t::template entity_t<DimensionType::Node> cellToNodes[8];
	    for (const auto cell : meshTopology_.template GetEntities<meshDimension>(cellType) )
	      {
		
		const auto numNodesInCell = meshTopology_.template GetEntityToEntities<meshDimension,DimensionType::Node>(cell,
															  cellToNodes);

		for (unsigned int i=0;i<numNodesInCell;++i)
		  {
		    auto nodeIndex = meshTopology_.template GetEntityIndex<DimensionType::Node>(cellToNodes[i]);
		    out_ << " " << nodeIndex + 1;
		  }
		out_ << " 0" << std::endl;
		//		out_ << cell << std::endl;
	      }
	    
	  }	
      }
    return out_;
  };




  
  

  
#if 0  
    
#if 0
  std::ostream& operator<< (std::ostream &out_,
			    const ReferenceCell::Volume::EnumType &volumeType_)    
    const unsigned int dimension = mesh_.GetGeometryDimension();
  out_ << "MeshVersionFormatted" << std::endl << "1" << std::endl << "Dimension" << std::endl << dimension << std::endl;      
  out_ << "Vertices" << std::endl << mesh_.GetNumVertices() << std::endl;      
    
  typedef typename CRTP_ReadOnlyMesh<_derivedClass>::PointType PointType;
  PointType p;
      
  for (unsigned int vertexIndex=0;vertexIndex<mesh_.GetNumVertices();++vertexIndex)
    {
      mesh_.GetVertex(vertexIndex,p);	 
      out_ << p.Get(0);
      for (unsigned int i=1;i<dimension;++i)
	{
	  out_ << " " << p.Get(i);
	}
      out_ << " 0" << std::endl;
    }
#endif

  //    using 
  //    typedef typename CRTP_ReadOnlyMesh<_derivedClass>::MeshTopology MeshTopology;
  //	unsigned int cnc[4];
  for (unsigned int topologyIndex=0;topologyIndex<mesh_.GetNumZones();++topologyIndex)
    {
      const MeshTopology * topology = mesh_.GetMeshTopology(topologyIndex);
#if 0
      template <typename _entityImpl> inline void GetEntity(const unsigned int&cellIndex_,CRTP_Entity<_entityImpl>&entity_) const
      {		
      };
#endif
	    
      const unsigned int numCells = topology->GetNumCells();
      out_ << topology->GetTypeShape() << std::endl;
      out_ << numCells << std::endl;

      typedef typename ReadOnlyMesh::MeshView MeshView;

      typedef typename MeshView::CellIterator CellIterator;
#if 1
      const MeshView view = mesh_.GetMeshView();
      CellIterator end_iterator = view.template end<0>();
      for (CellIterator iterator = view.template begin<0>();
	   iterator!=end_iterator;
	   ++iterator)
	{
		
	  //		cout << "allo"  << endl;
	}

#endif
	    
	    
#if 0
      int c = mesh_.template begin<int>();
      cout << c << std::endl;
      long int cc = mesh_.template begin<long int>();
      cout << cc << std::endl;

#endif
#if 0
      for (Meshentity_iterator it = topology->begin();it != topology->end();++it)
	{
		
	}
#endif

#if 0

      if (topology->GetTypeShape()==Mesh::QUADRILATERAL)
	{		
	  for (unsigned int cellIndex=0;cellIndex < numCells;++cellIndex)
	    {
	      topology->GetCellToNodes(cellIndex,cnc);
	      out_ << 1+cnc[0] << " " << 1+cnc[1] << " " << 1+cnc[3] << " " << 1+cnc[2] << " 0" << std::endl;
	    }
	}
      else
	{
	  for (unsigned int cellIndex=0;cellIndex < numCells;++cellIndex)
	    {
	      topology->GetCellToNodes(cellIndex,cnc);
	      out_ << 1+cnc[0] << " " << 1+cnc[1] << " " << 1+cnc[2] <<  " 0" << std::endl;
	    }
	}

#endif

#if 0
      const unsigned int numCells = topology->GetNumCells();
      out_ << topology->GetTypeShape() << std::endl;
      out_ << numCells << std::endl;
      if (topology->GetTypeShape()==Mesh::QUADRILATERAL)
	{
	  for (unsigned int cellIndex=0;cellIndex < numCells;++cellIndex)
	    {
	      topology->GetCellToNodes(cellIndex,cnc);
	      out_ << 1+cnc[0] << " " << 1+cnc[1] << " " << 1+cnc[3] << " " << 1+cnc[2] << " 0" << std::endl;
	    }
	}
      else
	{
	  for (unsigned int cellIndex=0;cellIndex < numCells;++cellIndex)
	    {
	      topology->GetCellToNodes(cellIndex,cnc);
	      out_ << 1+cnc[0] << " " << 1+cnc[1] << " " << 1+cnc[2] <<  " 0" << std::endl;
	    }
	}
#endif
    }

  out_ << "End" << std::endl;      
  return out_;
};
#endif
};
