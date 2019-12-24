#include <valarray>

#include <iostream>
#include <ostream>
#include "Program.hpp"
#include <limits>
#include <array>

#include "DimensionType.hpp"
#include "ReferenceShape.hpp"


#include "ReferenceShapeNode.hpp"
#include "ReferenceShapeEdge.hpp"
#include "ReferenceShapeFaceTriangle.hpp"
#include "ReferenceShapeFaceQuadrilateral.hpp"
#include "ReferenceShapeVolumeTetrahedron.hpp"
#include "ReferenceShapeVolumePyramid.hpp"
#include "ReferenceShapeVolumeWedge.hpp"
#include "ReferenceShapeVolumeHexahedron.hpp"

#include "Config.hpp"
#include "Point.hpp"
#include "CellToNodes.hpp"


#include "Input/Medit.hpp"
//#include "Mesh.hpp"
#include "Output/Medit.hpp"
#include "Output/Vtk.hpp"

#include "FiniteElement/Lagrange.hpp"
#include "FiniteElement/LagrangeNode.hpp"
#include "FiniteElement/LagrangeEdge.hpp"
#include "FiniteElement/LagrangeFaceTriangle.hpp"
#include "FiniteElement/LagrangeFaceQuadrilateral.hpp"
#include "FiniteElement/LagrangeVolumeHexahedron.hpp"
#include "FiniteElement/LagrangeVolumeTetrahedron.hpp"
#include "FiniteElement/LagrangeVolumeWedge.hpp"


#include "EdgeCompare.hpp"
#include "FaceCompare.hpp"
#include "AHF/Mesh.hpp"

#include "generate_finite_element_space.hpp"
#include "Treilli.hpp"
#include "CellsToCellsCalculator.hpp"



template <int i>
static int comp(const void * a,const void * b)
{
  const int_t * a_ = (const int_t * )a;
  const int_t * b_ = (const int_t * )b;
  if (a_[i] < b_[i]) return -1;
  else
    if (a_[i] > b_[i]) return 1;
  else
    return 0;
}

unsigned long long int hilbert(double crd[3], double box[6], int itr)
{
  unsigned long long int IntCrd[3], m=1LL<<63, cod;
  int i, j, b, GeoWrd, NewWrd, BitTab[3] = {1,2,4};
  double TmpCrd[3];
  int rot[8], GeoCod[8]={0,3,7,4,1,2,6,5};  /* Z curve = {5,4,7,6,1,0,3,2} */
  int HilCod[8][8] = {{0,7,6,1,2,5,4,3}, {0,3,4,7,6,5,2,1}, {0,3,4,7,6,5,2,1}, {2,3,0,1,6,7,4,5}, \
		      {2,3,0,1,6,7,4,5}, {6,5,2,1,0,3,4,7}, {6,5,2,1,0,3,4,7}, {4,3,2,5,6,1,0,7}};
  
  /* Convert double precision coordinates to integers */
  
  
  for(j=0;j<3;j++)
    {
      TmpCrd[j] = (crd[j] - box[j]) * box[j+3];
      double a = TmpCrd[j];
      unsigned long long int * pp = (unsigned long long int*)&a;
      IntCrd[j] = pp[0];
#if 0
      printf("[%d]=%Lu\n",j,IntCrd[j]);
#endif
    }

  /* Binary hilbert renumbering loop */
  
  cod = 0;
  
  for(j=0;j<8;j++)
    rot[j] = GeoCod[j];
  
  for(b=0;b<itr;b++)
    {
      GeoWrd = 0;
      
      for(j=0;j<3;j++)
	{
	  if(IntCrd[j] & m)
	    GeoWrd |= BitTab[j];
	  
	  IntCrd[j] = IntCrd[j]<<1;
	}
      
      NewWrd = rot[ GeoWrd ];
      
      cod = cod<<3 | NewWrd;

      for(j=0;j<8;j++)
	rot[j] = HilCod[ NewWrd ][ rot[j] ];
      
    }
#if 0
  printf("cod %Lu\n",cod);  
#endif
  return(cod);
}

static void coordinates3d(int_t 				N_, 
			  double*				box,
			  const double * vertex_,
			  int_t vertexoff_,
			  unsigned long long int * cod2_,
			  int_t 				codoff_)
{
  double loc[3];
  for(int_t i=0; i<N_; i++)
    {
      loc[0] = vertex_[i*vertexoff_+0];
      loc[1] = vertex_[i*vertexoff_+1];
      loc[2] = vertex_[i*vertexoff_+2];
      cod2_[codoff_*i+0] 	= hilbert(loc,box,21);
      cod2_[codoff_*i+1] 	= i;
    }
}

#if 0

static void coordinates2d(int_t				N_, 
			  double				box_[],
			  double 				vertex_[],
			  int_t                         vertexoff_,
			  unsigned long long int *  	cod2_,
			  int_t				codoff_)
{
  static const unsigned long long m = 1LL<<62;
  static const int BitTab[2] 	= {1,2};
  static const int GeoCod[4]	= {1,2,0,3};
  static const int HilCod[4][4] = {{0,3,2,1}, {0,1,2,3}, {0,1,2,3}, {2,1,0,3}};
  static const double len 	= 4.611686018427387904e+18;
  double box[4],dbl;
  unsigned long long int IntCrd[2],cod;
  int rot[4];
  box[0] 	= box_[0];
  box[1] 	= box_[1];
  box[2] 	= len / (box_[2] - box_[0]);
  box[3] 	= len / (box_[3] - box_[1]);
  double loc[2];
  { int i;
    for(i=0; i<N_; i++)
      {
	loc[0] = vertex_[vertexoff_*i+0];
	loc[1] = vertex_[vertexoff_*i+1];
	/* Convert double precision coordinates to integers */
	dbl 	= (loc[0] - box[0]) * box[0+2];
	IntCrd[0] = (unsigned long long int)dbl;
	dbl 	= (loc[1] - box[1]) * box[1+2];
	IntCrd[1] = (unsigned long long int)dbl;
	/* Binary hilbert renumbering loop */
	cod = 0;


	rot[0] = GeoCod[0];
	rot[1] = GeoCod[1];
	rot[2] = GeoCod[2];
	rot[3] = GeoCod[3];
	{ int b;
	  for(b=0;b<31;b++)
	    {
	      int GeoWrd = 0;

	      if(IntCrd[0] & m)
		GeoWrd |= BitTab[0];
	      IntCrd[0] = IntCrd[0]<<1;

	      if(IntCrd[1] & m)
		GeoWrd |= BitTab[1];
	      IntCrd[1] = IntCrd[1]<<1;

	      const int NewWrd = rot[ GeoWrd ];
	      cod = cod<<2 | NewWrd;
	      rot[0] = HilCod[ NewWrd ][ rot[0] ];
	      rot[1] = HilCod[ NewWrd ][ rot[1] ];
	      rot[2] = HilCod[ NewWrd ][ rot[2] ];
	      rot[3] = HilCod[ NewWrd ][ rot[3] ];
	    } }
	cod2_[codoff_*i+0] 	= cod;
	cod2_[codoff_*i+1] 	= i;
      } }

}
#endif

class MnsFiniteElementSpace3D : public Program
{
  
private: int_t m_degree;
public:
  
  MnsFiniteElementSpace3D(int 	 argc,
			  char * argv[]) : Program(argc, argv, true)
  {
    this->AddOption(new Option<int_t>(&m_degree,
				      "-k",
				      2,
				      true,
				      "Lagrange"));
  };

  template <typename _derivedClass>  void extract (const CRTP_MeshTopology<_derivedClass> &meshTopology_,
						   int_t * cnc,
						   int_t cncLd)
  {
    using this_t = CRTP_MeshTopology<_derivedClass>;
    using entitykind_t = typename this_t::entitykind_t;
    static constexpr DimensionType::enum_t meshDimension = this_t::Dimension;
    for( const auto cellType : entitykind_t::All)
      {
	const auto numCellsOfType = meshTopology_.template GetNumEntities<meshDimension>(cellType);
	if ( numCellsOfType > 0 )
	  {
	    typename this_t::template entity_t<DimensionType::Node> cellToNodes[8];
	    for (const auto cell : meshTopology_.template GetEntities<meshDimension>(cellType) )
	      {
		const auto numNodesInCell = meshTopology_.template GetEntityToEntities<meshDimension,DimensionType::Node>(cell,
															  cellToNodes);
		const auto cellIndex = meshTopology_.template GetEntityIndex<DimensionType::Volume>(cell);		
		for (unsigned int i=0;i<numNodesInCell;++i)
		  {
		    const auto nodeIndex = meshTopology_.template GetEntityIndex<DimensionType::Node>(cellToNodes[i]);
		    cnc[cellIndex * cncLd + i] = nodeIndex;
		  }
	      }
	  }	
      }
  };


  template <typename _derivedClass>  void extract (const CRTP_MeshGeometry<_derivedClass> &meshGeometry_,
						   double * coo,
						   int_t cooLd)
  {
    auto meshTopology_ = meshGeometry_.GetMeshTopology();
    int_t numPoints = meshGeometry_.GetNumPoints();
    for (const auto node : meshTopology_->template GetEntities<DimensionType::Node>(NodeType::Node) )
      {
	int_t index = meshTopology_->template GetEntityIndex<DimensionType::Node>(node);		
	auto p = meshGeometry_.GetPoint(node);
	coo[cooLd*index+0] = p[0];
	coo[cooLd*index+1] = p[1];
	coo[cooLd*index+2] = p[2];
      }
  };

  



  
  virtual void Main()
  {
    

    if (this->HasVerbose())
      {
	LogMessage  << "degree "
		    << m_degree
		    << std::endl;
	
	LogMessage  << "ofilename "
		    << this->GetOfilename()
		    << std::endl;
      }
    
    if (Mpi::IsMaster())
      {
	{
	  std::cout << "master pre-processing"  << std::endl;
	  AHF::Mesh3D mesh(this->GetInputFilename(0).c_str());	  
#if 0
	  Output::Medit outputMedit("spaceAHF");
	  outputMedit << mesh;
	  
	  Output::Vtk::Writer outputVtk("spaceVtk");
	  outputVtk << mesh;
#endif	  
	  const auto topology 		= mesh.GetTopology();
	  const auto geometry 		= mesh.GetGeometry();

	  const auto numPyramids 	= topology->GetNumEntities<DimensionType::Volume>(VolumeType::Pyramid);
	  const auto numWedges 		= topology->GetNumEntities<DimensionType::Volume>(VolumeType::Wedge);
	  const auto numTetrahedra 	= topology->GetNumEntities<DimensionType::Volume>(VolumeType::Tetrahedron);
	  const auto numHexahedra 	= topology->GetNumEntities<DimensionType::Volume>(VolumeType::Hexahedron);
	  const auto numNodes 		= topology->GetNumEntities<DimensionType::Node>();
	  const auto numCells 		= numHexahedra +  numWedges + numTetrahedra + numPyramids;
#if 0	  
	  std::cout << "num nodes " << numNodes << std::endl;
	  std::cout << "num cells " << numCells << std::endl;
	  std::cout << " -  " << VolumeType::Tetrahedron << " " << numTetrahedra << std::endl;
	  std::cout << " -  " << VolumeType::Pyramid << " " << numPyramids << std::endl;
	  std::cout << " -  " << VolumeType::Wedge << " " << numWedges << std::endl;
	  std::cout << " -  " << VolumeType::Hexahedron << " " << numHexahedra << std::endl;
	  
#endif
	  if (numTetrahedra > 0)
	    {
	      static constexpr int_t cncLd = 5;
	      static constexpr int_t c2cLd = 10;
	      int_t * cnc = new int_t[numCells * 5];
	      int_t * c2c = new int_t[numCells * 10];
	      std::cout << "extract  ..." << std::endl;
	      extract(*topology,cnc,cncLd);
	      std::cout << "cells to cells ..." << std::endl;

	      double * mesh_coo = new double[3*numNodes];
	      extract(*geometry,mesh_coo,3);

	      //
	      // Compute hilbert coordinates of the nodes.
	      //
	      if (0)
	      {
		std::cout << "hilbert coordinates"  << std::endl;
		double * mesh_coo2 = new double[5*numNodes];
		extract(*geometry,mesh_coo2,5);
		double box_[]={-1.2,-1.2,-1.2,1.2,1.2,1.2};

		coordinates3d(numNodes,
			      box_,
			      mesh_coo2,
			      5,
			      (unsigned long long int*)mesh_coo2+3,
			      5);
		std::cout << "hilbert coordinates done"  << std::endl;
		qsort(mesh_coo2,numNodes,5*sizeof(int_t),comp<3>);
		auto a = (unsigned long long int*)mesh_coo2+3;
		int_t * perm = new int_t[numNodes];
		for (int_t i=0;i<numNodes;++i)
		  {
#if 0
		    std::cout << mesh_coo2[5*i+0] << " " <<
		      mesh_coo2[5*i+1] << " " << 
		      mesh_coo2[5*i+2] << " " <<
		      a[5*i+0] << " " <<
		      a[5*i+1] << " " << std::endl;
#endif
		    perm[a[5*i+1]] = i;
		  }
		
		for (int_t i=0;i<numNodes;++i)
		  {
		    mesh_coo[3*i+0] = mesh_coo2[3*i+0];
		    mesh_coo[3*i+1] = mesh_coo2[3*i+1];
		    mesh_coo[3*i+2] = mesh_coo2[3*i+2];
		  }
		
		for (int_t i=0;i<numCells;++i)
		  {
		    for(int j=0;j<4;++j)
		      {
			auto k = cnc[cncLd*i+j];
			cnc[cncLd*i+j] = perm[k];
		      }
		  }
		
		delete [] perm;
		perm = nullptr;

		

		


		
		// int_t p = i;
		// int_t q = a[5*i+0];		
	      }
	      
	      
#if 1
	      //
	      // Calculate the cells-to-cells.
	      //
	      
	      unsigned long long int outNumFaces;
	      unsigned long long int outNumInteriorFaces;
	      unsigned long long int outNumBoundaryFaces;
	      
	      calculator_cells_to_cells_t<VolumeType::Tetrahedron,FaceType::Triangle,int_t>
		::calculate(numCells,cnc,cncLd,c2c,c2cLd,&outNumFaces,&outNumInteriorFaces,&outNumBoundaryFaces);
	      
	      std::cout << "numFaces " << outNumFaces << std::endl;
	      std::cout << "numInteriorFaces " << outNumInteriorFaces << std::endl;
	      std::cout << "numBoundaryFaces " << outNumBoundaryFaces << std::endl;

	      //
	      // Transform into a sparse graph for METIS.
	      //

	      //
	      // METIS.
	      //

	      
	      exit(1);
	      
#endif
	      

	      
	      const int_t numDofsPerCell = FiniteElementSpaceSize<3,VolumeType::Tetrahedron,int_t>();	  
	      int_t * cncP = new int_t[numDofsPerCell * numCells];	  
	      int_t numDofs = 0;



	      
	      GenerateFiniteElementSpace<3,VolumeType::Tetrahedron,int_t>(numCells,
									  cnc,
									  cncLd,
									  &numDofs,
									  cncP,
									  numDofsPerCell);

	      Uniform3d<VolumeType::Tetrahedron,3> uniform;
	      


	      auto lnnodes = uniform.nnodes();
	      std::cout << "MeshVersionFormatted"
			<< std::endl
			<< "1"
			<< std::endl
			<< "Dimension"
			<< std::endl
			<< "3"
			<< std::endl
			<< "Vertices"
			<< std::endl
			<< numDofs
			<< std::endl;

	      double * coo = new double[3*numDofs];
	      extract(*geometry,coo,3);
	      bool *marker = new bool[numDofs];
	      for (int_t i=0;i<numDofs;++i)marker[i]=false;

	      double cooel[512];
	      for (int_t i=0;i<numCells;++i)
		{
		  //
		  // Get the coordinates of the cell.
		  //
		  
		  //
		  // Interpolate the coordinates from uniform.
		  //
		  for (unsigned int j=0;j<4;++j)
		    {
		      cooel[j*3+0] = coo[3*cncP[numDofsPerCell * i + j] + 0];
		      cooel[j*3+1] = coo[3*cncP[numDofsPerCell * i + j] + 1];
		      cooel[j*3+2] = coo[3*cncP[numDofsPerCell * i + j] + 2];		      
		    }

		  double ll[8]; // 8 being max nodes
		  for (unsigned int j=0;j<numDofsPerCell;++j)
		    {
		      auto dof = cncP[numDofsPerCell * i + j];
		      if (false == marker[dof])
			{
			  auto r = uniform.GetCoordinate<double>(j,0);
			  auto s = uniform.GetCoordinate<double>(j,1);
			  auto t = uniform.GetCoordinate<double>(j,2);
			  
			  ll[0] = 1.0 - r - s - t;
			  ll[1] = r;
			  ll[2] = s;
			  ll[3] = t;

			  double res = 0.0;
			  for (unsigned int k=0;k<4;++k)
			    {
			      res += ll[k] * cooel[k*3+0];
			    }
			  
			  coo[3*dof+0] = res;
			  res = 0.0;
			  for (unsigned int k=0;k<4;++k)
			    {
			      res += ll[k] * cooel[k*3+1];
			    }
			  coo[3*dof+1] = res;
			  res = 0.0;
			  for (unsigned int k=0;k<4;++k)
			    {
			      res += ll[k] * cooel[k*3+2];
			    }
			  coo[3*dof+2] = res;
			  
			  marker[cncP[numDofsPerCell * i + j]] = true;
			}
		    }
		}
	      
	      for (int_t i=0;i<numDofs;++i)
		{
		  std::cout << coo[3*i+0] << " " << coo[3*i+1] << " " << coo[3*i+2] << " 0" <<std::endl;
		}
	      
	      auto lnsubcells = uniform.nsubcells();
	      auto lnnodesincell = uniform.nnodesincell();
	  
	      std::cout << ((lnnodesincell==4)?"Tetrahedra" : "Hexahedra")
			<< std::endl
			<< lnsubcells * numCells
			<< std::endl;
	      for (int_t i=0;i<numCells;++i)
		{
		  for (unsigned int isub=0;isub<lnsubcells;++isub)
		    {
		      for (unsigned int j=0;j<lnnodesincell;++j)
			{
			  std::cout << " " << cncP[numDofsPerCell * i + uniform.GetNodeIndex(isub,j)]+1;
			}
		      std::cout << " 0"
				<< std::endl;
		    }
		  std::cout << std::endl;
		}
	      std::cout << "End"  << std::endl;
	      
	      //
	      // Generate coordinate on edges.
	      //

	      //
	      // Generate coordinate on faces.
	      //

	      //
	      // Generate coordinate on volumes.
	      //


	      
	    }
#if 1
#if 1
	  
	  //
	  // 
	  //
	  
	  std::cout << "FINITE ELEMENT SPACE #####" << std::endl;
#endif
#endif

	  


	  
//	  const int_t numCells_,
//									     const int_t*__restrict__ cellsToNodes_,
//									     int_t*numDofs_)
//	  
	  //
	  // Generate the finite element space
	  //
	  

	}



	
#if 0
	{	  
	  Input::Medit inputMedit(this->GetInputFilename(0).c_str());
	  
	  AHF::MeshTopology3D ahfMeshTopology3D(inputMedit);
	  
	  Output::Medit outputMedit("spaceAHFTopology");
	  outputMedit << ahfMeshTopology3D;
	  
	}
#endif


#if 0
	Mesh<3,double,DimensionType::Volume> mesh(this->GetInputFilename(0).c_str());
	{
	  Output::Medit outputMedit("space");
	  outputMedit << *mesh.m_topology;//VolumeType::Tetrahedron;
	}
	
	{
	  MeshEntity<DimensionType::Volume> meshEntity(VolumeType::Hexahedron,4);
	  
	  std::cout << "encoding " << meshEntity.m_encoding << std::endl;
	  std::cout << "type " << meshEntity.GetType() << std::endl;
	  std::cout << "index " << meshEntity.GetIndex() << std::endl;
	}
	
	{
	  MeshEntity<DimensionType::Face> meshEntity(FaceType::Quadrilateral, 4);
	  
	  std::cout << "encoding " << meshEntity.m_encoding << std::endl;
	  std::cout << "type " << meshEntity.GetType() << std::endl;
	  std::cout << "index " << meshEntity.GetIndex() << std::endl;
	}

	{
	  MeshEntity<DimensionType::Face> meshEntity(FaceType::Quadrilateral, 5);
	  
	  std::cout << "encoding " << meshEntity.m_encoding << std::endl;
	  std::cout << "type " << meshEntity.GetType() << std::endl;
	  std::cout << "index " << meshEntity.GetIndex() << std::endl;
	}
#endif
	
	
#if 0
	static constexpr unsigned int _degree = 29;
	int_t numDofs;
	int_t* cellsToDofs = GenerateFiniteElementSpace<_degree,VolumeType::Hexahedron,int_t>(10,
											       NULL,
											       &numDofs);

#endif	
#if 0
	
	// const unsigned int numIfilenames = this->GetNumInputFiles();
	pMedit medit = Medit_new(this->GetInputFilename(0).c_str());

	eVolume convert[] = {__eVolume_TETRAHEDRON,
			     __eVolume_PYRAMID,
			     __eVolume_WEDGE,
			     __eVolume_HEXAHEDRON};
	
	
	I numCells[4]{0};
	for(const auto volumeKind : VolumeType::AllValues)
	  {
	    numCells[volumeKind] = Medit_get_nbVolumes	(medit,
							 convert[volumeKind]);
	    std::cout << "# "  << volumeKind << " " << numCells[volumeKind] << std::endl;
	  }
	I totalNumCells = numCells[0] + numCells[1] + numCells[2] + numCells[3];
	I numDifferentCells = (numCells[0] == 0) ? 0 : 1;
	numDifferentCells += (numCells[1] == 0) ? 0 : 1;
	numDifferentCells += (numCells[2] == 0) ? 0 : 1;
	numDifferentCells += (numCells[3] == 0) ? 0 : 1;
	pI begin = new I[4+1];
	begin[0] = 0;
	
	const I nn[]= {4,5,6,8};
	{
	  unsigned int i=1;
	  for(const auto volumeKind : VolumeType::AllValues)
	    {
	      begin[i] = begin[i-1] + numCells[volumeKind] * nn[volumeKind];
	      ++i;
	    }
	}
	
	
 	pI tetcod = new I[totalNumCells];
	pI cellsToNodes = new I[4*numCells[0] + 5*numCells[1] + 6*numCells[2] + 8*numCells[3]];
	for(const auto volumeKind : VolumeType::AllValues)
	  {
	    if (numCells[volumeKind]>0)
	      {		
		Medit_get_cncVolume(medit,
				    &cellsToNodes[begin[volumeKind]],
				    nn[volumeKind],
				    //
				    //
				    // Come and correct tetcod, 
				    //
				    //
				    tetcod,
				    1,
				    convert[volumeKind]);
	      }
	    
	  }
#if 0
 	pI tetcod = new I[numTetrahedrons];
	pI cellsToNodes = new I[4*numTetrahedrons];
	
	
	
	I numTetrahedrons = Medit_get_nbVolumes	(medit,
						 __eVolume_TETRAHEDRON);
#endif
	I numVertices = Medit_get_nbVertices	(medit);
	
	std::cout << "1##############################" << std::endl;
	pPoints xyz = Points_new(__eDim_3,
				 numVertices);
	
	pI cod = new I[numVertices];
	I n1 = 1;
	Medit_get_Points(medit,
			 xyz,
			 cod,
			 &n1);

	I numTetrahedrons = numCells[0];
#if 0
	std::cout << "1##############################" << std::endl;
 	pI tetcod = new I[numTetrahedrons];
	pI cellsToNodes = new I[4*numTetrahedrons];
	
	Medit_get_cncVolume(medit,
			    cellsToNodes,
			    4,
			    tetcod,
			    1,
			    __eVolume_TETRAHEDRON);
	
#endif
#if 0	
	std::cout << "numtet " << numTetrahedrons << std::endl;
	for (I i=0;i<numTetrahedrons;++i)
	  {
	    std::cout << " " << cellsToNodes[4*i+0]
		      << " " << cellsToNodes[4*i+1]
		      << " " << cellsToNodes[4*i+2]
		      << " " << cellsToNodes[4*i+3]
		      << std::endl;
	  }
#endif
	
	static constexpr unsigned int _degree = 7;
	I numDofs;
	pI cellsToDofs = GenerateFiniteElementSpace<_degree,VolumeType::Tetrahedron,I>(numTetrahedrons,
										       cellsToNodes,
										       &numDofs);
	
	





	
	std::cout << "NumDofs " << numDofs << std::endl;
	double*coo = new double[numDofs * 3];
	pTreilliVolume treilli = TreilliVolume_new(__eVolume_TETRAHEDRON,
						   _degree);
	I nbv = TreilliVolume_get_nbVertices(treilli);
	
	double*rst = new double[nbv*3];
	TreilliVolume_get_coo_double(treilli,
				     rst,
				     nbv);
	double x0[3];
	double x1[3];
	double x2[3];
	double x3[3];
	bool *flag=new bool[numDofs];
	for (I i=0;i<numDofs;++i)
	  {
	    flag[i]=false;
	  }
	for (I i = 0 ;i<numTetrahedrons;++i)
	  {
	    Points_get		(xyz,
				 cellsToNodes[4*i+0],
				 x0);
	    Points_get		(xyz,
				 cellsToNodes[4*i+1],
				 x1);
	    Points_get		(xyz,
				 cellsToNodes[4*i+2],
				 x2);
	    Points_get		(xyz,
				 cellsToNodes[4*i+3],
				 x3);

	    for (I j=0;j<nbv;++j)
	      {
		I k = cellsToDofs[nbv*i+j];
		
		const double
		  r = rst[j],
		  s = rst[nbv+j],
		  t = rst[2*nbv+j];

		const double
		  l0 = 1.0 -r -s-t,
		  l1 = r,
		  l2= s,
		  l3 = t;

		const double
		  x = l0 * x0[0] + l1*x1[0] + l2*x2[0] + l3*x3[0],
		  y = l0 * x0[1] + l1*x1[1] + l2*x2[1] + l3*x3[1],
		  z = l0 * x0[2] + l1*x1[2] + l2*x2[2] + l3*x3[2];
		
		if (k>=numDofs)
		  {
		    std::cerr << "ddd" << std::endl;
		    exit(1);
		  }
		flag[k]=true;
		coo[3*k+0] = x;
		coo[3*k+1] = y;
		coo[3*k+2] = z;
	      }
	  }

	bool invalid = false;
	for (I i=0;i<numDofs;++i)
	  {
	    if (flag[i]==false)
	      {
		invalid = true;
		std::cout << "invalid " << i << std::endl;
	      }
	  }

	if (invalid)
	  {
	    exit(1);
	  }
	FILE * out = fopen("out.mesh","w");
	fprintf(out,"MeshVersionFormatted\n1\nDimension\n3\nVertices\n" ifmt "\n",numDofs);
#if 1
	for (I i =0;i<numDofs;++i)
	  {
	    fprintf(out,"%e %e %e 0\n",coo[3*i+0],coo[3*i+1],coo[3*i+2]);
	  }
#endif	
	I nbSubVolumes = TreilliVolume_get_nbSubvolumes	(treilli);
	fprintf(out,"Tetrahedra\n" ifmt "\n",numTetrahedrons * nbSubVolumes);
	I* subcnc=new I [4*nbSubVolumes];
	TreilliVolume_get_cnc		(treilli,
					 subcnc,
					 4);
	//	printf("allo "ifmt"\n",nbv);
	//	exit(1);
	for (I i = 0 ;i<numTetrahedrons;++i)
	  {	    
	    for (I j = 0 ;j<nbSubVolumes;++j)
	      {
		for (I k= 0;k<4;++k)
		  {
		    //	    std::cout << "yo " << 4*j+k<<std::endl;
		    //		    std::cout << " " << subcnc[4*j+k];
		    fprintf(out," " ifmt "",cellsToDofs[nbv*i+subcnc[4*j+k]]+1);
		  }
		fprintf(out," 0\n");
		//	std::cout << std::endl;
	      }
	    //	    exit(1);
	  }

	fprintf(out,"End");
	fclose(out);
#endif
      }

  };
};


#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <assert.h>





template <typename int_t> struct MeshEntityHalfFace
{
private: int_t m_value;

public: using encoder_t = GenericEncoding<int_t,6>;
public: inline constexpr MeshEntityHalfFace(const int_t volumeId_,
					    const int_t localFaceIndex_) noexcept
  : m_value(encoder_t::Encod(volumeId_,localFaceIndex_))
  {    
  };
  
public: inline constexpr MeshEntityHalfFace(const int_t value_) noexcept
  : m_value(value_)
  {
  };
  
public: inline constexpr int_t GetId() noexcept
  {
    return m_value;
  };

public: inline constexpr int_t GetLow() noexcept
  {
    return encoder_t::Low(m_value);
  };

public: inline constexpr int_t GetUp() noexcept
  {
    return encoder_t::Up(m_value);
  };
  
};





using HalfFace = MeshEntityHalfFace<long long int>;
int main(int 	argc_, 
	 char* argv_[])
{
  
#if 0
  int h = 10;
  HalfFace::encoder_t::Info();
  HalfFace e(10,1);
  std::cout << e.GetId() << std::endl;
  std::cout << e.GetLow() << std::endl;
  std::cout << e.GetUp() << std::endl;
  std::cout << sizeof(HalfFace) << std::endl;
  std::cout << sizeof(int) << std::endl;
  std::cout << sizeof(long int) << std::endl;
  std::cout << sizeof(long long int) << std::endl;
#endif  
  MnsFiniteElementSpace3D application(argc_,argv_);
  application.Run();
  return 0;
}


#if 0

  
  std::cout << ReferenceCell::Tetrahedron::GetNumEntities<DimensionType::Node>() << std::endl;
  std::cout << ReferenceCell::Tetrahedron::GetNumEntities<DimensionType::Edge>() << std::endl;
  std::cout << ReferenceCell::Tetrahedron::GetNumEntities<DimensionType::Face>() << std::endl;
  
  for(const auto& edgeToNodes : ReferenceCell::Tetrahedron::EdgesToNodes)
    {
      for (const auto localNodeIndex : edgeToNodes)
	{	  
	  std::cout << " " << localNodeIndex;	  
	}
      std::cout << std::endl;
    }
  
  std::cout << ReferenceCell::Wedge::GetNumEntities<DimensionType::Node>() << std::endl;
  std::cout << ReferenceCell::Wedge::GetNumEntities<DimensionType::Edge>() << std::endl;
  std::cout << ReferenceCell::Wedge::GetNumEntities<DimensionType::Face>() << std::endl;

#endif






