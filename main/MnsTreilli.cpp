// #include "Quadrature/Edge/Legendre.hpp"

#include "Treilli.hpp"

int main()
{
#if 0
  {
    Uniform3d<VolumeType::Wedge,3> uniform;
    auto nnodes = uniform.nnodes();
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
	      << nnodes
	      << std::endl;
    
    for (unsigned int i=0;i<nnodes;++i)
      {
      
	std::cout << uniform.GetCoordinate<double>(i,0)
		  << " "
		  << uniform.GetCoordinate<double>(i,1)
		  << " "
		  << uniform.GetCoordinate<double>(i,2)
		  << " 0" << std::endl;
      
      }
    auto nsubcells = uniform.nsubcells();
    auto n = uniform.nnodesincell();
    std::cout << ((n==4)?"Tetrahedra" : "Hexahedra")
	      << std::endl
	      << nsubcells
	      << std::endl;
    for (unsigned int i=0;i<nsubcells;++i)
      {
	for (unsigned int j=0;j<n;++j)
	  {
	    std::cout << " " << uniform.GetNodeIndex(i,j)+1;
	  }
	std::cout << " 0"
		  << std::endl;
      }
    std::cout << "End"  << std::endl;
    return 0;
  }
#endif  

  //  FiniteElement::Uniform<FaceType::Quadrilateral,21> uniform;
  double v[] = {
    -0.96816023950762608983557620290367287004940480049192,
    -0.83603110732663579429942978806973487654410671812468,
    -0.61337143270059039730870203934147418478572060494056,
    -0.32425342340380892903853801464333660857195626073698,
    0.0,
    0.32425342340380892903853801464333660857195626073697,
    0.61337143270059039730870203934147418478572060494056,
    0.83603110732663579429942978806973487654410671812468,
    0.96816023950762608983557620290367287004940480049192 };
  FiniteElement::Generator<FaceType::Quadrilateral,10> uniform(v);

  auto nnodes = uniform.nnodes();
  std::cout << "MeshVersionFormatted"
	    << std::endl
	    << "1"
	    << std::endl
	    << "Dimension"
	    << std::endl
	    << "2"
	    << std::endl
	    << "Vertices"
	    << std::endl
	    << nnodes
	    << std::endl;

  for (unsigned int i=0;i<nnodes;++i)
    {
      std::cout << uniform.GetCoordinate<double>(i,0) << " " << uniform.GetCoordinate<double>(i,1) << " 0" << std::endl;
    }
  auto nsubcells = uniform.nsubcells();
  auto n = uniform.nnodesincell();
  std::cout << ((n==3)?"Triangles" : "Quadrilaterals")
	    << std::endl
	    << nsubcells
	    << std::endl;
  for (unsigned int i=0;i<nsubcells;++i)
    {
      for (unsigned int j=0;j<n;++j)
	{
	  std::cout << " " << uniform.GetNodeIndex(i,j)+1;
	}
      std::cout << " 0"
		<< std::endl;
    }
  std::cout << "End"  << std::endl;

  
  
  
  return 0;  
}
