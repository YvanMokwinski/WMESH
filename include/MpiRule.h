#ifndef __HEADER_MpiRule_h__
#define __HEADER_MpiRule_h__
#include "Mpi.h"




typedef enum __eDistributedType { LOCAL = 0,   
				  DISTRIBUTED } DistributedType;  

template <DistributedType _distributedType = DISTRIBUTED> struct MpiRule
{
  static long long int SplitLength(const long long int& length_);
  static long long int ComputeStartIndex(const long long int& length_);
};


template <> struct MpiRule<DISTRIBUTED>
{
  static long long int SplitLength(const long long int& length_)
  {
#ifdef SIMLAB_MPI
    int numCpu = Mpi::GetNumCpus();
    if (Mpi::IsMaster())
      {	
	return length_/numCpu + length_%numCpu;
      }
    else
      {
	return length_/numCpu;
      }
#else
    return length_;
#endif
  };

  static long long int ComputeStartIndex(const long long int& length_)
  {
#ifdef SIMLAB_MPI
    int numCpu = Mpi::GetNumCpus();
    if (Mpi::IsMaster())
      {	
	return 0;
      }
    else
      {
	return length_/numCpu + length_%numCpu + length_/numCpu * (Mpi::GetCpuId()-1);
      }
#else 
    return 0;
#endif
  };
};

template <> struct MpiRule<LOCAL>
{
  static long long int SplitLength(const long long int& length_)
  {
    return length_;
  };

  static long long int ComputeStartIndex(const long long int& length_)
  {
    return 0;
  };
};

#endif
