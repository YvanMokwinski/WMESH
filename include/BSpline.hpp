#pragma once




template <unsigned int _dim> class BSpline
{

public:  BSpline(const int_t		numPoints_);

public:  BSpline(const int_t		numPoints_,
		 cst_pR 		const 	ctrlpts_,
		 const int_t 			ctrlptsoff_);
  
public: int_t 		GetNumPoints() const noexcept;
public: int_t 		GetNumEdges() const noexcept;
public: unsigned int 	GetDimension() const noexcept;
public: double 		GetLength() const noexcept;

public: void 		GetCtrlPoint(const int_t ctrlPointIndex_,
				     double coordinates_[]) const noexcept;
  
public: void 		GetPoint(const int_t pointIndex_,
				 double coordinates_[]) const noexcept;

public: void 		GetTangent(const int_t pointIndex_,
				   double coordinates_[]) const noexcept;
  
public: void		Eval(const double	s_,
			     const int_t	edgeIndex_,
			     double 	position_[]) const noexcept;
  
public: void 		GetDofGeometry		(const int_t 		edgeIndex_,
						 pR 			dofGeometry_,
						 cst_pI		 	dofGeometryOff_) const noexcept;
  
  
private:  int_t			m_numPoints;
private:  double		m_length;  
private:  Point<_dim,double>*  	m_pts;
private:  Point<_dim,double>*  	m_ctrlpts;
private:  double* 		m_derivatives;
private:  double*		m_normalized_length;
};



template <unsigned int _dim>
double BSpline<_dim>::GetLength() const noexcept
{
  return this->m_length;  
};

template <unsigned int _dim>
int_t 	BSpline<_dim>::GetNumPoints() const noexcept
{
  return this->m_numPoints;
}

template <unsigned int _dim>
int_t 	BSpline<_dim>::GetNumEdges() const noexcept
{
  return this->m_numPoints-1;
}

template <unsigned int _dim>
unsigned int BSpline<_dim>::GetDimension() const noexcept
{
  return _dim;
}

template <unsigned int _dim>
void	BSpline<_dim>::GetCtrlPoint(const int_t 			ctrlPointIndex_,
				    pR 			const 	ctrlPointCoordinates_)
{
  for (unsigned int dimIndex = 0;dimIndex < _dim;++dimIndex)
    {
      ctrlPointCoordinates_[dimIndex] = this->m_ctrlpts[ctrlPointIndex_*_dim + dimIndex];
    }
}



template <unsigned int _dim>
void	BSpline<_dim>::GetTangent(const int_t 			pointIndex_,
				  pR 			const 	coordinates_)
{
  for (unsigned int dimIndex = 0;dimIndex < _dim;++dimIndex)
    {
      coordinates_[dimIndex] = this->m_derivatives[ctrlPointIndex_*_dim + dimIndex];
    }
}

template <unsigned int _dim>
void	BSpline<_dim>::GetPoint(const int_t 			pointIndex_,
				pR 			const 	coordinates_)
{
  for (unsigned int dimIndex = 0;dimIndex < _dim;++dimIndex)
    {
      coordinates_[dimIndex] = this->m_pts[ctrlPointIndex_*_dim + dimIndex];
    }
}

#include "BSpline.h"
#include "Blas.h"
#include <math.h>
#include "Matrix.h"
#include "FiniteElementSegmentHermite.h"
#include "SparseSymbolic.h"
extern long int random();


template <unsigned int _dim>
void	BSpline<_dim>::GetDofGeometry(const int_t edgeIndex_,
				      pR 	const	dofGeometry_,
				      const int_t 	dofGeometryOff_)
{
  for (unsigned int dimIndex = 0;dimIndex < _dim;++dimIndex)
    {      
      dofGeometry_[dofGeometryOff_ * dimIndex + 0] = self_->m_pts[edgeIndex*_dim + dimIndex];
      dofGeometry_[dofGeometryOff_ * dimIndex + 1] = self_->m_pts[(edgeIndex+1)*_dim + dimIndex];
      dofGeometry_[dofGeometryOff_ * dimIndex + 2] = self_->m_derivatives[edgeIndex*_dim + dimIndex];
      dofGeometry_[dofGeometryOff_ * dimIndex + 3] = self_->m_derivatives[(edgeIndex+1)*_dim + dimIndex];
    } 
}

template <unsigned int _dim>
void	BSpline<_dim>::~BSpline()
{
      if (self_->m_ctrlpts)
	{
	  free(self_->m_ctrlpts);
	  self_->m_ctrlpts = NULL;
	}
      if (self_->m_pts)
	{
	  free(self_->m_pts);
	  self_->m_pts = NULL;
	}
      if (self_->m_derivatives)
	{
	  free(self_->m_derivatives);
	  self_->m_derivatives = NULL;
	}
      if (self_->m_normalized_length)
	{
	  free(self_->m_normalized_length);
	  self_->m_normalized_length = NULL;
	}
}

L BSplineReadOnly_findEdge(pBSplineReadOnly 	const 	self_,
			   const R 			t_,
			   pI 			const	iedge_,
			   pR 			const	s_)
{ 
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(iedge_);
  DebugVerif(s_);
#endif
  const I numPoints = BSplineReadOnly_get_numPoints(self_);
  { I k;
    for (k=1;k<numPoints;++k)
      {
	if (self_->m_normalized_length[k]>=t_)
	  {
	    break;	    
	  }
      }
    if (k<numPoints)
      {
	iedge_[0] 	= k-1;
	s_[0]		= ( (t_ - self_->m_normalized_length[k-1])/(self_->m_normalized_length[k]-self_->m_normalized_length[k-1]));	
	return __emnsYES;
      }
    else
      {
	return __emnsNO;
      } }
}


#define hermite_interval_dtt(_n,_r,_roff,_p,_poff)			\
  { I _i;								\
    for (_i=0;_i<(_n);++_i)						\
      { const R s  = ((_p)[_i*(_poff)]+1.0)*0.5;			\
	(_r)[(_roff)*_i+0] = (s*12.0-6.0);				\
	(_r)[(_roff)*_i+1] = (6.0-s*12.0);				\
	(_r)[(_roff)*_i+2] = (s*6.0-4.0);				\
	(_r)[(_roff)*_i+3] = (s*6.0-2.0);				\
      }									\
  }

#if 0
#define hermite_interval_dt(_n,_r,_roff,_p,_poff)	\
  { I _i;						\
    for (_i=0;_i<(_n);++_i)				\
      {							\
	const R s  = ((_p)[_i*(_poff)]+1.0)*0.5;	\
	(_r)[(_roff)*_i+0] = (s-1.0)*s*6.0*0.5;		\
	(_r)[(_roff)*_i+1] = (1.0-s)*s*6.0*0.5;		\
	(_r)[(_roff)*_i+2] = (s*(s*3.0-4.0)+1.0)*0.5;	\
	(_r)[(_roff)*_i+3] = s*(s*3.0-2.0)*0.5;		\
      }							\
  }
#endif

#define hermite_interval_dt(_n,_r,_roff,_p,_poff)	\
  { I _i;						\
    for (_i=0;_i<(_n);++_i)				\
      {							\
	const R s  = ((_p)[_i*(_poff)]+1.0)*0.5;	\
	(_r)[(_roff)*_i+0] = (s-1.0)*s*6.0;		\
	(_r)[(_roff)*_i+1] = (1.0-s)*s*6.0;		\
	(_r)[(_roff)*_i+2] = (s*(s*3.0-4.0)+1.0);	\
	(_r)[(_roff)*_i+3] = s*(s*3.0-2.0);		\
      }							\
  }


#define hermite_interval(_n,_r,_roff,_p,_poff)				\
  { I _i;								\
    for (_i=0;_i<(_n);++_i)						\
      {									\
	const R s  		= ((_p)[_i*(_poff)]+1.0)*0.5;		\
	const R s2 		= s*s;					\
	(_r)[(_roff)*_i+0] 	= (s*2.0-3.0)*s2+1.0;			\
	(_r)[(_roff)*_i+1] 	= (s*(-2.0)+3.0)*s2;			\
	(_r)[(_roff)*_i+2] 	= s*(s-1.0)*(s-1.0);			\
	(_r)[(_roff)*_i+3] 	= s2*(s-1.0);				\
      }									\
  }







void Eval	(const double 	s_,
		 const I	iedge_,
		 pR 		position_) const noexcept
{

  R hermite[4];
  hermite_interval(1,hermite,1,s_,1);   
  /* calcul de la tangente */
  for (unsigned int dimIndex = 0;dimIndex < _dim;++dimIndex)
    {      
      position_[dimIndex_]  = 
	hermite[0]   * this->m_pts[iedge_*dim+idim]
	+ hermite[1] * this->m_pts[(iedge_+1)*dim+idim]
	+ hermite[2] * this->m_derivatives[iedge_*dim+idim]
	+ hermite[3] * this->m_derivatives[(iedge_+1)*dim+idim];
    }   
}

void EvalSecondDerivative(const double s_,
			  const int_t			iedge_,
			  pR 			const	acc_)
{
  if (iedge_>=self_->m_numPoints-1)
    {
      printf("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n");
      exit(1);
    }
  cst_eDim dim = BSplineReadOnly_get_spatialDimension(self_);
  R hermite[4];
  hermite_interval_dtt(1,hermite,1,s_,1);   
  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim;++idim)
      {
	acc_[idim]  = 
	  hermite[0]   * self_->m_pts[iedge_*dim+idim]
	  + hermite[1] * self_->m_pts[(iedge_+1)*dim+idim]
	  + hermite[2] * self_->m_derivatives[iedge_*dim+idim]
	  + hermite[3] * self_->m_derivatives[(iedge_+1)*dim+idim];
      } } 
}

void EvalFirstDerivative(cst_pR		const	s_,
			 const I 		iedge_,
			 pR 		const	normal_)
{
  R hermite[4];
  hermite_interval_dt(1,hermite,1,s_,1);   
  
  /* calcul de la tangente */
  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim;++idim)
      {

	normal_[idim]  = 
	  hermite[0]   * self_->m_pts[iedge_*dim+idim]
	  + hermite[1] * self_->m_pts[(iedge_+1)*dim+idim]
	  + hermite[2] * self_->m_derivatives[iedge_*dim+idim]
	  + hermite[3] * self_->m_derivatives[(iedge_+1)*dim+idim];
      } } 

}

void BSpline<_dim>::BSpline(const I 		numPoints_,
			    cst_pR 	const 	ctrlpts_,
			    const I 		ctrlptsoff_)
{
  static constexpr const	R 	oneOnSix 		= ((R)1.0)/((R)6.0);
  static constexpr const 	R	twoOnThree 		= ((R)2.0)/((R)3.0);
  static constexpr const 	R 	oneOnTwo		= ((R)0.5);
  static constexpr const 	I 	nbCubaturePoints 	= ((I)10);
  static constexpr const 	I 	nequal1			= ((I)1);
  static constexpr const 	I 	nequal3 		= ((I)3);
  static constexpr const 	I 	nequal4 		= ((I)4);
  static constexpr const 	R 	requal1 		= ((R)1.0);
  static constexpr const 	R 	requal0 		= ((R)0.0);
  const 	I 	numEdges		= numPoints_-1;
  const 	I 	numEdges_X_dim 		= numEdges * dim_;
  R 			coo[__eDim_ALL];  
  R 			tmp[nbCubaturePoints];

  BSpline_clear(self_);

  self_->m_dim 					= dim_;
  self_->m_numPoints 				= numPoints_;
  self_->m_pts 					= (pR)malloc(sizeof(R)*numPoints_*dim_);
  self_->m_derivatives				= (pR)malloc(sizeof(R)*numPoints_*dim_);
  self_->m_ctrlpts					= (pR)malloc(sizeof(R)*numPoints_*dim_);
  self_->m_normalized_length			= (pR)malloc(sizeof(R)*numPoints_);
  self_->m_normalized_length[0] 			= ((R)0.0);
  self_->m_normalized_length[numEdges] 		= ((R)1.0);
  self_->m_length				= ((R)0.0);

  pMatrix dofValues 	= Matrix_malloc		(nequal4,
						 (numEdges+1)*dim_);
  
  pMatrix evalCurve 	= Matrix_malloc		(nbCubaturePoints,
						 numEdges_X_dim);
  
  pMatrix evalDtHermite = Matrix_malloc		(nbCubaturePoints,
						 nequal4);
  



  /*########################################################################*/

  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim_;++idim)
      {
	Blas_dcopy(&numPoints_,&ctrlpts_[idim],&ctrlptsoff_,&self_->m_ctrlpts[idim],&dimension);
      } }

  /*########################################################################*/
  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim_;++idim)
      {
	self_->m_pts[idim] = ctrlpts_[idim];
      } }
  /**/
  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim_;++idim)
      {
	self_->m_derivatives[idim] = (ctrlpts_[1*ctrlptsoff_+idim]-ctrlpts_[idim] );
      } }
  /**/
  { I i;
    for (i=1;i<numEdges;++i)
      {
	{ eDim idim = __eDim_ERROR;
	  for (idim = __eDim_ERROR;idim<dim_;++idim)
	    {
	      self_->m_pts[i*dim_+idim] 		= oneOnSix * ( ctrlpts_[(i-1)*ctrlptsoff_+idim]  + ctrlpts_[(i+1)*ctrlptsoff_+idim] ) + twoOnThree * ctrlpts_[i*ctrlptsoff_+idim];
	    } }
      } }
  { I i;
    for (i=1;i<numEdges;++i)
      {
	{ eDim idim = __eDim_ERROR;
	  for (idim = __eDim_ERROR;idim<dim_;++idim)
	    {
	      self_->m_derivatives[i*dim_+idim] 	= (ctrlpts_[(i+1)*ctrlptsoff_+idim]-ctrlpts_[(i-1)*ctrlptsoff_+idim])*oneOnTwo;
	    } }
      } }
  
  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim_;++idim)
      {
	self_->m_pts[numEdges*dim_+idim] = ctrlpts_[numEdges*ctrlptsoff_+idim];
      } }
  
  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim_;++idim)
      {
	self_->m_derivatives[numEdges*dim_+idim] = (ctrlpts_[(numPoints_-1)*ctrlptsoff_+idim]-ctrlpts_[(numPoints_-2)*ctrlptsoff_+idim]);
      } }
  /*########################################################################*/
  
  static R w1d[]={
    0.295524224714752870173892994651338329421046717026853601354308029755995938217152329270356595793754216722717164401252558386818490789552005826001936342494186966609562718648884168043231305061535867409083051270663865287483901746874726597515954450775158914556548308329986393605934912382356670244,
    0.295524224714752870173892994651338329421046717026853601354308029755995938217152329270356595793754216722717164401252558386818490789552005826001936342494186966609562718648884168043231305061535867409083051270663865287483901746874726597515954450775158914556548308329986393605934912382356670244,
    0.2692667193099963550912269215694693528597599384608837958005632762421534323191792767642266367092527607555958114503686983086929234693811452415564658846634423711656014432259960141729044528030344411297902977067142537534806284608399276575006911686749842814086288868533208042150419508881916391898,
    0.2692667193099963550912269215694693528597599384608837958005632762421534323191792767642266367092527607555958114503686983086929234693811452415564658846634423711656014432259960141729044528030344411297902977067142537534806284608399276575006911686749842814086288868533208042150419508881916391898,
    0.2190863625159820439955349342281631924587718705226770898809565436351999106529512812426839931772021927865912168728128876347666269080669475688309211843316656677105269915322077536772652826671027878246851010208832173320064273483254756250668415885349420711613410227291565477768928313300688702802,
    0.2190863625159820439955349342281631924587718705226770898809565436351999106529512812426839931772021927865912168728128876347666269080669475688309211843316656677105269915322077536772652826671027878246851010208832173320064273483254756250668415885349420711613410227291565477768928313300688702802,
    0.1494513491505805931457763396576973324025566396694273678354772687532386547266300109459472646347319519140057525610454363382344517067454976014713716011937109528798134828865118770953566439639333773939909201690204649083815618779157522578300343427785361756927642128792412282970150172590842897331,
    0.1494513491505805931457763396576973324025566396694273678354772687532386547266300109459472646347319519140057525610454363382344517067454976014713716011937109528798134828865118770953566439639333773939909201690204649083815618779157522578300343427785361756927642128792412282970150172590842897331,
    0.066671344308688137593568809893331792857864834320158145128694881613412064084087101776785509685058877821090054714520419331487507126254403762139304987316994041634495363706400187011242315504393526242450629832718198718647480566044117862086478449236378557180717569208295026105115288152794421677,
    0.066671344308688137593568809893331792857864834320158145128694881613412064084087101776785509685058877821090054714520419331487507126254403762139304987316994041634495363706400187011242315504393526242450629832718198718647480566044117862086478449236378557180717569208295026105115288152794421677};

  static R p1d[]={
    -0.1488743389816312108848260011297199846175648594206916957079892535159036173556685213711776297994636912300311608052553388261028901818643765402316761969968090913050737827720371059070942475859422743249837177174247346216914852902942929003193466659082433838094355075996833570230005003837280634351,
    0.1488743389816312108848260011297199846175648594206916957079892535159036173556685213711776297994636912300311608052553388261028901818643765402316761969968090913050737827720371059070942475859422743249837177174247346216914852902942929003193466659082433838094355075996833570230005003837280634351,
    -0.4333953941292471907992659431657841622000718376562464965027015131437669890777035012251027579501177212236829350409989379472742247577232492051267741032822086200952319270933462032011328320387691584063411149801129823141488787443204324766414421576788807708483879452488118549797039287926964254222,
    0.4333953941292471907992659431657841622000718376562464965027015131437669890777035012251027579501177212236829350409989379472742247577232492051267741032822086200952319270933462032011328320387691584063411149801129823141488787443204324766414421576788807708483879452488118549797039287926964254222,
    -0.6794095682990244062343273651148735757692947118348094676648171889525585753950749246150785735704803794998339020473993150608367408425766300907682741718202923543197852846977409718369143712013552962837733153108679126932544954854729341324727211680274268486617121011712030227181051010718804444161,
    0.6794095682990244062343273651148735757692947118348094676648171889525585753950749246150785735704803794998339020473993150608367408425766300907682741718202923543197852846977409718369143712013552962837733153108679126932544954854729341324727211680274268486617121011712030227181051010718804444161,
    -0.8650633666889845107320966884234930485275430149653304525219597318453747551380555613567907289460457706944046310864117651686783001614934535637392729396890950011571349689893051612072435760480900979725923317923795535739290595879776956832427702236942765911483643714816923781701572597289139322313,
    0.8650633666889845107320966884234930485275430149653304525219597318453747551380555613567907289460457706944046310864117651686783001614934535637392729396890950011571349689893051612072435760480900979725923317923795535739290595879776956832427702236942765911483643714816923781701572597289139322313,
    -0.9739065285171717200779640120844520534282699466923821192312120666965952032346361596257235649562685562582330425187742112150221686014344777799205409587259942436704413695764881258799146633143510758737119877875210567067452435368713683033860909388311646653581707125686970668737259229449284383797,
    0.9739065285171717200779640120844520534282699466923821192312120666965952032346361596257235649562685562582330425187742112150221686014344777799205409587259942436704413695764881258799146633143510758737119877875210567067452435368713683033860909388311646653581707125686970668737259229449284383797};

  { I ierr;
    FiniteElementSegmentHermite_EvalDt("N",
				       &nbCubaturePoints,
				       p1d,
				       &nequal1,
				       Matrix_at(evalDtHermite,0,0),
				       &nequal4,
				       &ierr); }

  //  hermite_interval_dt(nbCubaturePoints,Matrix_at(evalDtHermite,0,0),4,p1d,1);     
  
  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim_;++idim)
      {
	Blas_dcopy(&numEdges,&self_->m_pts[idim],&dimension,Matrix_at(dofValues,0,idim),&nequal4);
      } }
  
  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim_;++idim)
      {
	Blas_dcopy(&numEdges,&self_->m_pts[dim_+idim],&dimension,Matrix_at(dofValues,1,idim),&nequal4);
      } }
  
  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim_;++idim)
      {
	Blas_dcopy(&numEdges,&self_->m_derivatives[idim],&dimension,Matrix_at(dofValues,2,idim),&nequal4);
      } }

  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim_;++idim)
      {
	Blas_dcopy(&numEdges,&self_->m_derivatives[dim_+idim],&dimension,Matrix_at(dofValues,3,idim),&nequal4);
      } }
  
  Matrix_MM		(evalDtHermite,
			 &requal1,
			 dofValues,
			 &requal0,
			 evalCurve);

  R length = ((R)0.0);

  { I iedge;
    for (iedge=0;iedge<numEdges;++iedge)
      {
	{ I i;
	  for (i=0;i<nbCubaturePoints;++i)
	    {
	      { eDim idim = __eDim_ERROR;
		for (idim = __eDim_ERROR;idim<dim_;++idim)
		  {
	            cst_pR c = Matrix_at(evalCurve,idim *nbCubaturePoints  + i,iedge);
		    coo[idim] = c[0];
		  } } 


	      { R s = ((R)0.0);
		{ eDim idim = __eDim_ERROR;
		  for (idim = __eDim_ERROR;idim<dim_;++idim)
		    {
		      s+=coo[idim]*coo[idim];
		    } }
		tmp[i] = nsSQRT(s);
	      }
	    } 
	  const R iedge_length 			= Blas_ddot(&nbCubaturePoints,w1d,&nequal1,tmp,&nequal1);
	  length 				+= iedge_length*((R)0.5);
	  self_->m_normalized_length[iedge+1] 	= length;
	}
      } }

  dofValues = Matrix_kill(dofValues);
  evalCurve = Matrix_kill(evalCurve);

  printf("length = "rfmt"\n",length);
  self_->m_length = length;
  { I iedge;
    for (iedge=1;iedge<numPoints_;++iedge)
      {
	self_->m_normalized_length[iedge] /= self_->m_length;
      } }

#if 1
  /* 
     POUR UNIFORMISER 
  */
  { I iedge;
    for (iedge=1;iedge<numPoints_;++iedge)
      {
	self_->m_normalized_length[iedge] = ((R)iedge)/((R)numPoints_-1);
      } }
#endif
  
#if 1

  /*
    -1  -1/2 0   0   0   0
    -1  -1   1   0   0   0
    0  -1  -1   1   0   0
    0   0  -1  -1   1   0
    0   0   0  -1  -1   1 
    0   0   0   0  1/2  1 


    P_{i-1} * B0(1) + P_{i} * B1(1) + T_{i-1} * B2(1) + T_{i} * B3(1)
    = 
    P_{i} * B0(0) + P_{i+1} * B1(0) + T_{i} * B2(0) + T_{i+1} * B3(0)
    
    T_{i-1} * B2(1) + T_{i} * (B3(1) -  B2(0)) - T_{i+1} * B3(0)
    = 
    - P_{i-1} * B0(1)  + P_{i} * ( B0(0) - B1(1) ) + P_{i+1} * B1(0) 

    T_{0} * B2(0) + T_{1} * B3(0) = - P_{0} * B0(0) - P_{1} * B1(0) 
    T_{n-1} * B2(1) + T_{n} * B3(1) = - P_{n-1} * B0(1) - P_{n} * B1(1) 
  */

  
  R ee[1];
  R coeff0[4];
  R coeff1[4];
  ee[0] = -1.0;


  { I ierr;
    FiniteElementSegmentHermite_EvalDDt("N",
					&nequal1,
					ee,
					&nequal1,
					coeff0,
					&nequal4,
					&ierr); }
  

  ee[0] = 1.0;

  { I ierr;
    FiniteElementSegmentHermite_EvalDDt("N",
					&nequal1,
					ee,
					&nequal1,
					coeff1,
					&nequal4,
					&ierr); }

  printf("coeff0 %e %e %e %e\n",coeff0[0],coeff0[1],coeff0[2],coeff0[3]);
  printf("coeff1 %e %e %e %e\n",coeff1[0],coeff1[1],coeff1[2],coeff1[3]);
  
  pR A 		= calloc(numPoints_*numPoints_,sizeof(R));
  pR B 		= calloc(numPoints_*3,sizeof(R));
  
  /* T_{0} * B2(0) + T_{1} * B3(0) = - P_{0} * B0(0) - P_{1} * B1(0) */
  A[0] 			= coeff0[2];
  A[numPoints_] 	= coeff0[3];
  B[0] 			= -self_->m_pts[0] * coeff0[0] - self_->m_pts[dim_+0] * coeff0[1];
  B[numPoints_+0] 	= -self_->m_pts[1] * coeff0[0] - self_->m_pts[dim_+1] * coeff0[1];
  B[2*numPoints_+0]	= -self_->m_pts[2] * coeff0[0] - self_->m_pts[dim_+2] * coeff0[1];
  
  { I i;
    for (i=1;i<numPoints_-1;++i)
      {
	/*T_{i-1} * B2(1) + T_{i} * (B3(1) -  B2(0)) - T_{i+1} * B3(0)*/
	A[(i-1)*numPoints_+i] 	= coeff1[2];
	A[i*numPoints_+i] 	= coeff1[3]-coeff0[2];
	A[(i+1)*numPoints_+i] 	= -coeff0[3];
	/* - P_{i-1} * B0(1)  + P_{i} * ( B0(0) - B1(1) ) + P_{i+1} * B1(0) */
	B[i] 		= -self_->m_pts[ (i-1)*dim_ + 0 ] * coeff1[0] + self_->m_pts[ (i+1)*dim_ + 0 ] * coeff0[1] + self_->m_pts[ i*dim_ + 0 ] * (coeff0[0]-coeff1[1]);
	B[numPoints_+i] 	= -self_->m_pts[ (i-1)*dim_ + 1 ] * coeff1[0] + self_->m_pts[ (i+1)*dim_ + 1 ] * coeff0[1] + self_->m_pts[ i*dim_ + 1 ] * (coeff0[0]-coeff1[1]);
	B[2*numPoints_+i]	= -self_->m_pts[ (i-1)*dim_ + 2 ] * coeff1[0] + self_->m_pts[ (i+1)*dim_ + 2 ] * coeff0[1] + self_->m_pts[ i*dim_ + 2 ] * (coeff0[0]-coeff1[1]);
      } }

  /*T_{n-1} * B2(1) + T_{n} * B3(1) = - P_{n-1} * B0(1) - P_{n} * B1(1) */  
  A[(numPoints_-2)*numPoints_+numPoints_-1] 	= coeff1[2];
  A[(numPoints_-1)*numPoints_+numPoints_-1] 	= coeff1[3];  

  B[numPoints_-1] 	= -self_->m_pts[(numPoints_-2)*dim_+0] * coeff1[0] - self_->m_pts[(numPoints_-1)*dim_+0] * coeff1[1];
  B[numPoints_+numPoints_-1] 	= -self_->m_pts[(numPoints_-2)*dim_+1] * coeff1[0] - self_->m_pts[(numPoints_-1)*dim_+1] * coeff1[1];
  B[2*numPoints_+numPoints_-1]	= -self_->m_pts[(numPoints_-2)*dim_+2] * coeff1[0] - self_->m_pts[(numPoints_-1)*dim_+2] * coeff1[1];


  { I i;
    for (i=0;i<numPoints_;++i)
      {
	{ I j;
	  for (j=0;j<numPoints_;++j)
	    {
	      printf(" %1.1f",A[j*numPoints_+i]);
	    } }
	printf("\n");
      } }


  I info_lapack  = 0;
  pI perm = malloc(sizeof(I)*numPoints_);
  { I i;
    for (i=0;i<numPoints_;++i)
      {
	printf("%e %e %e\n",B[i],B[numPoints_+i],B[2*numPoints_+i]);
      } }
  
  { I i;
    for (i=0;i<numPoints_-1;++i)
      {
	R c1;
	R c2;
	c1 = coeff0[0]*self_->m_pts[i*dim_+0]
	  +coeff0[1]*self_->m_pts[(i+1)*dim_+0]
	  +coeff0[2]*self_->m_derivatives[(i)*dim_+0]
	  +coeff0[3]*self_->m_derivatives[(i+1)*dim_+0];
	c2 = coeff1[0]*self_->m_pts[i*dim_+0]
	  +coeff1[1]*self_->m_pts[(i+1)*dim_+0]
	  +coeff1[2]*self_->m_derivatives[(i)*dim_+0]
	  +coeff1[3]*self_->m_derivatives[(i+1)*dim_+0];
	printf("AVANT %e %e\n",c1,c2);
      } }

  dgesv(&numPoints_,
	&nequal3,
	A,
	&numPoints_,
	perm,
	B,
	&numPoints_,
	&info_lapack);   

  printf("info lapack "ifmt"\n",info_lapack);

  { I i;
    for (i=0;i<numPoints_;++i)
      {
	self_->m_derivatives[i*dim_+0] = B[i];
	self_->m_derivatives[i*dim_+1] = B[numPoints_+i];
	self_->m_derivatives[i*dim_+2] = B[2*numPoints_+i];
	printf("%e %e %e\n",B[i],B[numPoints_+i],B[2*numPoints_+i]);
      } }

  { I i;
    for (i=0;i<numPoints_-1;++i)
      {
	R c1;
	R c2;
	c1 = coeff0[0]*self_->m_pts[i*dim_+0]
	  +coeff0[1]*self_->m_pts[(i+1)*dim_+0]
	  +coeff0[2]*self_->m_derivatives[(i)*dim_+0]
	  +coeff0[3]*self_->m_derivatives[(i+1)*dim_+0];
	c2 = coeff1[0]*self_->m_pts[i*dim_+0]
	  +coeff1[1]*self_->m_pts[(i+1)*dim_+0]
	  +coeff1[2]*self_->m_derivatives[(i)*dim_+0]
	  +coeff1[3]*self_->m_derivatives[(i+1)*dim_+0];
	printf("APRES %e %e\n",c1,c2);
      } }

  free(perm);
  free(A);
  free(B);
#endif  

}



void BSpline_defInterpolatory(pBSpline 	const	self_,
			      cst_eDim 		dim_,
			      const I 		numPoints_,
			      cst_pR 	const 	ctrlpts_,
			      const I 		ctrlptsoff_)
{
#ifndef NDEBUG
  DebugVerif(numPoints_>1);
  DebugVerif(ctrlpts_);
  DebugVerif(ctrlptsoff_>=dim_);
#endif
  static const 	I 	nequal1			= ((I)1);
  static const 	I 	nequal3 		= ((I)3);
  const 	I 	dimension 		= dim_;
  const 	I 	numEdges		= numPoints_-1;


  BSpline_clear(self_);

  self_->m_dim 					= dim_;
  self_->m_numPoints 				= numPoints_;
  self_->m_pts 					= (pR)malloc(sizeof(R)*numPoints_*dim_);
  self_->m_derivatives				= (pR)malloc(sizeof(R)*numPoints_*dim_);
  self_->m_ctrlpts				= (pR)malloc(sizeof(R)*numPoints_*dim_);
  self_->m_normalized_length			= (pR)malloc(sizeof(R)*numPoints_);
  self_->m_normalized_length[0] 			= ((R)0.0);
  self_->m_normalized_length[numEdges] 		= ((R)1.0);
  self_->m_length				= ((R)0.0);

  /*########################################################################*/

  { eDim idim = __eDim_ERROR;
    for (idim = __eDim_ERROR;idim<dim_;++idim)
      {
	Blas_dcopy(&numPoints_,&ctrlpts_[idim],&ctrlptsoff_,&self_->m_ctrlpts[idim],&dimension);
      } }
  
  { I arrayLength = numPoints_ * dim_;
    Blas_dcopy(&arrayLength,self_->m_ctrlpts,&nequal1,self_->m_pts,&nequal1); } 
  
  /*########################################################################*/
  pR A = calloc(numPoints_*numPoints_,sizeof(R));
  A[0] = ((R)2.0);
  A[numPoints_] = ((R)1.0);
  { I pointIndex = 0;
    const I numPoints_minus1 = numPoints_-1;
    for (pointIndex=1;pointIndex < numPoints_minus1;++pointIndex)
      {
	A[(pointIndex-1) * numPoints_ + pointIndex] = ((R)1.0);
	A[pointIndex * numPoints_ + pointIndex] = ((R)4.0);
	A[(pointIndex+1) * numPoints_ + pointIndex] = ((R)1.0);
      } }
  A[(numPoints_-2) * numPoints_ + numPoints_-1] = ((R)1.0);
  A[(numPoints_-1) * numPoints_ + numPoints_-1] = ((R)2.0);

  pI perm = malloc(sizeof(I)*numPoints_);
  pR B 	  = calloc(numPoints_*3,sizeof(R));

  { I pointIndex = 0;
    const I numPoints_minus1 = numPoints_-1;
    I coordinateIndex;
    for (coordinateIndex=0;coordinateIndex<dim_;++coordinateIndex)
      {
	pointIndex = 0;
	B[numPoints_*coordinateIndex + pointIndex] = 3.0*( self_->m_pts[dim_*(pointIndex+1) + coordinateIndex] - self_->m_pts[dim_*pointIndex + coordinateIndex]  );

	for (pointIndex=1;pointIndex < numPoints_minus1;++pointIndex)
	  {
	    B[numPoints_*coordinateIndex + pointIndex] = 3.0*( self_->m_pts[dim_*(pointIndex+1) + coordinateIndex] - self_->m_pts[dim_*(pointIndex-1) + coordinateIndex]  );
	  }
	
	pointIndex = numPoints_-1;
	B[numPoints_*coordinateIndex + pointIndex] = 3.0*( self_->m_pts[dim_*(pointIndex) + coordinateIndex] - self_->m_pts[dim_*(pointIndex-1) + coordinateIndex]  );
	
      } }

  I info_lapack;
  dgesv(&numPoints_,
	&nequal3,
	A,
	&numPoints_,
	perm,
	B,
	&numPoints_,
	&info_lapack);   

  {
    
    I coordinateIndex,pointIndex;
    for (pointIndex=0;pointIndex < numPoints_;++pointIndex)
      {
	for (coordinateIndex=0;coordinateIndex<dim_;++coordinateIndex)
	  {
	    self_->m_derivatives[dim_*pointIndex + coordinateIndex] = B[numPoints_*coordinateIndex + pointIndex];
	  }
      }
  }
  
}

