
include $(DISTRIB)/config.mk

#HEADERS=-I$(DISTRIB)/../Common/include -I$(DISTRIB)/../../include -I$(DISTRIB)/include  -I$(DISTRIB)/libExtern/metis/Lib/ #-I$(DISTRIB)/libGrid/include

HEADERS= -I$(DISTRIB)/include -I$(DISTRIB)/libMeshb-master/build/include
LIBMETIS=  -L$(DISTRIB)/libExtern/metis/ -lmetis
LIBSPARSEKIT2=  #-L$(DISTRIB)/libExtern/sparsekit2/ -lskit

#LIBS=-L$(DISTRIB)/lib/ -lMns_$(CC)$(DEBUG)  $(LIBMETIS) $(LIBSPARSEKIT2) # -L$(DISTRIB)/libGrid/lib -lGrid_$(CC)$(DEBUG)	

#LIBS= -L$(DISTRIB)/lib/ -lMns_$(CC)$(DEBUG) #-lMeshb.7

LIBS= -lMeshb.7

#-L$(DISTRIB)/libMeshb-master/build/lib -lMeshb.7

#-L$(DISTRIB)/../Common/lib/ -lSLCommon_$(CC)$(DEBUG)

OBJ_MAIN=	\
		main/MnsBuildManpage.o \
		main/MnsSplit.o \
		main/MnsClean.o  \
		main/Mns.o \
		main/MnsLoad.o \
		\
		main/MnsReader.o \
		main/MnsMeditMerge.o \
		main/MnsCheckMesh.o \
		main/TestPartitioning.o \
		main/TestThreadPool.o \
		main/TestDg.o \
		main/TestExpression.o \
		main/TestOctree.o \
		main/TestGeometryReader.o \
		main/TestQuadtree.o \
		main/TestDelaunay.o \
		main/TestVoronoi.o \
		main/TestHighDefinition.o

# working dirs
TARGETDIR=main
SRCDIR=src
INCDIR=include
OBJDIR=$(PLATFORM)_$(CC)$(DEBUG)
DIRDIR=$(OBJDIR)
VPATH=$(SRCDIR)


targetsrc0=$(wildcard main/*.cpp) 
targetobjs0=$(patsubst $(TARGETDIR)%,$(OBJDIR)%,$(targetsrc0:.cpp=.o)) 
target0=$(targetobjs0:.o=.exe) 

#$(wildcard main/*.cpp) 
targetsrc=$(filter-out $(targetmpisrc), $(wildcard main/*.cpp))
targetobjs=$(patsubst $(TARGETDIR)%,$(OBJDIR)%,$(targetsrc:.c=.o)) 
target=$(targetobjs:.o=.exe) $(target0)

# objects list


##################################################################################
## FLEX/BISON C
wildflex_c=$(wildcard $(SRCDIR)/*.flex.c.lex)
srcflex_c=$(patsubst $(SRCDIR)%,$(OBJDIR)%,$(wildflex_c:.flex.c.lex=.flex.c))
objflex_c=$(srcflex_c:.c=.o)

wildbison_c=$(wildcard $(SRCDIR)/*.bison.c.y)
srcbison_c=$(patsubst $(SRCDIR)%,$(OBJDIR)%,$(wildbison_c:.bison.c.y=.bison.c))
objbison_c=$(srcbison_c:.c=.o)

##################################################################################
## FLEX/BISON CPP
srcflex=$(wildcard $(SRCDIR)/*.flex.lex)
cppflex=$(patsubst $(SRCDIR)%,$(OBJDIR)%,$(srcflex:.flex.lex=.flex.cpp))
objflex=$(cppflex:.cpp=.o)

srcbison=$(wildcard $(SRCDIR)/*.bison.y)
cppbison=$(patsubst $(SRCDIR)%,$(OBJDIR)%,$(srcbison:.bison.y=.bison.cpp))
objbison=$(cppbison:.cpp=.o)




c_src    = $(wildcard $(SRCDIR)/*.c)
cpp_src  = $(wildcard $(SRCDIR)/*.cpp)
f_src    = $(wildcard $(SRCDIR)/*.f)
header = $(wildcard $(INCDIR)/*.h $(SRCDIR)/*.h)
objs   = $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(c_src:.c=.o)) $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(cpp_src:.cpp=.o)) $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(f_src:.f=.o))  $(objflex) $(objbison) $(objbison_c) $(objflex_c) $(qtobj)
srcdep = $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(c_src:.c=.c)) 


#$(objflex) $(objbison) 


all:$(DIRDIR) $(cppflex) $(cppbison) $(srcflex_c) $(srcbison_c)  $(OBJDIR)/libMns_$(CC)$(DEBUG).a $(targetmpiobjs) $(targetmpi) $(targetobjs) $(target)   
#all: $(targetmpiobjs) $(targetmpi)



#$(cppflex) $(cppbison) 

#$(OBJ_MAIN) $(OBJ_MAIN:%.o=%.exe) mns_meshbuilder.exe mns_spacebuilder.exe #test.exe test_surface.exe 



dep:
	makedepend -f $(DISTRIB)/Makefile.version $(HEADERS) $(srcdep)


$(OBJDIR)/%.o: $(OBJDIR)/%.c
	@echo $@
	@$(CC) $(CFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 




$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@echo $@
	@$(CC) $(CFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo $@
	$(CPP) $(CPPFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 


$(OBJDIR)/%.o: $(TARGETDIR)/%.c
	@echo $@
	@$(CC) $(CFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 

$(OBJDIR)/%.o: $(TARGETDIR)/%.cpp
	@echo $@
	$(CPP) $(CPPFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 

$(OBJDIR)/%.o: $(SRCDIR)/%.f
	@echo $@
	@$(F77) $(FFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 

$(OBJDIR)/%.mpi.o: $(TARGETDIR)/%.mpi.cpp
	@echo $@
	mpic++ $(CPPFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 

$(OBJDIR)/%.o: $(OBJDIR)/%.cpp
	@echo $@
	@$(CPP) $(CPPFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 



$(OBJDIR)/%.moc:$(SRCDIR)/%.cpp
	@echo $@
	@moc-qt4 $< -o $@ -DQ_MOC_RUN -I$(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)


$(OBJDIR)/%.exe:$(OBJDIR)/%.o
	@echo $@
	$(CPP) $(LDFLAGS) $(FLAGS) $< -o $@ $(LIBS)
	ln -fs $(DISTRIB)/$@ $(DISTRIB)/bin/$*$(DEBUG)

$(OBJDIR)/%.mpi.exe:$(OBJDIR)/%.mpi.o
	@echo $@
	mpic++ $(LDFLAGS) $(FLAGS) $< -o $@ $(LIBS)
	ln -fs $(DISTRIB)/$@ $(DISTRIB)/bin/$*$(DEBUG)


#-L$(DISTRIB)/lib  -lMns_$(CC)$(DEBUG) -lcadnaC $(LIBS)  $(LIBMETIS)  $(LIBBLAS) 
#-L$(DISTRIB)/lib  -lMns_$(CC)$(DEBUG) -lcadnaC $(LIBS)  $(LIBMETIS)  $(LIBBLAS) 

$(OBJDIR)/%.bison.cpp : $(SRCDIR)/%.bison.y $(SRCDIR)/%.flex.lex
	printf "install:C:$(MK_DEBUG_TOKEN):bison\t"`basename $<`"\n" 
	bison -k -d -p $*_  $< -o $*.tab.cpp
	mv $*.tab.hpp $*.h
	mv $*.tab.cpp $@
	mv $*.h $(DISTRIB)/include


$(OBJDIR)/%.bison.c : $(SRCDIR)/%.bison.c.y $(SRCDIR)/%.flex.c.lex
#	printf "install:C:$(MK_DEBUG_TOKEN):bison\t"`basename $<`"\n" 
	bison -k -d -p $*_  $< -o $*.tab.c
	mv $*.tab.h $*.h
	mv $*.tab.c $@
	mv $*.h $(DISTRIB)/include


$(OBJDIR)/%.flex.cpp : $(SRCDIR)/%.flex.lex $(SRCDIR)/Reader.bison.y
	flex -P$*_  -o$@ $< 


$(OBJDIR)/%.flex.c : $(SRCDIR)/%.flex.c.lex #$(SRCDIR)/Reader.bison.c.y
	flex -P$*_  -o$@ $< 

#


# -I/opt/Qt/5.3/gcc_64/include/QtWidgets/  -I/opt/Qt/5.3/gcc_64/include/ -I/opt/Qt/5.3/gcc_64/include/QtGui -I/opt/Qt/5.3/gcc_64/include/QtCore -fPIC



#$(objs): $(header)
$(target): $(header) $(objs)

$(DIRDIR):
	@[ -d $@ ] || mkdir -p $@


$(OBJDIR)/libMns_$(CC)$(DEBUG).a:$(objs)
	@echo $@
	@\rm -f $@
	@ar crsu $@ $(objs)
	@ranlib $@
	@ln -fs $(DISTRIB)/$(PLATFORM)_$(CC)$(DEBUG)/libMns_$(CC)$(DEBUG).a  $(DISTRIB)/lib/libMns_$(CC)$(DEBUG).a

#clean:
#	-rm $(objs) $(EXEDIR)/$(prog)
#tar:$(DIRDIR)
#	tar czf $(ARCDIR)/mns.`date +"%Y.%m.%d"`.tgz sources makefile
#printf "install:C:$(MK_DEBUG_TOKEN):flex\t"`basename $<`"\n" 
#	mv $*.yy.cpp $@
#$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
#	$(CPP) $(CPPFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 
#target: $(prog)


