#  Declaration
CXXFLAGS = -Wall -v -O3 -std=c++0x
PARFLAGS = -fopenmp
CC	=	gcc
C++	=	g++
F77    	=	f77

.SUFFIXES:

.SUFFIXES: .f .F .cpp .c .o .h

%.o: %.f
	${F77} ${FFLAGS} $< -c
	
%.o: %.F
	${F77} ${FFLAGS} $< -c

%.o: %.cpp %.h
	${C++} ${CXXFLAGS} $< -c

%.o: %.c %.h
	${CC} ${CFLAGS} $< -c


objects =  App_Content_Dist.o App_Network_3D.o Backbone_Network.o Collision_detection.o Contact_grid.o Cutoff_Wins.o Direct_Electrifying.o \
           Electrical_analysis.o Generate_Network.o Geometry_3D.o Hns.o Hoshen_Kopelman.o Input_Reader.o \
           MathMatrix.o Printer.o Shells.o Triangulation.o VTK_Export.o MainPro.o \
	                   
necn : $(objects)        
	${C++} ${PARFLAGS} -o necn $(objects)
	              
MainPro.o : MainPro.cpp
	${C++} ${CXXFLAGS} MainPro.cpp -c

# clean
.PHONY : clean
clean :
	-rm necn $(objects)
