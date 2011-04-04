.DEFAULT_GOAL = help

.PHONY : dirs clean volume_tools ezminc minc oobicpl n3 bicpl ebtks conglomerate mincmorph mincblob pve ray_trace list

# some important environment variables
TEMP   ?= ../../temp
PREFIX ?= ../../quarantine
# specify number of parallel tasks in building (number of cores you can use for compiling)
PARALLEL_BUILD ?= -j4 


#CC=gcc-4.1  

#get absolute paths
TEMP   := $(shell realpath $(TEMP) )
PREFIX := $(shell realpath $(PREFIX) )

dirs := ${TEMP}/src ${PREFIX}

WGET := wget -c 
CMAKE := $(PREFIX)/bin/cmake

#package versions
NETCDF_VER       := 4.0.1
HDF5_VER         := 1.8.6
MINC_VER         := 2.1.00
BICPL_VER        := 1.4.6
EBTKS_VER        := 1.6.4
OOBICPL_VER      := 0.4.4
N3_VER           := 1.12.0
CLASSIFY_VER     := 1.1.0
CONGLOMERATE_VER := 1.6.6
MNI_AUTOREG_VER  := 0.99.6
MNI_PERLLIB_VER  := 0.08
MINCBLOB_VER     := 1.2.1
MINCMORPH_VER    := 1.4
RAYTRACE_VER     := 1.0.3
GLIMIMAGE_VER    := 1.2
ITK_VER          := 3.20.0
FFTW_VER         := 3.2.2
CMAKE_VER        := 2.8.4
GSL_VER          := 1.14
GETOPT_TABULAR_VER := 0.3
EZMINC_VER       := r83
NETPBM_VER       := 10.35.74
VTK_VER          := 5.6.1
FLTK_VER         := 1.3.x-r7725
REGISTER_VER     := 1.4.0
DISPLAY_VER      := 1.5.0


archives := \
${TEMP}/src/netcdf-$(NETCDF_VER).tar.gz \
${TEMP}/src/hdf5-$(HDF5_VER).tar.gz \
${TEMP}/src/minc-$(MINC_VER).tar.gz \
${TEMP}/src/bicpl-$(BICPL_VER).tar.gz \
${TEMP}/src/ebtks-$(EBTKS_VER).tar.gz \
${TEMP}/src/oobicpl-$(OOBICPL_VER).tar.gz \
${TEMP}/src/N3-$(N3_VER).tar.gz \
${TEMP}/src/classify-$(CLASSIFY_VER).tar.gz \
${TEMP}/src/conglomerate-$(CONGLOMERATE_VER).tar.gz \
${TEMP}/src/mni_autoreg-$(MNI_AUTOREG_VER).tar.gz \
${TEMP}/src/mni_perllib-$(MNI_PERLLIB_VER).tar.gz \
${TEMP}/src/mincblob-$(MINCBLOB_VER).tar.gz \
${TEMP}/src/mincmorph-$(MINCMORPH_VER).tar.gz \
${TEMP}/src/ray_trace-$(RAYTRACE_VER).tar.gz \
${TEMP}/src/glim_image-$(GLIMIMAGE_VER).tar.gz \
${TEMP}/src/InsightToolkit-$(ITK_VER).tar.gz \
${TEMP}/src/fftw-$(FFTW_VER).tar.gz \
${TEMP}/src/cmake-$(CMAKE_VER).tar.gz \
${TEMP}/src/gsl-$(GSL_VER).tar.gz \
${TEMP}/src/Getopt-Tabular-$(GETOPT_TABULAR_VER).tar.gz \
${TEMP}/src/ezminc-$(EZMINC_VER).tar.gz  \
${TEMP}/src/netpbm-$(NETPBM_VER).tar.gz \
${TEMP}/src/Register-$(REGISTER_VER).tar.gz \
${TEMP}/src/Display-$(DISPLAY_VER).tar.gz \

sources := $(archives:.tar.gz=)

#1. download the sources
${TEMP}/src/netcdf-$(NETCDF_VER).tar.gz : 
	$(WGET) ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-$(NETCDF_VER).tar.gz \
 -O ${TEMP}/src/netcdf-$(NETCDF_VER).tar.gz

${TEMP}/src/hdf5-$(HDF5_VER).tar.gz : 
	$(WGET) http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-$(HDF5_VER).tar.gz \
 -O ${TEMP}/src/hdf5-$(HDF5_VER).tar.gz

${TEMP}/src/minc-$(MINC_VER).tar.gz : 
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/minc-$(MINC_VER).tar.gz \
 -O ${TEMP}/src/minc-$(MINC_VER).tar.gz

${TEMP}/src/bicpl-$(BICPL_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/bicpl-$(BICPL_VER).tar.gz  \
  -O ${TEMP}/src/bicpl-$(BICPL_VER).tar.gz

${TEMP}/src/ebtks-$(EBTKS_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/ebtks-$(EBTKS_VER).tar.gz  \
 -O ${TEMP}/src/ebtks-$(EBTKS_VER).tar.gz

${TEMP}/src/oobicpl-$(OOBICPL_VER).tar.gz : 
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/oobicpl-$(OOBICPL_VER).tar.gz \
 -O ${TEMP}/src/oobicpl-$(OOBICPL_VER).tar.gz

${TEMP}/src/N3-$(N3_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/N3-$(N3_VER).tar.gz \
 -O ${TEMP}/src/N3-$(N3_VER).tar.gz

${TEMP}/src/classify-$(CLASSIFY_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/classify-$(CLASSIFY_VER).tar.gz \
 -O ${TEMP}/src/classify-$(CLASSIFY_VER).tar.gz

${TEMP}/src/conglomerate-$(CONGLOMERATE_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/conglomerate-$(CONGLOMERATE_VER).tar.gz \
 -O ${TEMP}/src/conglomerate-$(CONGLOMERATE_VER).tar.gz

${TEMP}/src/mni_autoreg-$(MNI_AUTOREG_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/mni_autoreg-$(MNI_AUTOREG_VER).tar.gz \
 -O ${TEMP}/src/mni_autoreg-$(MNI_AUTOREG_VER).tar.gz

${TEMP}/src/mni_perllib-$(MNI_PERLLIB_VER).tar.gz :
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/mni_perllib-$(MNI_PERLLIB_VER).tar.gz \
 -O ${TEMP}/src/mni_perllib-$(MNI_PERLLIB_VER).tar.gz

${TEMP}/src/mincblob-$(MINCBLOB_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/mincblob-$(MINCBLOB_VER).tar.gz \
 -O ${TEMP}/src/mincblob-$(MINCBLOB_VER).tar.gz

${TEMP}/src/mincmorph-$(MINCMORPH_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/mincmorph-$(MINCMORPH_VER).tar.gz \
 -O ${TEMP}/src/mincmorph-$(MINCMORPH_VER).tar.gz

${TEMP}/src/ray_trace-$(RAYTRACE_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/ray_trace-$(RAYTRACE_VER).tar.gz \
 -O ${TEMP}/src/ray_trace-$(RAYTRACE_VER).tar.gz

${TEMP}/src/glim_image-$(GLIMIMAGE_VER).tar.gz : 
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/glim_image-$(GLIMIMAGE_VER).tar.gz \
 -O ${TEMP}/src/glim_image-$(GLIMIMAGE_VER).tar.gz

${TEMP}/src/InsightToolkit-$(ITK_VER).tar.gz : 
	$(WGET) http://downloads.sourceforge.net/project/itk/itk/$(ITK_VER:.0=)/InsightToolkit-$(ITK_VER).tar.gz \
 -O ${TEMP}/src/InsightToolkit-$(ITK_VER).tar.gz

${TEMP}/src/vtk-$(VTK_VER).tar.gz : 
	$(WGET) http://www.vtk.org/files/release/5.6/vtk-$(VTK_VER).tar.gz \
 -O ${TEMP}/src/vtk-$(VTK_VER).tar.gz

${TEMP}/src/fltk-$(FLTK_VER).tar.bz2 : 
	$(WGET) http://ftp2.easysw.com/pub/fltk/snapshots/fltk-$(FLTK_VER).tar.bz2 \
 -O ${TEMP}/src/fltk-$(FLTK_VER).tar.bz2

${TEMP}/src/fftw-$(FFTW_VER).tar.gz :
	$(WGET)  http://www.fftw.org/fftw-$(FFTW_VER).tar.gz \
 -O ${TEMP}/src/fftw-$(FFTW_VER).tar.gz

${TEMP}/src/cmake-$(CMAKE_VER).tar.gz :
	$(WGET)  http://www.cmake.org/files/v2.8/cmake-$(CMAKE_VER).tar.gz \
 -O ${TEMP}/src/cmake-$(CMAKE_VER).tar.gz

${TEMP}/src/gsl-$(GSL_VER).tar.gz :
	$(WGET) http://mirrors.ibiblio.org/pub/mirrors/gnu/ftp/gnu/gsl/gsl-$(GSL_VER).tar.gz \
 -O ${TEMP}/src/gsl-$(GSL_VER).tar.gz

${TEMP}/src/Getopt-Tabular-$(GETOPT_TABULAR_VER).tar.gz : 
	$(WGET) http://search.cpan.org/CPAN/authors/id/G/GW/GWARD/Getopt-Tabular-$(GETOPT_TABULAR_VER).tar.gz \
 -O ${TEMP}/src/Getopt-Tabular-$(GETOPT_TABULAR_VER).tar.gz

${TEMP}/src/ezminc-$(EZMINC_VER).tar.gz : 
	$(WGET) http://www.bic.mni.mcgill.ca/~vfonov/software/ezminc-$(EZMINC_VER).tar.gz \
 -O ${TEMP}/src/ezminc-$(EZMINC_VER).tar.gz

${TEMP}/src/netpbm-$(NETPBM_VER).tar.gz : 
	$(WGET) 'http://downloads.sourceforge.net/project/netpbm/super_stable/$(NETPBM_VER)/netpbm-$(NETPBM_VER).tgz?use_mirror=cdnetworks-us-1' \
             -O ${TEMP}/src/netpbm-$(NETPBM_VER).tar.gz

${TEMP}/src/Display-$(DISPLAY_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/Display-$(DISPLAY_VER).tar.gz \
 -O ${TEMP}/src/Display-$(DISPLAY_VER).tar.gz

${TEMP}/src/Register-$(REGISTER_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/Register-$(REGISTER_VER).tar.gz \
 -O ${TEMP}/src/Register-$(REGISTER_VER).tar.gz



# unpacking archives:

${TEMP}/src/ezminc :  ${TEMP}/src/ezminc-$(EZMINC_VER).tar.gz  
	cd ${TEMP}/src && tar zxf  $<
	touch $@

${TEMP}/src/VTK : ${TEMP}/src/vtk-$(VTK_VER).tar.gz
	cd ${TEMP}/src && tar zxf  $<
	touch $@

${TEMP}/src/fltk-$(FLTK_VER) :  ${TEMP}/src/fltk-$(FLTK_VER).tar.bz2
	cd ${TEMP}/src && tar jxf  $<
	touch $@

$(sources) : % : %.tar.gz
	cd ${TEMP}/src && tar zxf  $<
	touch $@


netcdf         := $(PREFIX)/lib/libnetcdf.a $(PREFIX)/include/netcdf.h
hdf5           := $(PREFIX)/lib/libhdf5.a $(PREFIX)/include/hdf5.h
minc           := $(PREFIX)/lib/libminc2.a $(PREFIX)/include/minc2.h
fftw           := $(PREFIX)/lib/libfftw3f.a $(PREFIX)/include/fftw3.h
gsl            := $(PREFIX)/lib/libgsl.a $(PREFIX)/include/gsl/gsl_math.h
bicpl          := $(PREFIX)/lib/libbicpl.a $(PREFIX)/include/bicpl.h
itk            := $(PREFIX)/lib/InsightToolkit/UseITK.cmake
ezminc         := $(PREFIX)/lib/libminc_io.a $(PREFIX)/include/minc_1_rw.h $(PREFIX)/lib/libminc4itk.a
cmake          := $(PREFIX)/bin/cmake
bicpl          := $(PREFIX)/include/bicpl.h $(PREFIX)/lib/libbicpl.a
netpbm         := $(PREFIX)/lib/libnetpbm.a
n3             := $(PREFIX)/bin/nu_correct
ebtks          := $(PREFIX)/lib/libEBTKS.a
oobicpl        := $(PREFIX)/lib/liboobicpl.a
conglomerate   := $(PREFIX)/bin/mincdefrag
glim_image     := $(PREFIX)/bin/glim_image
mni_autoreg    := $(PREFIX)/bin/minctracc
mincblob       := $(PREFIX)/bin/mincblob
mincmorph      := $(PREFIX)/bin/mincmorph
ray_trace      := $(PREFIX)/bin/ray_trace
getopt_tabular := $(PREFIX)/perl/Getopt/Tabular.pm
mni_perllib    := $(PREFIX)/perl/MNI.pm
register       := $(PREFIX)/bin/register
Display        := $(PREFIX)/bin/Display
vtk            := $(PREFIX)/lib/vtk-5.6/UseVTK.cmake
fltk           := $(PREFIX)/lib/FLTK-1.3/UseFLTK.cmake


help : 
	@echo Usage: to compile all packages \"make all\"
	@echo        to compile core packages \"make core\"
	@echo        to compile visual programs \(Display, register\) \"make visual\"
	@echo        important variables:
	@echo           PREFIX - where to install packages
	@echo           TEMP  - where to compile everything
	@echo           PARALLEL_BUILD - how to execute parallel build \( -j4 \)
	@echo Prerequisits: realpath bison 


all : $(volume_tools) $(bicpl) $(n3) $(conglomerate) $(glim_image) $(mni_autoreg) $(mincblob) $(mni_perllib) $(ray_trace) $(ezminc) $(mincbet) $(pve) | dirs


core : $(dirs) minc itk ezminc 


clean : 
	rm -rf $(sources)

minc : $(dirs) $(minc) $(hdf5) $(netcdf)

ezminc : $(dirs)  $(ezminc)  

volume_tools : $(dirs)  $(volume_tools) 

itk : $(dirs)  $(itk)

vtk : $(dirs) $(vtk)

fltk : $(dirs) $(fltk)

bicpl : $(dirs) $(bicpl)

n3 : $(n3)

oobicpl : $(oobicpl)

conglomerate : $(conglomerate)

mincmorph : $(mincmorph)

mincblob : $(mincblob)

mincbet : $(mincbet)

pve : $(pve)

ray_trace : $(ray_trace)

getopt_tabular : $(getopt_tabular)

mni_perllib : $(mni_perllib)

register : $(dirs) $(register)

Display : $(dirs) $(Display)

cmake : $(dirs) $(cmake)

#3 configure & compile
$(netcdf) : ${TEMP}/src/netcdf-$(NETCDF_VER)
	cd ${TEMP}/src/netcdf-$(NETCDF_VER) && \
	./configure --prefix=$(PREFIX) --disable-netcdf4 --disable-hdf4 --disable-dap --disable-shared --disable-cxx --disable-f77 --disable-f90 --disable-examples --enable-v2 --disable-docs &&\
	make clean && \
	make && \
	make install

$(hdf5) : ${TEMP}/src/hdf5-$(HDF5_VER)
	cd ${TEMP}/src/hdf5-$(HDF5_VER) && \
	./configure --prefix=$(PREFIX) --disable-shared --disable-cxx --disable-fortran && \
	make clean && \
	make $(PARALLEL_BUILD) && \
	make install

$(fftw) : ${TEMP}/src/fftw-$(FFTW_VER)
	cd ${TEMP}/src/fftw-$(FFTW_VER) && \
	./configure --prefix=$(PREFIX) --disable-shared --enable-float --enable-threads  && \
	make clean  && \
	make $(PARALLEL_BUILD)  && \
	make install 

$(cmake) : ${TEMP}/src/cmake-$(CMAKE_VER)
	cd ${TEMP}/src/cmake-$(CMAKE_VER) && \
	./configure --prefix=$(PREFIX) && \
	make clean && \
	make && \
	make install && \
	touch $(PREFIX)/bin/cmake

$(gsl) : ${TEMP}/src/gsl-$(GSL_VER)
	cd ${TEMP}/src/gsl-$(GSL_VER) && \
	./configure --prefix=$(PREFIX) --disable-shared && \
	make clean && \
	make $(PARALLEL_BUILD) && \
	make install

$(minc) : ${TEMP}/src/minc-$(MINC_VER) $(netcdf) $(hdf5) $(getopt_tabular) $(mni_perllib)
	cd ${TEMP}/src/minc-$(MINC_VER) && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --disable-shared && make clean  && \
	make $(PARALLEL_BUILD) && \
	make install 

$(ebtks) : ${TEMP}/src/ebtks-$(EBTKS_VER) $(minc)
	cd ${TEMP}/src/ebtks-$(EBTKS_VER)  && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --disable-shared --with-minc2 && make clean  && \
	make $(PARALLEL_BUILD) && \
	make install

$(netpbm) : ${TEMP}/src/netpbm-$(NETPBM_VER)
	cd ${TEMP}/src/netpbm-$(NETPBM_VER)  && \
	wget http://www.bic.mni.mcgill.ca/~vfonov/software/netpbm.config -O Makefile.config && \
	make clean && \
	make && \
	rm -rf ${TEMP}/src/netpbm-$(NETPBM_VER)/pkg && \
	make package pkgdir=${TEMP}/src/netpbm-$(NETPBM_VER)/pkg 
	cp ${TEMP}/src/netpbm-$(NETPBM_VER)/pkg/link/libnetpbm.a $(PREFIX)/lib 
	cp ${TEMP}/src/netpbm-$(NETPBM_VER)/pkg/include/*.h $(PREFIX)/include

$(bicpl) : ${TEMP}/src/bicpl-$(BICPL_VER) $(netpbm) $(minc) 
	cd ${TEMP}/src/bicpl-$(BICPL_VER) && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --disable-shared --with-minc2 --with-image-netpbm  && \
	make clean   && \
	make $(PARALLEL_BUILD) && make install

$(oobicpl) : ${TEMP}/src/oobicpl-$(OOBICPL_VER) $(bicpl)
	cd ${TEMP}/src/oobicpl-$(OOBICPL_VER) && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --disable-shared --with-minc2 --with-image-netpbm  && \
	make clean   && \
	make $(PARALLEL_BUILD) && \
	make install

$(n3) : ${TEMP}/src/N3-$(N3_VER) $(ebtks) $(minc) 
	cd ${TEMP}/src/N3-$(N3_VER) && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --with-minc2 && make clean  && \
	make $(PARALLEL_BUILD) && \
	make install

$(conglomerate) : ${TEMP}/src/conglomerate-$(CONGLOMERATE_VER) 
	cd ${TEMP}/src/conglomerate-$(CONGLOMERATE_VER)  && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --with-minc2  && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install

$(glim_image) : ${TEMP}/src/glim_image-$(GLIMIMAGE_VER) $(minc)
	cd ${TEMP}/src/glim_image-$(GLIMIMAGE_VER) && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --with-minc2 && \
	make clean  && \
	make $(PARALLEL_BUILD) && \
	make install 

$(mni_autoreg) : ${TEMP}/src/mni_autoreg-$(MNI_AUTOREG_VER) $(minc)
	cd ${TEMP}/src/mni_autoreg-$(MNI_AUTOREG_VER) && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX) --with-minc2 && \
	make clean  && \
	make $(PARALLEL_BUILD) && \
	make install
	

$(mni_perllib) : ${TEMP}/src/mni_perllib-$(MNI_PERLLIB_VER) 
	cd ${TEMP}/src/mni_perllib-$(MNI_PERLLIB_VER) && \
	echo $(PREFIX)>yes && \
	perl ./Makefile.PL SITEPREFIX=$(PREFIX) INSTALLDIRS=site INSTALLSITELIB=$(PREFIX)/perl INSTALLSITEBIN=$(PREFIX)/bin \
		INSTALLSITEMAN1DIR=$(PREFIX) INSTALLSITEMAN3DIR=$(PREFIX) && < yes  && \
	make && \
	make install
	touch $(mni_perllib)
  
$(getopt_tabular) : ${TEMP}/src/Getopt-Tabular-$(GETOPT_TABULAR_VER)
	cd ${TEMP}/src/Getopt-Tabular-$(GETOPT_TABULAR_VER) && \
	perl ./Makefile.PL SITEPREFIX=$(PREFIX) INSTALLDIRS=site INSTALLSITELIB=$(PREFIX)/perl INSTALLSITEBIN=$(PREFIX)/bin \
		INSTALLSITEMAN1DIR=$(PREFIX) INSTALLSITEMAN3DIR=$(PREFIX) && \
	make && \
	make install
	touch $(getopt_tabular)

$(mincblob) : ${TEMP}/src/mincblob-$(MINCBLOB_VER) $(minc)
	cd ${TEMP}/src/mincblob-$(MINCBLOB_VER) && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --with-minc2 && \
	make clean  && \
	make $(PARALLEL_BUILD)  && \
	make install

$(mincmorph) : ${TEMP}/src/mincmorph-$(MINCMORPH_VER) $(minc)
	cd ${TEMP}/src/mincmorph-$(MINCMORPH_VER) && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --with-minc2  && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install

$(ray_trace) : ${TEMP}/src/ray_trace-$(RAYTRACE_VER) $(minc)
	cd ${TEMP}/src/ray_trace-$(RAYTRACE_VER) && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --with-minc2  && \
	make clean  && \
	make $(PARALLEL_BUILD)  && \
	make install

$(classify) : ${TEMP}/src/classify-$(CLASSIFY_VER) minc
	cd ${TEMP}/src/classify-$(CLASSIFY_VER) && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --with-minc2  && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install

$(register) : ${TEMP}/src/Register-$(REGISTER_VER) $(bicpl)
	cd ${TEMP}/src/Register-1.4.0 && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --with-minc2 --with-x && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 

$(Display) :${TEMP}/src/Display-$(DISPLAY_VER) $(bicpl)
	cd ${TEMP}/src/Display-1.5.0 && \
	./configure --prefix=$(PREFIX) --with-build-path=$(PREFIX)  --with-minc2 --with-x && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 
	

$(itk) : ${TEMP}/src/InsightToolkit-$(ITK_VER) $(cmake)
	cd  ${TEMP}/src && \
	mkdir -p build_itk && \
	cd build_itk && \
	$(CMAKE) \
			-D BUILD_EXAMPLES:BOOL=OFF \
			-D BUILD_SHARED_LIBS:BOOL=OFF \
			-D BUILD_TESTING:BOOL=OFF \
			-D CMAKE_BUILD_TYPE:STRING=Release \
			-D CMAKE_INSTALL_PREFIX:PATH=$(PREFIX) \
			-D ITK_USE_REVIEW:BOOL=ON \
				${TEMP}/src/InsightToolkit-3.20.0  && \
	make $(PARALLEL_BUILD) && \
	make install && \
	touch $(PREFIX)/lib/InsightToolkit/UseITK.cmake

$(fltk) : ${TEMP}/src/fltk-$(FLTK_VER) $(cmake)
	cd ${TEMP}/src/fltk-$(FLTK_VER) && \
	mkdir -p build_fltk && \
	cd build_fltk && \
	$(CMAKE) \
	-D BUILD_EXAMPLES:BOOL=OFF \
	-D BUILD_SHARED_LIBS:BOOL=OFF \
	-D BUILD_TESTING:BOOL=OFF \
	-D CMAKE_BUILD_TYPE:STRING=Release \
	-D CMAKE_INSTALL_PREFIX:PATH=$(PREFIX) \
	-D CMAKE_USE_PTHREADS:BOOL=1 \
	-D OPTION_USE_SYSTEM_LIBJPEG:BOOL=OFF \
	-D OPTION_USE_SYSTEM_LIBPNG:BOOL=OFF \
	-D OPTION_USE_SYSTEM_ZLIB:BOOL=OFF \
	-D USE_OPENGL:BOOL=ON \
	${TEMP}/src/fltk-1.3.x-r7725 && \
	make $(PARALLEL_BUILD) && \
	make install && \
	touch $(fltk)
	

$(vtk) : ${TEMP}/src/VTK
	cd ${TEMP}/src/VTK && \
	mkdir -p build_vtk && \
	cd build_vtk && \
	$(CMAKE) \
	-D BUILD_EXAMPLES:BOOL=OFF \
	-D BUILD_SHARED_LIBS:BOOL=OFF \
	-D BUILD_TESTING:BOOL=OFF \
	-D CMAKE_BUILD_TYPE:STRING=Release \
	-D CMAKE_INSTALL_PREFIX:PATH=$(PREFIX) \
	${TEMP}/src/VTK && \
	make $(PARALLEL_BUILD) && make install && \
	touch $(vtk)

$(ezminc) : ${TEMP}/src/ezminc $(minc) $(itk) $(cmake) $(fftw) $(gsl)
	cd ${TEMP}/src && \
	rm -rf   build_ezminc && \
	mkdir -p build_ezminc && \
	cd build_ezminc && \
	$(CMAKE)  \
      -D BUILD_EXAMPLES:BOOL=ON \
      -D BUILD_TESTS:BOOL=ON \
      -D CMAKE_BUILD_TYPE:STRING=Release \
      -D CMAKE_INSTALL_PREFIX:PATH=$(PREFIX) \
      -D ITK_DIR:PATH=$(PREFIX)/lib/InsightToolkit \
      -D MINC_ROOT:PATH=$(PREFIX) \
      -D USE_MINC2:BOOL=ON \
      -D FFTW3_ROOT:PATH=$(PREFIX) \
      -D HAVE_GSL:BOOL=ON \
      -D HAVE_FFTW3:BOOL=ON \
      -D GSL_ROOT:PATH=$(PREFIX) \
      -D BUILD_DISTORTION_CORRECTION:BOOL=ON \
      -D BUILD_MINCNLM:BOOL=ON \
      -D BUILD_TOOLS:BOOL=ON \
      ${TEMP}/src/ezminc && \
	make $(PARALLEL_BUILD) && \
	make install && \
	touch $(ezminc)


${TEMP}/src : 
	mkdir -p ${TEMP}/src
  
${PREFIX} :
	mkdir -p ${PREFIX}

dirs : 
	mkdir -p ${TEMP}/src
	mkdir -p ${PREFIX} ${PREFIX}/lib ${PREFIX}/include ${PREFIX}/bin


visual :  $(dirs) register Display


environment : ${PREFIX}/environment.sh ${PREFIX}/environment.csh


${PREFIX}/environment.sh : 
	echo -e export PATH=${PREFIX}/bin:'$$PATH' \\nexport PERL5LIB=${PREFIX}/perl:'$$PERL5LIB' >${PREFIX}/environment.sh

${PREFIX}/environment.csh : 
	echo -e setenv PATH ${PREFIX}/bin:'$${PATH}' \\nexport PERL5LIB=${PREFIX}/perl:'$${PERL5LIB}' >${PREFIX}/environment.csh

list : 
	@echo $(subst $(TEMP)/src/,,$(sources)) fltk-$(FLTK_VER) vtk-$(VTK_VER)
