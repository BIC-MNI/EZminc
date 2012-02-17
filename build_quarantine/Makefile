# =============================================================================
# Vlad's minc tools builder script
#
# PURPOSE: Build ALL or a subset of the MNI-BIC tools 
#
# DETAILS:
#     This script is supposed to download one or more packages, and the build em.
# Now since one may not always want to compile the kitchen sink, Vlad has
# conveniently provided an assortment of tagets which will permit you to
# build various subsets of the totality.  Here are some options:
#
# (a) make minc-full
#     -> this makes a bunch of stuff (but not ALL, ... not really)
#     -> for those of you who don't like reading Makefiles, "minc-full" 
#        includes: netcdf, hdf5, minc,
#                  bicpl, N3, conglomerate, glim_image, mni_autoreg,
#                  mincblob, mni_perllib, ray_trace, mincbet
#
# (b) make minc-only
#     -> just makes minc and dependencies (netcdf, hdf5)
#
# (c) make minc-with-itk
#     -> just like it sounds: minc + itk + ezminc
#      
# (d) make visual
#     -> make the visualization programs: register + Display + postf
#     -> I suspect that these have been segregated from the "all" build
#        since these programs are trickier to build ... i.e. more 
#        dependencies, therefore more to go wrong.
#
# (e) make minc-niak
#     -> all the tools necessary to use the Neuroimaging Analysis Kit 
#        http://www.nitrc.org/projects/niak
#     includes: create_output_dirs bicpl n3 conglomerate mni_autoreg 
#     mincblob mni_perllib ray_trace classify mincmorph inormalize ebtks 
#     arguments getopt_tabular 
#
#
# EXAMPLES:
#
# [A] >> make minc-only INSTALL_DIR=`realpath install_dir`  BUILD_DIR=`realpath build_dir`
#
# Notes --
# (1) we are building a relatively small install comprising minc (and dependencies) only 
# (2) INSTALL_DIR = the directory root to which you expect the built packages
#     to be installed. For example, this could be "/usr/local/bic" ... but if
#     it is, you had better be running as root. Default is the current dir.
# (3) we use the "realpath" program to expand out my local dir (eg. "install_dir")
#     into a full absolute path name. It would appear as if this Makefile does not
#     like relative paths. If on Debian/Ubuntu, type "sudo apt-get realpath".
# (4) BUILD_DIR is a temporary build directory. Feel free to wipe it out when 
#     the build is complete. Default is the current dir.
#
#
# [B] >>make minc-with-itk INSTALL_DIR=`realpath install_dir`  BUILD_DIR=`realpath build_dir`
#
# Notes -- 
# (1) the same as above, but add ezminc, itk, and their dependencies
#
#
# [C] >>make minc-full INSTALL_DIR=`realpath install_dir`  BUILD_DIR=`realpath build_dir` PARALLEL_BUILD=-j2
#     >>make visual INSTALL_DIR=`realpath install_dir`  BUILD_DIR=`realpath build_dir`
#     >>make models  INSTALL_DIR=`realpath install_dir`  BUILD_DIR=`realpath build_dir`
#
# Notes -- 
# (1) make all ... and use up to 2 cores
# (2) make register, Display and postf (needs the previous make to finish first)
#     NB: currently (20May2011) Display needs to have its final link line corrected in order
#         to produce an executable
# (3) install a whack of MNI-BIC models for use with mritotal, etc.
#
#
# DEPENDENCIES:
#  Ubuntu Lucid 10.04 LTS
#     (i) You do not want nor need netcdf or hdf5 from the repos, as we build  
#         these from source. If these are around, get rid of them if you can,
#         else linking could end of grabbing the wrong version
#    (ii) apt-get the following before you start:
#         -> build-essential (for compilers, make, etc)
#         -> realpath (see above)
#         -> zlib1g-dev (minc)
#         -> bison (minc)
#         -> flex (minc)
#         -> libx11-dev (ray-trace)
#         -> glutg3-dev (ray-trace)
#         -> libxmu-dev, libxi-dev (ray-trace)
#
#  sudo apt-get install build-essential realpath \
#                       zlib1g-dev bison flex \
#                       libx11-dev glutg3-dev libxmu-dev libxi-dev
#
#  
#  Ubuntu 11.04 
#  
#  sudo apt-get install build-essential realpath \
#                       zlib1g-dev bison flex \
#                       libx11-dev glutg3-dev libxmu-dev libxi-dev \
#                       automake libtool
#  
# =============================================================================
#

.DEFAULT_GOAL = help

.PHONY : create_output_dirs output_dirs clean list models \
         minc-only minc-full minc-with-itk ezminc  \
         bicpl ebtks conglomerate n3 mincmorph mincblob oobicpl ray_trace \
         mni-models_average305-lin mni-models_colin27-lin mni-models_icbm152-lin \
         mni-models_icbm152-nl-6gen mni-models_icbm152-nl-40gen mni_autoreg_model \
         help

help :
	@echo      This script is supposed to download one or more packages, and the build em.
	@echo  Now since one may not always want to compile the kitchen sink, Vlad has
	@echo  conveniently provided an assortment of tagets which will permit you to
	@echo  build various subsets of the totality.  Here are some options\:
	@echo 
	@echo  \(a\) make minc-full
	@echo      -\> this makes a bunch of stuff \(but not ALL, ... not really\)
	@echo      -\> for those of you who don\'t like reading Makefiles, \"minc-full\"
	@echo         includes: netcdf, hdf5, minc,
	@echo                   bicpl, N3, conglomerate, glim_image, mni_autoreg,
	@echo                   mincblob, mni_perllib, ray_trace, mincbet
	@echo 
	@echo  \(b\) make minc-only
	@echo      -\> just makes minc and dependencies \(netcdf, hdf5\)
	@echo 
	@echo  \(c\) make minc-with-itk
	@echo      -\> just like it sounds: minc + itk + ezminc
	@echo       
	@echo  \(d\) make visual
	@echo      -\> make the visualization programs: register, Display, and postf 
	@echo      -\> I suspect that these have been segregated from the \"all\" build
	@echo         since these programs are trickier to build ... i.e. more 
	@echo         dependencies, therefore more to go wrong.
	@echo 
	@echo EXAMPLE\:
	@echo make minc-only INSTALL_DIR=\`realpath install_dir\`  BUILD_DIR=\`realpath build_dir\` PARALLEL_BUILD=-j2



# some important environment variables
# default to current dir
BUILD_DIR?=$(shell pwd)
# default to current dir
INSTALL_DIR?=$(shell pwd)
# specify number of parallel tasks in building 
PARALLEL_BUILD?=-j1
#                                      (number of cores you can use for compiling)
WGET := wget -c                      # wget command
CMAKE := $(INSTALL_DIR)/bin/cmake    # where to find CMAKE (allow for over-ride)

# define output directories used by the build/install process
output_dirs := $(BUILD_DIR)/src \
			$(INSTALL_DIR)/include \
			$(INSTALL_DIR)/lib \
			$(INSTALL_DIR)/bin \
			$(INSTALL_DIR)/share/mni-models \
			$(INSTALL_DIR)/share/mni_autoreg \
			$(INSTALL_DIR)/share/man/man1 \
			$(INSTALL_DIR)/share/man/man3




# *****************************************************************************
# Package versions
# *****************************************************************************
#
NETCDF_VER       := 4.0.1
HDF5_VER         := 1.8.7
MINC_VER         := 9e69692
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
CMAKE_VER        := 2.8.6
GSL_VER          := 1.15
GETOPT_TABULAR_VER := 0.3
EZMINC_VER       := 24ce294
NETPBM_VER       := 10.35.74
VTK_VER          := 5.6.1
FLTK_VER         := 1.3.0
REGISTER_VER     := 1.4.0
DISPLAY_VER      := 1.5.0
POSTF_VER        := 1.0.03
INORMALIZE_VER   := 1.0.2
ARGUMENTS_VER    := 0.2.1
PCRE_VER         := 8.12
PCREPP_VER       := 0.9.5
MNI_MODELS_AVERAGE305_LIN_VER  := 1.1
MNI_MODELS_COLIN27_LIN_VER     := 1.1
MNI_MODELS_ICBM_LIN_VER        := 1.1
MNI_MODELS_ICBM_NL_6GEN_VER    := 1.0
MNI_MODELS_ICBM_NL_40GEN_VER   := 1.0
MNI_AUTOREG_MODEL_VER          := 1.1


# *****************************************************************************
# Name and path of compressed source packages following download 
# ... that is, where they are expected to be unless not yet downloaded
# *****************************************************************************
#
archives := \
$(BUILD_DIR)/src/netcdf-$(NETCDF_VER).tar.gz \
$(BUILD_DIR)/src/hdf5-$(HDF5_VER).tar.gz \
$(BUILD_DIR)/src/andrewjanke-minc-$(MINC_VER).tar.gz \
$(BUILD_DIR)/src/bicpl-$(BICPL_VER).tar.gz \
$(BUILD_DIR)/src/ebtks-$(EBTKS_VER).tar.gz \
$(BUILD_DIR)/src/oobicpl-$(OOBICPL_VER).tar.gz \
$(BUILD_DIR)/src/N3-$(N3_VER).tar.gz \
$(BUILD_DIR)/src/classify-$(CLASSIFY_VER).tar.gz \
$(BUILD_DIR)/src/conglomerate-$(CONGLOMERATE_VER).tar.gz \
$(BUILD_DIR)/src/mni_autoreg-$(MNI_AUTOREG_VER).tar.gz \
$(BUILD_DIR)/src/mni_perllib-$(MNI_PERLLIB_VER).tar.gz \
$(BUILD_DIR)/src/mincblob-$(MINCBLOB_VER).tar.gz \
$(BUILD_DIR)/src/mincmorph-$(MINCMORPH_VER).tar.gz \
$(BUILD_DIR)/src/ray_trace-$(RAYTRACE_VER).tar.gz \
$(BUILD_DIR)/src/glim_image-$(GLIMIMAGE_VER).tar.gz \
$(BUILD_DIR)/src/InsightToolkit-$(ITK_VER).tar.gz \
$(BUILD_DIR)/src/fftw-$(FFTW_VER).tar.gz \
$(BUILD_DIR)/src/cmake-$(CMAKE_VER).tar.gz \
$(BUILD_DIR)/src/gsl-$(GSL_VER).tar.gz \
$(BUILD_DIR)/src/Getopt-Tabular-$(GETOPT_TABULAR_VER).tar.gz \
$(BUILD_DIR)/src/vfonov-EZminc-$(EZMINC_VER).tar.gz  \
$(BUILD_DIR)/src/netpbm-$(NETPBM_VER).tar.gz \
$(BUILD_DIR)/src/Register-$(REGISTER_VER).tar.gz \
$(BUILD_DIR)/src/Display-$(DISPLAY_VER).tar.gz \
$(BUILD_DIR)/src/postf-$(POSTF_VER).tar.gz \
$(BUILD_DIR)/src/inormalize-${INORMALIZE_VER}.tar.gz \
$(BUILD_DIR)/src/arguments-$(ARGUMENTS_VER).tar.gz \
$(BUILD_DIR)/src/pcre-$(PCRE_VER).tar.gz \
$(BUILD_DIR)/src//pcre++-$(PCREPP_VER).tar.gz \
$(BUILD_DIR)/src/mni-models_average305-lin-$(MNI_MODELS_AVERAGE305_LIN_VER).tar.gz  \
$(BUILD_DIR)/src/mni-models_colin27-lin-$(MNI_MODELS_COLIN27_LIN_VER).tar.gz     \
$(BUILD_DIR)/src/mni-models_icbm152-lin-$(MNI_MODELS_ICBM_LIN_VER).tar.gz     \
$(BUILD_DIR)/src/mni-models_icbm152-nl-$(MNI_MODELS_ICBM_NL_6GEN_VER).tar.gz  \
$(BUILD_DIR)/src/mni-models_icbm152-nl-2009-$(MNI_MODELS_ICBM_NL_40GEN_VER).tar.gz  \
$(BUILD_DIR)/src/mni_autoreg_model-$(MNI_AUTOREG_MODEL_VER).tar.gz \



# *****************************************************************************
# Install target files
# ... these are the files that Make uses to check whether the target
#     needs to be rebuilt or not. 
# ... Make will look for these files and :
#     (1) if non-existent, will rebuild the package
#     (2) if existent, but out-of-date, will also rebuild  
# *****************************************************************************
#
netcdf         := $(INSTALL_DIR)/lib/libnetcdf.a $(INSTALL_DIR)/include/netcdf.h
hdf5           := $(INSTALL_DIR)/lib/libhdf5.a $(INSTALL_DIR)/include/hdf5.h
minc           := $(INSTALL_DIR)/lib/libminc2.a $(INSTALL_DIR)/include/minc2.h
fftw           := $(INSTALL_DIR)/lib/libfftw3f.a $(INSTALL_DIR)/include/fftw3.h
gsl            := $(INSTALL_DIR)/lib/libgsl.a $(INSTALL_DIR)/include/gsl/gsl_math.h
bicpl          := $(INSTALL_DIR)/lib/libbicpl.a $(INSTALL_DIR)/include/bicpl.h
itk            := $(INSTALL_DIR)/lib/InsightToolkit/UseITK.cmake
ezminc         := $(INSTALL_DIR)/lib/libminc_io.a $(INSTALL_DIR)/include/minc_1_rw.h $(INSTALL_DIR)/lib/libminc4itk.a
cmake          := $(INSTALL_DIR)/bin/cmake
bicpl          := $(INSTALL_DIR)/include/bicpl.h $(INSTALL_DIR)/lib/libbicpl.a
netpbm         := $(INSTALL_DIR)/lib/libnetpbm.a
n3             := $(INSTALL_DIR)/bin/nu_correct
ebtks          := $(INSTALL_DIR)/lib/libEBTKS.a
oobicpl        := $(INSTALL_DIR)/lib/liboobicpl.a
conglomerate   := $(INSTALL_DIR)/bin/mincdefrag
glim_image     := $(INSTALL_DIR)/bin/glim_image
mni_autoreg    := $(INSTALL_DIR)/bin/minctracc
mincblob       := $(INSTALL_DIR)/bin/mincblob
mincmorph      := $(INSTALL_DIR)/bin/mincmorph
ray_trace      := $(INSTALL_DIR)/bin/ray_trace
getopt_tabular := $(INSTALL_DIR)/perl/Getopt/Tabular.pm
mni_perllib    := $(INSTALL_DIR)/perl/MNI.pm
register       := $(INSTALL_DIR)/bin/register
Display        := $(INSTALL_DIR)/bin/Display
postf          := $(INSTALL_DIR)/bin/postf
vtk            := $(INSTALL_DIR)/lib/vtk-5.6/UseVTK.cmake
fltk           := $(INSTALL_DIR)/lib/FLTK-1.3/UseFLTK.cmake
classify       := $(INSTALL_DIR)/bin/classify
inormalize     := $(INSTALL_DIR)/bin/inormalize
arguments      := $(INSTALL_DIR)/lib/libarguments.a
pcre           := $(INSTALL_DIR)/lib/libpcre.a
pcrepp         := $(INSTALL_DIR)/lib/libpcre++.a
mni-models_average305-lin   := $(INSTALL_DIR)/share/mni-models/average305_t1_tal_lin.mnc.gz
mni-models_colin27-lin      := $(INSTALL_DIR)/share/mni-models/colin27_t1_tal_lin.mnc.gz
mni-models_icbm152-lin      := $(INSTALL_DIR)/share/mni-models/icbm_avg_152_t2_tal_lin.mnc.gz
mni-models_icbm152-nl-6gen  := $(INSTALL_DIR)/share/mni-models/icbm_avg_152_t1_tal_nlin_symmetric_VI.mnc.gz
mni-models_icbm152-nl-40gen := $(INSTALL_DIR)/share/mni-models/mni_icbm152_t1_tal_nlin_sym_09a.mnc.gz
mni_autoreg_model           := $(INSTALL_DIR)/share/mni_autoreg/average_305.mnc.gz





# *****************************************************************************
#
# Top-level targets
#
# *****************************************************************************
#
minc-full : create_output_dirs bicpl n3 conglomerate glim_image mni_autoreg mincblob mni_perllib ray_trace ezminc mincbet classify mincmorph inormalize oobicpl

minc-niak : create_output_dirs bicpl n3 conglomerate mni_autoreg mincblob mni_perllib ray_trace classify mincmorph inormalize ebtks arguments getopt_tabular

minc-with-itk : $(output_dirs) minc itk ezminc 

netcdf : $(output_dirs) $(netcdf) 
hdf5 : $(output_dirs) $(hdf5) 
minc : $(output_dirs) $(minc) $(hdf5) $(netcdf)
minc-only : $(output_dirs) $(minc) $(hdf5) $(netcdf)

models : $(output_dirs) mni-models_average305-lin \
                        mni-models_colin27-lin \
                        mni-models_icbm152-lin \
                        mni-models_icbm152-nl-6gen \
                        mni-models_icbm152-nl-40gen \
                        mni_autoreg_model

clean : 
	rm -rf $(source_dirs)

ezminc : $(output_dirs)  $(ezminc)  
itk : $(output_dirs)  $(itk)
vtk : $(output_dirs) $(vtk)
fltk : $(output_dirs) $(fltk)
bicpl : $(output_dirs) $(bicpl)
ebtks : $(output_dirs) $(ebtks)
n3 : $(output_dirs) $(n3)
oobicpl : $(output_dirs) $(oobicpl)
conglomerate : $(output_dirs)  $(conglomerate)
mincmorph : $(output_dirs) $(mincmorph)
mincblob : $(output_dirs) $(mincblob)
mincbet : $(output_dirs) $(mincbet)
ray_trace : $(output_dirs) $(ray_trace)
getopt_tabular : $(output_dirs) $(getopt_tabular)
mni_perllib : $(output_dirs) $(mni_perllib)
visual : $(output_dirs)  $(postf) $(register) $(Display)
register : $(output_dirs) $(register)
Display : $(output_dirs) $(Display)
postf : $(output_dirs) $(postf)
cmake : $(output_dirs) $(cmake)
classify : $(output_dirs) $(classify) 
fftw : $(output_dirs) $(fftw)
glim_image : $(output_dirs) $(glim_image)
mni_autoreg : $(output_dirs) $(mni_autoreg)
inormalize : $(output_dirs)  $(inormalize)
arguments : $(output_dirs) $(arguments) 
pcre      : $(output_dirs) $(pcre)
pcrepp    : $(output_dirs) $(pcrepp)
mni-models_average305-lin : $(mni-models_average305-lin)
mni-models_colin27-lin : $(mni-models_colin27-lin)
mni-models_icbm152-lin : $(mni-models_icbm152-lin)
mni-models_icbm152-nl-6gen : $(mni-models_icbm152-nl-6gen)
mni-models_icbm152-nl-40gen : $(mni-models_icbm152-nl-40gen)
mni_autoreg_model : $(mni_autoreg_model)


# *****************************************************************************
#
# Targets involved in the unpacking of downloaded tarballs
#
# *****************************************************************************
#
$(BUILD_DIR)/src/ezminc :  $(BUILD_DIR)/src/ezminc-$(EZMINC_VER).tar.gz  
	cd $(BUILD_DIR)/src && tar zxf  $<
	touch $@

$(BUILD_DIR)/src/VTK : $(BUILD_DIR)/src/vtk-$(VTK_VER).tar.gz
	cd $(BUILD_DIR)/src && tar zxf  $<
	touch $@

$(BUILD_DIR)/src/fltk-$(FLTK_VER) :  $(BUILD_DIR)/src/fltk-$(FLTK_VER)-source.tar.gz
	cd $(BUILD_DIR)/src && tar zxf  $<
	touch $@

# strip the .tar.gz suffix from the tarball to create the source directory name 
source_dirs := $(archives:.tar.gz=)
#
# here we see the use of a static pattern rule in which the target names
# ... are used to generate the dependency names (one for each target)
# ... In this case, we append ".tar.gz" to each source dir name ... which really
# ... should equal the value in $(archives), yes?
$(source_dirs) : % : %.tar.gz
	$(shell echo cd $(BUILD_DIR)/src && echo tar zxf  $< )
	cd $(BUILD_DIR)/src && tar zxf  $<
	touch $@



# *****************************************************************************
#
# Targets involved in downloading the compressed sources
#
# *****************************************************************************
#
$(BUILD_DIR)/src/netcdf-$(NETCDF_VER).tar.gz : 
	$(WGET) ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-$(NETCDF_VER).tar.gz \
 -O $(BUILD_DIR)/src/netcdf-$(NETCDF_VER).tar.gz

$(BUILD_DIR)/src/hdf5-$(HDF5_VER).tar.gz : 
	$(WGET) http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-$(HDF5_VER)/src/hdf5-$(HDF5_VER).tar.gz \
 -O $(BUILD_DIR)/src/hdf5-$(HDF5_VER).tar.gz

$(BUILD_DIR)/src/andrewjanke-minc-$(MINC_VER).tar.gz : 
	$(WGET) --no-check-certificate https://github.com/andrewjanke/minc/tarball/9e69692c15a1644457144d18c37d55c31f2fe107 \
	-O $(BUILD_DIR)/src/andrewjanke-minc-$(MINC_VER).tar.gz

$(BUILD_DIR)/src/bicpl-$(BICPL_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/bicpl-$(BICPL_VER).tar.gz  \
  -O $(BUILD_DIR)/src/bicpl-$(BICPL_VER).tar.gz

$(BUILD_DIR)/src/ebtks-$(EBTKS_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/ebtks-$(EBTKS_VER).tar.gz  \
 -O $(BUILD_DIR)/src/ebtks-$(EBTKS_VER).tar.gz

$(BUILD_DIR)/src/ebtks-1.6.4-VF.patch.gz : 
	$(WGET) http://www.bic.mni.mcgill.ca/~vfonov/software/ebtks-1.6.4-VF.patch.gz \
 -O $(BUILD_DIR)/src/ebtks-1.6.4-VF.patch.gz

$(BUILD_DIR)/src/oobicpl-$(OOBICPL_VER).tar.gz : 
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/oobicpl-$(OOBICPL_VER).tar.gz \
 -O $(BUILD_DIR)/src/oobicpl-$(OOBICPL_VER).tar.gz

$(BUILD_DIR)/src/N3-$(N3_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/N3-$(N3_VER).tar.gz \
 -O $(BUILD_DIR)/src/N3-$(N3_VER).tar.gz

$(BUILD_DIR)/src/classify-$(CLASSIFY_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/classify-$(CLASSIFY_VER).tar.gz \
 -O $(BUILD_DIR)/src/classify-$(CLASSIFY_VER).tar.gz

$(BUILD_DIR)/src/conglomerate-$(CONGLOMERATE_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/conglomerate-$(CONGLOMERATE_VER).tar.gz \
 -O $(BUILD_DIR)/src/conglomerate-$(CONGLOMERATE_VER).tar.gz

$(BUILD_DIR)/src/mni_autoreg-$(MNI_AUTOREG_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/mni_autoreg-$(MNI_AUTOREG_VER).tar.gz \
 -O $(BUILD_DIR)/src/mni_autoreg-$(MNI_AUTOREG_VER).tar.gz

$(BUILD_DIR)/src/mni_perllib-$(MNI_PERLLIB_VER).tar.gz :
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/mni_perllib-$(MNI_PERLLIB_VER).tar.gz \
 -O $(BUILD_DIR)/src/mni_perllib-$(MNI_PERLLIB_VER).tar.gz

$(BUILD_DIR)/src/mincblob-$(MINCBLOB_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/mincblob-$(MINCBLOB_VER).tar.gz \
 -O $(BUILD_DIR)/src/mincblob-$(MINCBLOB_VER).tar.gz

$(BUILD_DIR)/src/mincmorph-$(MINCMORPH_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/mincmorph-$(MINCMORPH_VER).tar.gz \
 -O $(BUILD_DIR)/src/mincmorph-$(MINCMORPH_VER).tar.gz

$(BUILD_DIR)/src/ray_trace-$(RAYTRACE_VER).tar.gz : 
	$(WGET)  http://packages.bic.mni.mcgill.ca/tgz/ray_trace-$(RAYTRACE_VER).tar.gz \
 -O $(BUILD_DIR)/src/ray_trace-$(RAYTRACE_VER).tar.gz

$(BUILD_DIR)/src/glim_image-$(GLIMIMAGE_VER).tar.gz : 
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/glim_image-$(GLIMIMAGE_VER).tar.gz \
 -O $(BUILD_DIR)/src/glim_image-$(GLIMIMAGE_VER).tar.gz

$(BUILD_DIR)/src/InsightToolkit-$(ITK_VER).tar.gz : 
	$(WGET) http://downloads.sourceforge.net/project/itk/itk/$(ITK_VER:.0=)/InsightToolkit-$(ITK_VER).tar.gz \
 -O $(BUILD_DIR)/src/InsightToolkit-$(ITK_VER).tar.gz

$(BUILD_DIR)/src/vtk-$(VTK_VER).tar.gz : 
	$(WGET) http://www.vtk.org/files/release/5.6/vtk-$(VTK_VER).tar.gz \
 -O $(BUILD_DIR)/src/vtk-$(VTK_VER).tar.gz

$(BUILD_DIR)/src/fltk-$(FLTK_VER)-source.tar.gz : 
	$(WGET) http://ftp2.easysw.com/pub/fltk/$(FLTK_VER)/fltk-$(FLTK_VER)-source.tar.gz \
 -O $(BUILD_DIR)/src/fltk-$(FLTK_VER)-source.tar.gz

$(BUILD_DIR)/src/fftw-$(FFTW_VER).tar.gz :
	$(WGET)  http://www.fftw.org/fftw-$(FFTW_VER).tar.gz \
 -O $(BUILD_DIR)/src/fftw-$(FFTW_VER).tar.gz

$(BUILD_DIR)/src/cmake-$(CMAKE_VER).tar.gz :
	$(WGET)  http://www.cmake.org/files/v2.8/cmake-$(CMAKE_VER).tar.gz \
 -O $(BUILD_DIR)/src/cmake-$(CMAKE_VER).tar.gz

$(BUILD_DIR)/src/gsl-$(GSL_VER).tar.gz :
	$(WGET) http://mirrors.ibiblio.org/pub/mirrors/gnu/ftp/gnu/gsl/gsl-$(GSL_VER).tar.gz \
 -O $(BUILD_DIR)/src/gsl-$(GSL_VER).tar.gz

$(BUILD_DIR)/src/Getopt-Tabular-$(GETOPT_TABULAR_VER).tar.gz : 
	$(WGET) http://search.cpan.org/CPAN/authors/id/G/GW/GWARD/Getopt-Tabular-$(GETOPT_TABULAR_VER).tar.gz \
 -O $(BUILD_DIR)/src/Getopt-Tabular-$(GETOPT_TABULAR_VER).tar.gz

$(BUILD_DIR)/src/arguments-$(ARGUMENTS_VER).tar.gz : 
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/arguments-$(ARGUMENTS_VER).tar.gz \
 -O $(BUILD_DIR)/src/arguments-$(ARGUMENTS_VER).tar.gz

$(BUILD_DIR)/src/pcre-$(PCRE_VER).tar.gz :
	$(WGET) http://downloads.sourceforge.net/project/pcre/pcre/$(PCRE_VER)/pcre-$(PCRE_VER).tar.gz \
 -O $(BUILD_DIR)/src/pcre-$(PCRE_VER).tar.gz

$(BUILD_DIR)/src//pcre++-$(PCREPP_VER).tar.gz :
	$(WGET) http://www.daemon.de/idisk/Apps/pcre++/pcre++-$(PCREPP_VER).tar.gz \
 -O $(BUILD_DIR)/src//pcre++-$(PCREPP_VER).tar.gz

$(BUILD_DIR)/src/vfonov-EZminc-$(EZMINC_VER).tar.gz : 
	$(WGET) --no-check-certificate https://github.com/vfonov/EZminc/tarball/$(EZMINC_VER) \
 -O $(BUILD_DIR)/src/vfonov-EZminc-$(EZMINC_VER).tar.gz

$(BUILD_DIR)/src/netpbm-$(NETPBM_VER).tar.gz : 
	$(WGET) 'http://downloads.sourceforge.net/project/netpbm/super_stable/$(NETPBM_VER)/netpbm-$(NETPBM_VER).tgz?use_mirror=cdnetworks-us-1' \
             -O $(BUILD_DIR)/src/netpbm-$(NETPBM_VER).tar.gz

$(BUILD_DIR)/src/Display-$(DISPLAY_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/Display-$(DISPLAY_VER).tar.gz \
 -O $(BUILD_DIR)/src/Display-$(DISPLAY_VER).tar.gz

$(BUILD_DIR)/src/Register-$(REGISTER_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/Register-$(REGISTER_VER).tar.gz \
 -O $(BUILD_DIR)/src/Register-$(REGISTER_VER).tar.gz

$(BUILD_DIR)/src/postf-$(POSTF_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/postf-$(POSTF_VER).tar.gz \
 -O $(BUILD_DIR)/src/postf-$(POSTF_VER).tar.gz
 
$(BUILD_DIR)/src/inormalize-${INORMALIZE_VER}.tar.gz : 
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/inormalize-$(INORMALIZE_VER).tar.gz \
 -O $(BUILD_DIR)/src/inormalize-$(INORMALIZE_VER).tar.gz

$(BUILD_DIR)/src/mni-models_average305-lin-$(MNI_MODELS_AVERAGE305_LIN_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/mni-models_average305-lin-$(MNI_MODELS_AVERAGE305_LIN_VER).tar.gz \
 -O $(BUILD_DIR)/src/mni-models_average305-lin-$(MNI_MODELS_AVERAGE305_LIN_VER).tar.gz

$(BUILD_DIR)/src/mni-models_colin27-lin-$(MNI_MODELS_COLIN27_LIN_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/mni-models_colin27-lin-$(MNI_MODELS_COLIN27_LIN_VER).tar.gz \
 -O $(BUILD_DIR)/src/mni-models_colin27-lin-$(MNI_MODELS_COLIN27_LIN_VER).tar.gz

$(BUILD_DIR)/src/mni-models_icbm152-lin-$(MNI_MODELS_ICBM_LIN_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/mni-models_icbm152-lin-$(MNI_MODELS_ICBM_LIN_VER).tar.gz \
 -O $(BUILD_DIR)/src/mni-models_icbm152-lin-$(MNI_MODELS_ICBM_LIN_VER).tar.gz

$(BUILD_DIR)/src/mni-models_icbm152-nl-$(MNI_MODELS_ICBM_NL_6GEN_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/mni-models_icbm152-nl-$(MNI_MODELS_ICBM_NL_6GEN_VER).tar.gz \
 -O $(BUILD_DIR)/src/mni-models_icbm152-nl-$(MNI_MODELS_ICBM_NL_6GEN_VER).tar.gz

$(BUILD_DIR)/src/mni-models_icbm152-nl-2009-$(MNI_MODELS_ICBM_NL_40GEN_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/mni-models_icbm152-nl-2009-$(MNI_MODELS_ICBM_NL_40GEN_VER).tar.gz \
 -O $(BUILD_DIR)/src/mni-models_icbm152-nl-2009-$(MNI_MODELS_ICBM_NL_40GEN_VER).tar.gz

$(BUILD_DIR)/src/mni_autoreg_model-$(MNI_AUTOREG_MODEL_VER).tar.gz :
	$(WGET) http://packages.bic.mni.mcgill.ca/tgz/mni_autoreg_model-$(MNI_AUTOREG_MODEL_VER).tar.gz \
 -O $(BUILD_DIR)/src/mni_autoreg_model-$(MNI_AUTOREG_MODEL_VER).tar.gz






# *****************************************************************************
# Targets for configure & compile
# ... the assumption is that source packages have already been downloaded 
#     and unpacked
# *****************************************************************************
#
$(netcdf) : $(BUILD_DIR)/src/netcdf-$(NETCDF_VER)
	cd $(BUILD_DIR)/src/netcdf-$(NETCDF_VER) && \
	./configure --prefix=$(INSTALL_DIR) \
	            --with-pic \
	            --disable-netcdf4 \
	            --disable-hdf4 \
	            --disable-dap \
	            --disable-shared \
	            --disable-cxx \
	            --disable-f77 \
	            --disable-f90 \
	            --disable-examples \
	            --enable-v2 \
	            --disable-docs &&\
	make clean && \
	make && \
	make install

$(hdf5) : $(BUILD_DIR)/src/hdf5-$(HDF5_VER)
	cd $(BUILD_DIR)/src/hdf5-$(HDF5_VER) && \
	./configure --prefix=$(INSTALL_DIR) \
	            --with-pic \
	            --disable-shared \
	            --disable-cxx \
	            --disable-fortran  \
	            --disable-hl && \
	make clean && \
	make $(PARALLEL_BUILD) && \
	make install

$(fftw) : $(BUILD_DIR)/src/fftw-$(FFTW_VER)
	cd $(BUILD_DIR)/src/fftw-$(FFTW_VER) && \
	./configure --prefix=$(INSTALL_DIR) --disable-shared --enable-float --enable-threads  && \
	make clean  && \
	make $(PARALLEL_BUILD)  && \
	make install 

# ====================================================================================================================
# existence of the package build directory is the dependency ...
# ... if it's not there, then go off and WGET and unzip/untar
$(cmake) : $(BUILD_DIR)/src/cmake-$(CMAKE_VER)
	$(warning In CMAKE recipe ... TGT is $(cmake) ... DEP is $(BUILD_DIR)/src/cmake-$(CMAKE_VER))
	cd $(BUILD_DIR)/src/cmake-$(CMAKE_VER) && \
	./configure --prefix=$(INSTALL_DIR) && \
	make clean && \
	make && \
	make install && \
	touch $(INSTALL_DIR)/bin/cmake

$(gsl) : $(BUILD_DIR)/src/gsl-$(GSL_VER)
	cd $(BUILD_DIR)/src/gsl-$(GSL_VER) && \
	./configure --prefix=$(INSTALL_DIR) --disable-shared && \
	make clean && \
	make $(PARALLEL_BUILD) && \
	make install

$(minc) : $(BUILD_DIR)/src/andrewjanke-minc-$(MINC_VER) $(netcdf) $(hdf5)
	cd $(BUILD_DIR)/src/andrewjanke-minc-$(MINC_VER) && \
	./autogen.sh && \
	./configure --prefix=$(INSTALL_DIR) \
	            --with-build-path=$(INSTALL_DIR)  \
	            --disable-shared \
	            --with-pic   && make clean  && \
	make $(PARALLEL_BUILD) && \
	make install 

$(ebtks) : $(BUILD_DIR)/src/ebtks-$(EBTKS_VER) $(BUILD_DIR)/src/ebtks-1.6.4-VF.patch.gz $(minc)
	cd $(BUILD_DIR)/src/ebtks-$(EBTKS_VER) && \
	zcat $(BUILD_DIR)/src/ebtks-1.6.4-VF.patch.gz|patch -p1 && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --disable-shared --with-minc2 && make clean  && \
	make $(PARALLEL_BUILD) && \
	make install

$(netpbm) : $(BUILD_DIR)/src/netpbm-$(NETPBM_VER)
	cd $(BUILD_DIR)/src/netpbm-$(NETPBM_VER)  && \
	wget http://www.bic.mni.mcgill.ca/~vfonov/software/netpbm.config -O Makefile.config && \
	make clean && \
	make && \
	rm -rf $(BUILD_DIR)/src/netpbm-$(NETPBM_VER)/pkg && \
	make package pkgdir=$(BUILD_DIR)/src/netpbm-$(NETPBM_VER)/pkg 
	cp $(BUILD_DIR)/src/netpbm-$(NETPBM_VER)/pkg/link/libnetpbm.a $(INSTALL_DIR)/lib 
	cp $(BUILD_DIR)/src/netpbm-$(NETPBM_VER)/pkg/include/*.h $(INSTALL_DIR)/include

$(bicpl) : $(BUILD_DIR)/src/bicpl-$(BICPL_VER) $(netpbm) $(minc) 
	cd $(BUILD_DIR)/src/bicpl-$(BICPL_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --disable-shared --with-minc2 --with-image-netpbm  && \
	make clean   && \
	make $(PARALLEL_BUILD) && make install

$(oobicpl) : $(BUILD_DIR)/src/oobicpl-$(OOBICPL_VER) $(bicpl) $(arguments) $(pcre) $(pcrepp)
	cd $(BUILD_DIR)/src/oobicpl-$(OOBICPL_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --disable-shared --with-minc2 --with-image-netpbm  && \
	make clean   && \
	make $(PARALLEL_BUILD) && \
	make install

$(n3) : $(BUILD_DIR)/src/N3-$(N3_VER) $(ebtks) $(minc) 
	cd $(BUILD_DIR)/src/N3-$(N3_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2 && make clean  && \
	make $(PARALLEL_BUILD) && \
	make install

$(conglomerate) : $(BUILD_DIR)/src/conglomerate-$(CONGLOMERATE_VER) 
	cd $(BUILD_DIR)/src/conglomerate-$(CONGLOMERATE_VER)  && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2  && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install

$(glim_image) : $(BUILD_DIR)/src/glim_image-$(GLIMIMAGE_VER) $(minc)
	cd $(BUILD_DIR)/src/glim_image-$(GLIMIMAGE_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2 && \
	make clean  && \
	make $(PARALLEL_BUILD) && \
	make install 

$(mni_autoreg) : $(BUILD_DIR)/src/mni_autoreg-$(MNI_AUTOREG_VER) $(minc)
	cd $(BUILD_DIR)/src/mni_autoreg-$(MNI_AUTOREG_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR) --with-minc2 && \
	make clean  && \
	make $(PARALLEL_BUILD) && \
	make install


$(mni_perllib) : $(BUILD_DIR)/src/mni_perllib-$(MNI_PERLLIB_VER) 
	$(warning In perllib recipe ...)
	ls -l
	cd $(BUILD_DIR)/src/mni_perllib-$(MNI_PERLLIB_VER) && \
	echo -e $(INSTALL_DIR) \\n\\n\\n > dummy_response_file && \
	echo perl ./Makefile.PL SITEPREFIX=$(INSTALL_DIR) \
	                   INSTALLDIRS=site \
	                   INSTALLSITELIB=$(INSTALL_DIR)/perl \
	                   INSTALLSITEBIN=$(INSTALL_DIR)/bin \
                           INSTALLSITEARCH=$(INSTALL_DIR)/perl \
                           INSTALLSITESCRIPT=$(INSTALL_DIR)/bin \
                           INSTALLARCHLIB=$(INSTALL_DIR)/perl \
	                   INSTALLSITEMAN1DIR=$(INSTALL_DIR)/share/man/man1 \
	                   INSTALLSITEMAN3DIR=$(INSTALL_DIR)/share/man/man3 < dummy_response_file  && \
	perl ./Makefile.PL SITEPREFIX=$(INSTALL_DIR) \
	                   INSTALLDIRS=site \
	                   INSTALLSITELIB=$(INSTALL_DIR)/perl \
	                   INSTALLSITEBIN=$(INSTALL_DIR)/bin \
                           INSTALLSITEARCH=$(INSTALL_DIR)/perl \
                           INSTALLSITESCRIPT=$(INSTALL_DIR)/bin \
                           INSTALLARCHLIB=$(INSTALL_DIR)/perl \
	                   INSTALLSITEMAN1DIR=$(INSTALL_DIR)/share/man/man1 \
	                   INSTALLSITEMAN3DIR=$(INSTALL_DIR)/share/man/man3 < dummy_response_file  && \
	make && \
	make install && \
	touch $(mni_perllib)

$(getopt_tabular) : $(BUILD_DIR)/src/Getopt-Tabular-$(GETOPT_TABULAR_VER)
	cd $(BUILD_DIR)/src/Getopt-Tabular-$(GETOPT_TABULAR_VER) && \
	perl ./Makefile.PL SITEPREFIX=$(INSTALL_DIR) INSTALLDIRS=site INSTALLSITELIB=$(INSTALL_DIR)/perl INSTALLSITEBIN=$(INSTALL_DIR)/bin \
		INSTALLSITEMAN1DIR=$(INSTALL_DIR) INSTALLSITEMAN3DIR=$(INSTALL_DIR) INSTALLSITEARCH=$(INSTALL_DIR)/perl INSTALLARCHLIB=$(INSTALL_DIR)/perl && \
	make && \
	make install
	touch $(getopt_tabular)

$(mincblob) : $(BUILD_DIR)/src/mincblob-$(MINCBLOB_VER) $(minc)
	cd $(BUILD_DIR)/src/mincblob-$(MINCBLOB_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2 && \
	make clean  && \
	make $(PARALLEL_BUILD)  && \
	make install

$(mincmorph) : $(BUILD_DIR)/src/mincmorph-$(MINCMORPH_VER) $(minc)
	cd $(BUILD_DIR)/src/mincmorph-$(MINCMORPH_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2  && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install

$(ray_trace) : $(BUILD_DIR)/src/ray_trace-$(RAYTRACE_VER) $(minc)
	cd $(BUILD_DIR)/src/ray_trace-$(RAYTRACE_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2  && \
	make clean  && \
	make $(PARALLEL_BUILD)  && \
	make install

$(classify) : $(BUILD_DIR)/src/classify-$(CLASSIFY_VER) minc
	cd $(BUILD_DIR)/src/classify-$(CLASSIFY_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2  && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install

$(register) : $(BUILD_DIR)/src/Register-$(REGISTER_VER) $(bicpl)
	cd $(BUILD_DIR)/src/Register-1.4.0 && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2 --with-x && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 

# build Display: add the LIBS argument so that we link against the netpbm lib (built earlier)
$(Display) :$(BUILD_DIR)/src/Display-$(DISPLAY_VER) $(bicpl)
	cd $(BUILD_DIR)/src/Display-1.5.0 && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2 --with-x LDFLAGS="-L$(INSTALL_DIR)/lib" LIBS="-lnetpbm" && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 

$(postf) : $(BUILD_DIR)/src/postf-$(POSTF_VER) $(bicpl)
	cd $(BUILD_DIR)/src/postf-$(POSTF_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2 --with-x && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 

$(inormalize) : $(BUILD_DIR)/src/inormalize-$(INORMALIZE_VER) $(ebtks)
	cd $(BUILD_DIR)/src/inormalize-$(INORMALIZE_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR)  --with-minc2 && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 


$(arguments) : $(BUILD_DIR)/src/arguments-$(ARGUMENTS_VER)
	cd $(BUILD_DIR)/src/arguments-$(ARGUMENTS_VER) && \
	./configure --prefix=$(INSTALL_DIR) --with-build-path=$(INSTALL_DIR) --disable-shared --disable-dependency-tracking && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 


$(pcre) : $(BUILD_DIR)/src/pcre-$(PCRE_VER)
	cd $(BUILD_DIR)/src/pcre-$(PCRE_VER) && \
	./configure --prefix=$(INSTALL_DIR) --disable-shared --disable-dependency-tracking && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 


$(pcrepp) : $(BUILD_DIR)/src//pcre++-$(PCREPP_VER)
	cd $(BUILD_DIR)/src//pcre++-$(PCREPP_VER) && \
	./configure --prefix=$(INSTALL_DIR) --disable-shared --disable-dependency-tracking  --with-pcre-include=$(INSTALL_DIR)/include --with-pcre-lib=$(INSTALL_DIR)/lib && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 



# build the MODELS
# ====================================================================================================================
# existence of the package build directory is the dependency ...
# ... if it's not there, then go off and WGET and unzip/untar
$(mni-models_average305-lin) : $(BUILD_DIR)/src/mni-models_average305-lin-$(MNI_MODELS_AVERAGE305_LIN_VER)
	$(warning In recipe for target $(mni-models_average305-lin) ...) 
	cd $(BUILD_DIR)/src/mni-models_average305-lin-$(MNI_MODELS_AVERAGE305_LIN_VER) && \
	./configure --prefix=$(INSTALL_DIR) && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 
	touch $@

$(mni-models_colin27-lin) :$(BUILD_DIR)/src/mni-models_colin27-lin-$(MNI_MODELS_COLIN27_LIN_VER)
	$(warning In recipe for target $(mni-models_colin27-lin) ...) 
	cd $(BUILD_DIR)/src/mni-models_colin27-lin-$(MNI_MODELS_COLIN27_LIN_VER) && \
	./configure --prefix=$(INSTALL_DIR) && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 
	touch $@

$(mni-models_icbm152-lin) :$(BUILD_DIR)/src/mni-models_icbm152-lin-$(MNI_MODELS_ICBM_LIN_VER)
	$(warning In recipe for target $(mni-models_icbm152-lin) ...) 
	cd $(BUILD_DIR)/src/mni-models_icbm152-lin-$(MNI_MODELS_ICBM_LIN_VER) && \
	./configure --prefix=$(INSTALL_DIR) && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 
	touch $@

$(mni-models_icbm152-nl-6gen) :$(BUILD_DIR)/src/mni-models_icbm152-nl-$(MNI_MODELS_ICBM_NL_6GEN_VER)
	$(warning In recipe for target $(mni-models_icbm152-nl-6gen) ...) 
	cd $(BUILD_DIR)/src/mni-models_icbm152-nl-$(MNI_MODELS_ICBM_NL_6GEN_VER) && \
	./configure --prefix=$(INSTALL_DIR) && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 
	touch $@

$(mni-models_icbm152-nl-40gen) :$(BUILD_DIR)/src/mni-models_icbm152-nl-2009-$(MNI_MODELS_ICBM_NL_40GEN_VER)
	$(warning In recipe for target $(mni-models_icbm152-nl-40gen) ...) 
	cd $(BUILD_DIR)/src/mni-models_icbm152-nl-2009-$(MNI_MODELS_ICBM_NL_40GEN_VER) && \
	./configure --prefix=$(INSTALL_DIR) && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 
	touch $@

$(mni_autoreg_model) :$(BUILD_DIR)/src/mni_autoreg_model-$(MNI_AUTOREG_MODEL_VER)
	cd $(BUILD_DIR)/src/mni_autoreg_model-$(MNI_AUTOREG_MODEL_VER) && \
	./configure --prefix=$(INSTALL_DIR) && \
	make clean   && \
	make $(PARALLEL_BUILD)  && \
	make install 
	touch $@


$(itk) : $(BUILD_DIR)/src/InsightToolkit-$(ITK_VER) $(cmake)
	cd  $(BUILD_DIR)/src && \
	mkdir -p build_itk && \
	cd build_itk && \
	$(CMAKE) \
			-D BUILD_EXAMPLES:BOOL=OFF \
			-D BUILD_SHARED_LIBS:BOOL=OFF \
			-D BUILD_TESTING:BOOL=OFF \
			-D CMAKE_BUILD_TYPE:STRING=Release \
			-D CMAKE_INSTALL_PREFIX:PATH=$(INSTALL_DIR) \
			-D ITK_USE_REVIEW:BOOL=ON \
				$(BUILD_DIR)/src/InsightToolkit-3.20.0  && \
	make $(PARALLEL_BUILD) && \
	make install && \
	touch $(INSTALL_DIR)/lib/InsightToolkit/UseITK.cmake

$(fltk) : $(BUILD_DIR)/src/fltk-$(FLTK_VER) $(cmake)
	cd $(BUILD_DIR)/src/fltk-$(FLTK_VER) && \
	mkdir -p build_fltk && \
	cd build_fltk && \
	$(CMAKE) \
	-D BUILD_EXAMPLES:BOOL=OFF \
	-D BUILD_SHARED_LIBS:BOOL=OFF \
	-D BUILD_TESTING:BOOL=OFF \
	-D CMAKE_BUILD_TYPE:STRING=Release \
	-D CMAKE_INSTALL_PREFIX:PATH=$(INSTALL_DIR) \
	-D CMAKE_USE_PTHREADS:BOOL=1 \
	-D OPTION_USE_SYSTEM_LIBJPEG:BOOL=OFF \
	-D OPTION_USE_SYSTEM_LIBPNG:BOOL=OFF \
	-D OPTION_USE_SYSTEM_ZLIB:BOOL=OFF \
	-D USE_OPENGL:BOOL=ON \
	$(BUILD_DIR)/src/fltk-$(FLTK_VER) && \
	make $(PARALLEL_BUILD) && \
	make install && \
	touch $(fltk)
	

$(vtk) : $(BUILD_DIR)/src/VTK
	cd $(BUILD_DIR)/src/VTK && \
	mkdir -p build_vtk && \
	cd build_vtk && \
	$(CMAKE) \
	-D BUILD_EXAMPLES:BOOL=OFF \
	-D BUILD_SHARED_LIBS:BOOL=OFF \
	-D BUILD_TESTING:BOOL=OFF \
	-D CMAKE_BUILD_TYPE:STRING=Release \
	-D CMAKE_INSTALL_PREFIX:PATH=$(INSTALL_DIR) \
	$(BUILD_DIR)/src/VTK && \
	make $(PARALLEL_BUILD) && make install && \
	touch $(vtk)

$(ezminc) : $(BUILD_DIR)/src/vfonov-EZminc-$(EZMINC_VER) $(minc) $(itk) $(cmake) $(fftw) $(gsl)
	cd $(BUILD_DIR)/src && \
	rm -rf   build_ezminc && \
	mkdir -p build_ezminc && \
	cd build_ezminc && \
	$(CMAKE)  \
      -D BUILD_EXAMPLES:BOOL=ON \
      -D BUILD_TESTS:BOOL=ON \
      -D CMAKE_BUILD_TYPE:STRING=Release \
      -D CMAKE_INSTALL_PREFIX:PATH=$(INSTALL_DIR) \
      -D ITK_DIR:PATH=$(INSTALL_DIR)/lib/InsightToolkit \
      -D MINC_ROOT:PATH=$(INSTALL_DIR) \
      -D USE_MINC2:BOOL=ON \
      -D FFTW3_ROOT:PATH=$(INSTALL_DIR) \
      -D HAVE_GSL:BOOL=ON \
      -D HAVE_FFTW3:BOOL=ON \
      -D GSL_ROOT:PATH=$(INSTALL_DIR) \
      -D BUILD_DISTORTION_CORRECTION:BOOL=ON \
      -D BUILD_MINCNLM:BOOL=ON \
      -D BUILD_TOOLS:BOOL=ON \
      $(BUILD_DIR)/src/vfonov-EZminc-$(EZMINC_VER) && \
	make $(PARALLEL_BUILD) && \
	make install && \
	touch $(ezminc)

# ensure that the output directories (build and install) exist
${output_dirs} : create_output_dirs 
  
create_output_dirs : 
	mkdir -p $(BUILD_DIR)/src
	mkdir -p $(INSTALL_DIR)/bin
	mkdir -p $(INSTALL_DIR)/lib
	mkdir -p $(INSTALL_DIR)/include
	mkdir -p $(INSTALL_DIR)/share/mni-models
	mkdir -p $(INSTALL_DIR)/share/mni_autoreg
	mkdir -p $(INSTALL_DIR)/share/man/man1
	mkdir -p $(INSTALL_DIR)/share/man/man3

# build the visual packages 
visual :  $(output_dirs) register Display


# permit 'make list'
# ... lists all of the package names buildable by this Makefile
list : 
	@echo $(subst $(BUILD_DIR)/src/,,$(source_dirs)) fltk-$(FLTK_VER) vtk-$(VTK_VER)


# append the intall dir to the PATH and PERL5LIB env variables
# ... does anyone actually call this target?
environment : $(INSTALL_DIR)/environment.sh $(INSTALL_DIR)/environment.csh

$(INSTALL_DIR)/environment.sh : 
	echo export PATH=$(INSTALL_DIR)/bin:'$$PATH' \\nexport PERL5LIB=$(INSTALL_DIR)/perl:'$$PERL5LIB' >$(INSTALL_DIR)/environment.sh

$(INSTALL_DIR)/environment.csh : 
	echo setenv PATH $(INSTALL_DIR)/bin:'$${PATH}' \\nsetenv PERL5LIB=$(INSTALL_DIR)/perl:'$${PERL5LIB}' >$(INSTALL_DIR)/environment.csh


# *****************************************************************************
# Custom functions
# *****************************************************************************

#define print_header
#	echo ***************************************************************
#	echo *
#	echo *    Building ... $1
#	echo *
#	echo *
#	echo *
#	echo ***************************************************************
#endef
