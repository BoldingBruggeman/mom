The document describes the MS Visual Studio 2010 configuration for MOM5.

For parallel computing support, Microsoft MPI is used.
See https://msdn.microsoft.com/en-us/library/bb524831(v=vs.85).aspx

For NetCDF support a custom VS2010-built version of NetCDF 3.6.3 is used, located at vs2010/netcdf.


==============================================
Basic things
==============================================

- define use_netCDF;use_netCDF3;use_libMPI

- use include directories ..\netcdf\win32\include;$(MSMPI_INC);..\..\src\mom4p1\ocean_param\gotm-4.0\include;..\..\src\shared\include;..\..\src\shared\mpp\include;..\..\src\shared\drifters;..\..\src\mom4p1\ocean_core (for 64-bit builds: replace netcdf\win32 by netcdf\x64)

- use library directories $(MSMPI_LIB32);$(NETCDF_PREFIX)\lib; (for 64-bit build: replace MSMPI_LIB32 with MSMPI_LIB64)

- use libraries netcdfs.lib msmpi.lib msmpifec.lib

==============================================
Special things for Windows compilation
==============================================

- define NO_DEV_NULL to prevent MOM4 of trying to open /dev/null for non-root output (in MPI mode). Output will be redirected to "._mpp.nonrootpe.msgs.NODEID" instead.

- switch ifort to lower case calling convention, because the C compiler will export all MOM4 C routines with lower case names (NB make sure the Fortran NetCDF library that is linked against uses the lower case convention as well!)

- due to the switch to the lower case convention, the symbols in the Intel libraries are no longer recognized and brought in during linking. This applies to the FLUSH and GETENV subroutines defined in the portability library. To fix this, add preprocessor defines GETENV=getenv_ifort and FLUSH=flush_ifort to the F90 project, along with the custom file ifort.F90.

- create empty file unistd.h, and reference its directory in the additional include directories of the C project. This UNIX/Linux header file does not exist on Windows. Substituing a dummy does not seem to have harmul effects.

- define preprocessor symbol _USE_MATH_DEFINES for the MOM4 C code. This ensures that M_PI is defined in math.h (M_PI is no longer a standard part of the C ISO standard, and is therefore not available by default).

- create integer function mld_id in threadloc.c, calling mld_id_ in turn. This is needed to resolve the external dependency mld_id in mpp_util.inc. This would not be needed if the Fortran compiler automatically appended underscores to external names (/assume:underscore in ifort), but I expect that this would break links to other external procedures [MPI].

- switch C project from dynamical linking in C library (DLL) to statical linking, because names in runtime libraries of MSVC and IFORT conflict.

- swap "if (fabs(det) < EPSLN15 ) return 0;" and "const long double deti = 1.0/det;" in shared/mosaic/mosaic_util.c to ensure all variables are declared at routine start (likely not an issue for GNU because they do not compile with -pedantic)

- insert return after "multiplier = 1.0" in mem_dump (src/shared/memutils/memutils.F90) and mpp_mem_dump (src/shared/memutils/mpp_memutils.F90) because these routines access /proc/self/status [not there on Windows]

- mppnccombine: commented out code related to POSIX-only sys/resource.h, moved min subroutine to precede main, moved declaration of several local variables to start of routine, set project property C/C++, Advanced, "Compile As" to C++ to enable inline keyword

==============================================
Issues while running
==============================================

If the number if PEs does not match that set in the input.nml configuration file:

in coupler_nml: update atmos_npes and ocean_npes
in ocean_model_nml: update layout


From the MOM4 FAQ:

Q: I run a new experiment and get this error:

FATAL: MPP_UPDATE_DOMAINS user stack overflow: call mpp_domains_set_stack_size( 70656) from all PEs.
 
A: The problem can be solved by adding namelist fms_nml in your runscripts.

    &fms_nml
       domains_stack_size = 70656/

where 70656 is the number appeared in the error message.
 
