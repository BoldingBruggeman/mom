This folder contains include files (.h, .inc, .mod) and static libraries (.lib) for NetCDF 3.6.3,
compiled from source with MS VisualStudio 2008 and Intel Visual Fortran 11.

The libraries include both the C and F90 interfaces.

The calling convention for external procedures differs from the Intel Visual Fortran default,
in that it is set to lower-case symbol names. This is needed for MOM4, which refers (in C) to
F90 symbols by lower-case names.

Preprocessor definitions for NetCDF build:

C library:

WIN32;NDEBUG;_WINDOWS;_USRDLL;VERSION=3.6.3;VISUAL_CPLUSPLUS;HAVE_STRERROR

WIN32;NDEBUG;_WINDOWS;_USRDLL;NETCDF_EXPORTS;VERSION=3.6.1-beta1;DLL_NETCDF;DLL_EXPORT;VISUAL_CPLUSPLUS;NC_DLL_EXPORT;HAVE_STRERROR
