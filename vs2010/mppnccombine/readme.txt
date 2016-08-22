Instructions for building mppnccombine on Windows with Visual Studio 2008

- add mppnccombine.c to the VC++ console project

- define DLL_NETCDF - only if linking dynamically to NetCDF!

- add $(ProjectDir)\..\shared to additional include directories - this contain a dummy (empty) unistd.h

- add $(ProjectDir)\..\netcdf\win32\include to additional include directories

- add $(ProjectDir)\..\netcdf\win32\lib to additional library directories

- add netcdfs.lib to additional libraries
