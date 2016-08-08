@echo off
SET BASE=C:\MechatronixBuild
SET LIBDIR=%BASE%\lib
SET INCDIR=%BASE%\include

if not exist C:\MechatronixBuild mkdir C:\MechatronixBuild
if not exist %LIBDIR% mkdir %LIBDIR%
if not exist %INCDIR% mkdir %INCDIR%

SET LIBNAME GenericContainer

copy vs2013_32\Debug\%LIBNAME%.lib   %LIBDIR%\%LIBNAME%_vs2013_x86_debug.lib
copy vs2013_32\Release\%LIBNAME%.lib %LIBDIR%\%LIBNAME%_vs2013_x86.lib
copy vs2013_64\Debug\%LIBNAME%.lib   %LIBDIR%\%LIBNAME%_vs2013_x64_debug.lib
copy vs2013_64\Release\%LIBNAME%.lib %LIBDIR%\%LIBNAME%_vs2013_x64.lib

xcopy /Y /I ..\src\*.hh               %INCDIR%
xcopy /Y /I ..\src_lua_interface\*.hh %INCDIR%
