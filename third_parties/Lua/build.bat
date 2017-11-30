@IF [%1] EQU [] (SET YEAR=2013) else (SET YEAR=%1)
@IF [%2] EQU [] (SET BITS=x64)  else (SET BITS=%2)

@SET DIR=lua-5.3.4
@SET FILE=%DIR%.tar.gz
@SET URL=http://www.lua.org/ftp/$DIR.tar.gz
@IF EXIST lib_debug\lib\lua.lib (
  @echo.
  @powershell -command write-host -foreground "red" -background "yellow" -nonewline "lua already downloaded"
  @echo.
) ELSE (
  @CALL ..\common-download.bat %URL% %FILE%
  @CALL ..\common-tgz.bat lua %DIR%
)

@copy /Y CMakeLists.txt lua\CMakeLists.txt

@IF EXIST ..\..\lib3rd\lua_vs%YEAR%_%BITS%.lib (
  @echo.
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "lua already compiled"
  @echo.
) else (
  @CALL ..\common-cmake.bat %YEAR% %BITS% lua
)
