SET YEAR=%1
SET BITS=%2
SET DIR=%3

@IF %YEAR% == 2010  (
  @set VSCMAKE=Visual Studio 10 2010
) ELSE IF %YEAR% == 2012 (
  @set VSCMAKE=Visual Studio 11 2012
) ELSE IF %YEAR% == 2013 (
  @set VSCMAKE=Visual Studio 12 2013
) ELSE IF %YEAR% == 2015 (
  @set VSCMAKE=Visual Studio 14 2015
) ELSE IF %YEAR% == 2017 (
  @set VSCMAKE=Visual Studio 15 2017
) ELSE IF %YEAR% == 2019 (
  @set VSCMAKE=Visual Studio 16 2019
) ELSE (
  @echo.
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "Unsupported Visual Studio %YEAR%"
  @echo.
  @EXIT
)

@IF "%BITS%" == "x64" (@set VSCMAKE=%VSCMAKE% Win64)

@SET VSDIR=vs%YEAR%_%BITS%
@RMDIR /S /Q %VSDIR%
@mkdir %VSDIR%

@echo.
powershell -command write-host -foreground "green" -background "black" -nonewline "compiling %DIR% in %VSDIR%"
powershell -command write-host -foreground "green" -background "black" -nonewline "cmake -G "%VSCMAKE%" -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..\%DIR%"
@echo.
@cd %VSDIR%
@cmake -G "%VSCMAKE%" -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..\%DIR%
@cmake --build . --config Release  --target install

@cmake -G "%VSCMAKE%" -DCMAKE_INSTALL_PREFIX:PATH=..\lib_debug ..\%DIR%
@cmake --build . --config Debug  --target install
@cd ..

@RMDIR /S /Q %VSDIR%
