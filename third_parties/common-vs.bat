@IF [%1] EQU [] (SET YEAR=2013) else (SET YEAR=%1)
@IF [%2] EQU [] (SET BITS=x64)  else (SET BITS=%2)

@IF %YEAR% == 2010  (
  "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" %BITS%
) ELSE IF %YEAR% == 2012 (
  "C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\vcvarsall.bat" %BITS%
) ELSE IF %YEAR% == 2013 (
  "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\vcvarsall.bat" %BITS%
) ELSE IF %YEAR% == 2015 (
  "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" %BITS%
) ELSE IF %YEAR% == 2017 (
  "C:\Program Files (x86)\Microsoft Visual Studio 15.0\VC\vcvarsall.bat" %BITS%
) ELSE IF %YEAR% == 2019 (
  "C:\Program Files (x86)\Microsoft Visual Studio 16.0\VC\vcvarsall.bat" %BITS%
) ELSE (
  @echo.
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "Unsupported Visual Studio %YEAR%"
  @echo.
  @EXIT
)
