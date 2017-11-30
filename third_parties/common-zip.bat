SET DIR=%1
SET DIR1=%2
SET ZIPFILE=%DIR1%.zip

@if EXIST %DIR% (
  @echo.
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "%DIR% already expanded"
  @echo.
) else (
  @PowerShell -NonInteractive -Command "Expand-Archive -Path %ZIPFILE% -DestinationPath ."
  @PowerShell -Command "move-item -path %DIR1% -destination %DIR%"
	)
