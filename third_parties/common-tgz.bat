SET DIR=%1
SET DIR1=%2
SET TARFILE=%DIR1%.tar
SET TGZFILE=%TARFILE%.gz

@if EXIST %DIR% (
  @echo.
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "%DIR% already expanded"
  @echo.
) else (
  @PowerShell -Command "if (-not (Get-Command Expand-7Zip -ErrorAction Ignore)) { Install-Package -Scope CurrentUser -Force 7Zip4PowerShell > $null } Expand-7Zip %TGZFILE% . ; Expand-7Zip %TARFILE% . ; Remove-Item %TARFILE%"
  @PowerShell -Command "move-item -path %DIR1% -destination %DIR%"
)
