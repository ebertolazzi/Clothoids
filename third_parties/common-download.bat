SET URL=%1
SET FILE=%2

@if EXIST %FILE% (
  @echo.
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "%FILE% already downloaded"
  @echo.
) else (
  PowerShell -NonInteractive -Command "Import-Module BitsTransfer; Start-BitsTransfer -Source \"%URL%\" -Destination %FILE%"
)
