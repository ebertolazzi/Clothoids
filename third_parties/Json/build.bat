@SET URL=https://github.com/Tencent/rapidjson.git

@IF EXIST rapidjson (
  @echo.
  @powershell -command write-host -foreground "red" -background "yellow" -nonewline "rapidjson already downloaded"
  @echo.
) else (
  @git clone --depth 1 %URL%
)

@xcopy rapidjson/include ../../lib3rd /e /i /h
