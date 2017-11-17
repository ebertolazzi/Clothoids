@echo off

@rmdir vs2013_32 /s /q
@rmdir vs2013_64 /s /q
@mkdir vs2013_32

@pushd vs2013_32
@cmake -G "Visual Studio 12 2013" ..\..
@popd

@mkdir vs2013_64

@pushd vs2013_64
@cmake -G "Visual Studio 12 2013 Win64" ..\..
@popd

@cmake --build vs2013_32 --config Release
@cmake --build vs2013_64 --config Release
@cmake --build vs2013_32 --config Debug
@cmake --build vs2013_64 --config Debug
