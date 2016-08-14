@echo off

@SET BASE=C:\MechatronixBuild
@SET LIBDIR=%BASE%\lib
@SET INCDIR=%BASE%\include

@if not exist C:\MechatronixBuild mkdir C:\MechatronixBuild
@if not exist %LIBDIR% mkdir %LIBDIR%
@if not exist %INCDIR% mkdir %INCDIR%
