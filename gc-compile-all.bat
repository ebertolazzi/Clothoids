@echo off
@pushd win_compile
@call gc-compile-vs2013
@call gc-compile-vs2015
@call gc-install-vs2013
@call gc-install-vs2015
@popd
