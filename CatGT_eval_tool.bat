
:: You can call CatGT three ways:
::
:: 1) > CatGT cmd-line-parameters
:: 2) > runit.bat cmd-line-parameters
:: 3a) Edit parameters in runit.bat, then call it ...
:: 3b) > runit.bat
::
:: This script effectively says:
:: "If there are no parameters sent to runit.bat, call CatGT
:: with the parameters hard coded here, else, pass all of the
:: parameters through to CatGT."
::

@echo off
@setlocal enableextensions
@cd /d "%~dp0"

set LOCALARGS=-dir=D:\SC048_in\ -run=SC048_122920_ex -g=0 -t=0,0 ^
-prb_fld -t_miss_ok ^
-ap -prb=0 ^
-apfilter=butter,12,300,9000 -gfix=0.40,0.10,0.02 -gblcar ^
-dest=D:\SC048_out

if [%1]==[] (set ARGS=%LOCALARGS%) else (set ARGS=%*)

%~dp0CatGT %ARGS%

