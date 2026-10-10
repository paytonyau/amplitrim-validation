@echo off
rem Overnight validation of AmpliTrim: every test the sequence data on hand can support (see validation\overnight.py).
rem (python validation\check_build.py --plan / --run does the same from amplitrim_paths.json, and also checks the machine and fetches new data.)
rem Run it from the prompt where the amplitrim environment is active and R and vsearch are on PATH.
rem
rem   run_overnight.bat list    shows which steps are ready and what is missing (runs nothing)
rem   run_overnight.bat quick   about 30 minutes: small slices of every step, to check that every path and tool works
rem   run_overnight.bat         the real run, about 9 to 10 hours; steps already finished are not repeated if you start it again
rem
rem Keep the computer awake for the night (Windows sleep stops the run): powercfg /change standby-timeout-ac 0
setlocal

rem ---- paths: edit the ones that differ on this machine ----
set DATA=C:\data
set OUT=%DATA%\overnight
set MOCKA_HIGH=%DATA%\results3\mockA_high\ready
set MOCKA_ALL=%DATA%\results2\mockA_ee26\ready
set RUNINFO=%DATA%\mockA\runs_info.tsv
set ASVTABLE=%DATA%\mockA\Mock_staggered-A_ASVs.csv
set REFERENCE=%USERPROFILE%\Downloads\mockA_reference.fasta
set ITS=%DATA%\results4\cryobank_its3\ready
rem PRJEB42528 (V1-V2): the "ready" folder of its run, for example C:\data\results\prjeb42528\ready  (leave empty to skip its steps)
set V1V2=
rem FIGARO: the path of figaro.py, for example C:\figaro\figaro\figaro.py  (leave empty to skip the FIGARO steps)
set FIGARO=

set ARGS=--out "%OUT%" --mockA-high "%MOCKA_HIGH%" --mockA-all "%MOCKA_ALL%" --runinfo "%RUNINFO%" --asv-table "%ASVTABLE%" --reference "%REFERENCE%" --its "%ITS%"
if defined V1V2 set ARGS=%ARGS% --v1v2 "%V1V2%"
if defined FIGARO set ARGS=%ARGS% --figaro-script "%FIGARO%"

set MODE=
if /i "%1"=="quick" set MODE=--quick
if /i "%1"=="list" set MODE=--list
if /i "%1"=="quick" set ARGS=%ARGS% --out "%OUT%_quick"

python "%~dp0overnight.py" %ARGS% %MODE%
echo.
echo Finished. The summary is OVERNIGHT_SUMMARY.md in the output folder.
endlocal
