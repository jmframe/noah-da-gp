@echo off

copy ostIn_parallel.txt ostIn.txt

REM replace C:\bin with location of OSTRICH installation
mpiexec -n 4 C:\bin\OstrichMPI.exe

pause
