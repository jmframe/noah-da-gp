@echo OFF

echo q > unpause.txt

GOTO nextdemo

cd Demo1
md test\serial
copy * test\serial
cd test\serial
REM CALL OSTRICH.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 1 Serial is Complete"
PAUSE

cd Demo1
md test\parallel
copy * test\parallel
cd test\parallel
REM CALL mpirun.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 1 MPI Parallel is Complete"
PAUSE

cd Demo1
md test\fparallel
copy * test\fparallel
cd test\fparallel
CALL fmpirun.bat
cd ..\..\..

echo "Demo 1 FMPI Parallel is Complete"
PAUSE

cd Demo2
md test\serial
copy * test\serial
cd test\serial
REM CALL OSTRICH.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 2 Serial is Complete"
PAUSE

cd Demo2
md test\parallel
copy * test\parallel
cd test\parallel
REM CALL mpirun.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 2 MPI Parallel is Complete"
PAUSE

cd Demo2
md test\fparallel
copy * test\fparallel
cd test\fparallel
CALL fmpirun.bat
cd ..\..\..

echo "Demo 2 FMPI Parallel is Complete"
PAUSE

cd Demo3
md test\serial
copy * test\serial
cd test\serial
REM CALL OSTRICH.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 3 Serial is Complete"
PAUSE

cd Demo3
md test\parallel
copy * test\parallel
cd test\parallel
REM CALL mpirun.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 3 MPI Parallel is Complete"
PAUSE

cd Demo3
md test\fparallel
copy * test\fparallel
cd test\fparallel
CALL fmpirun.bat
cd ..\..\..

echo "Demo 3 FMPI Parallel is Complete"
PAUSE

cd Demo4
md test\serial
copy * test\serial
cd test\serial
REM CALL OSTRICH.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 4 Serial is Complete"
PAUSE

cd Demo4
md test\parallel
copy * test\parallel
cd test\parallel
REM CALL mpirun.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 4 MPI Parallel is Complete"
PAUSE

cd Demo4
md test\fparallel
copy * test\fparallel
cd test\fparallel
CALL fmpirun.bat
cd ..\..\..

echo "Demo 4 FMPI Parallel is Complete"
PAUSE

cd Demo5
md test\serial
copy * test\serial
cd test\serial
REM CALL OSTRICH.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 5 Serial is Complete"
PAUSE

cd Demo5
md test\parallel
copy * test\parallel
cd test\parallel
REM CALL mpirun.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 5 MPI Parallel is Complete"
PAUSE

cd Demo5
md test\fparallel
copy * test\fparallel
cd test\fparallel
CALL fmpirun.bat
cd ..\..\..

echo "Demo 5 FMPI Parallel is Complete"
PAUSE

cd Demo6
md test\serial
copy * test\serial
cd test\serial
REM CALL OSTRICH.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 6 Serial is Complete"
PAUSE

cd Demo6
md test\parallel
copy * test\parallel
cd test\parallel
REM CALL mpirun.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 6 MPI Parallel is Complete"
PAUSE

cd Demo6
md test\fparallel
copy * test\fparallel
cd test\fparallel
CALL fmpirun.bat
cd ..\..\..

echo "Demo 6 FMPI Parallel is Complete"
PAUSE

cd Demo7
md test\serial
copy * test\serial
cd test\serial
REM CALL OSTRICH.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 7 Serial is Complete"
PAUSE

cd Demo7
md test\parallel
copy * test\parallel
cd test\parallel
REM CALL mpirun.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 7 MPI Parallel is Complete"
PAUSE

cd Demo7
md test\fparallel
copy * test\fparallel
cd test\fparallel
CALL fmpirun.bat
cd ..\..\..

echo "Demo 7 FMPI Parallel is Complete"
PAUSE

cd Demo8
md test\serial
copy * test\serial
cd test\serial
REM CALL OSTRICH.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 8 Serial is Complete"
PAUSE

cd Demo8
md test\parallel
copy * test\parallel
cd test\parallel
REM CALL mpirun.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 8 MPI Parallel is Complete"
PAUSE

:nextdemo

cd Demo8
md test\fparallel
copy * test\fparallel
cd test\fparallel
CALL fmpirun.bat
cd ..\..\..

echo "Demo 8 FMPI Parallel is Complete"
PAUSE

cd Demo9
md test\serial
copy * test\serial
cd test\serial
REM CALL OSTRICH.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 9 Serial is Complete"
PAUSE

cd Demo9
md test\parallel
copy * test\parallel
cd test\parallel
REM CALL mpirun.bat < ..\..\..\unpause.txt
cd ..\..\..

echo "Demo 9 MPI Parallel is Complete"
PAUSE

cd Demo9
md test\fparallel
copy * test\fparallel
cd test\fparallel
CALL fmpirun.bat
cd ..\..\..

echo "Demo 9 FMPI Parallel is Complete"
PAUSE

PAUSE
