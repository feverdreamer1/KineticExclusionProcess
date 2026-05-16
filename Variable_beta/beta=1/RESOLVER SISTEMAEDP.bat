@echo off
echo HYDRO SOLVER - Build Script
echo.

echo Compiling hydro_solver.f90...
gfortran -o hydro_solver.exe hydro_solver.f90
if %errorlevel% neq 0 (
    echo ERROR: Compilation failed!
    pause
    exit /b 1
)

hydro_solver.exe
if %errorlevel% neq 0 (
    echo ERROR: Program execution failed!
    pause
    exit /b 1
)

pause

