@echo off
echo ========================================
echo HYDRO SOLVER - Build Script
echo ========================================
echo.

echo Compiling hydro_solver.f90...
gfortran -o hydro_solver.exe hydro_solver.f90
if %errorlevel% neq 0 (
    echo ERROR: Compilation failed!
    pause
    exit /b 1
)

echo.
echo ========================================
echo Build successful!
echo Executable created: hydro_solver.exe
echo ========================================
echo.

echo Running hydro_solver.exe...
echo.
hydro_solver.exe
if %errorlevel% neq 0 (
    echo ERROR: Program execution failed!
    pause
    exit /b 1
)

echo.
echo ========================================
echo Execution completed!
echo Output files created:
echo   - solution.dat (CSV format with headers)
echo   - solution.csv (CSV format for plotting)
echo ========================================
echo.
pause

