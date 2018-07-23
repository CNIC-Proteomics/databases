@ECHO off

:: get current path
REM SET DB_HOME=%~dp0
REM SET DB_HOME=%DB_HOME:~0,-1%
SET DB_HOME=%CD%

:: input parameters
REM :: Unit test
REM SET InFile=%DB_HOME%\test\test2.fa
REM SET OutFile=%DB_HOME%\test\test2.decoy.fa
SET /p InFile="Enter the target file (in FASTA format): "
SET /p OutFile="Enter the decoy file (in FASTA format): "

:: execute decoyPYrat
ECHO **
ECHO ** execute decoyPYrat
CMD /C " "%DB_HOME%/venv_win64/venv_win64_py3x/Scripts/activate.bat" && python "%DB_HOME%/src/decoyPYrat.v2.py" --output_fasta "%OutFile%" "%InFile%" "

SET /P DUMMY=Hit ENTER to continue...