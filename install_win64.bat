@ECHO OFF

:: check env varibles are defined
IF "%PYTHON3x_HOME%"=="" GOTO :EndProcess1

:: get the python executable files
SET DB_HOME=%CD%
ECHO **
ECHO ** use the following paths for python3x:
ECHO %PYTHON3x_HOME%

:: go to home
CD %DB_HOME%

REM :: install the PIP packages
ECHO **
ECHO ** install the 'pip' package for python3x
CMD /C " "%PYTHON3x_HOME%/python" "%DB_HOME%/venv_win64/get-pip.py" "
ECHO " "%PYTHON3x_HOME%/python" "%DB_HOME%/venv_win64/get-pip.py" "

:: install virtualenv packages
ECHO **
ECHO ** install the 'virtualenv' packages for python3x
CMD /C " "%PYTHON3x_HOME%/Scripts/pip" install virtualenv"
CMD /C " "%PYTHON3x_HOME%/Scripts/pip" install virtualenvwrapper-win"

:: create virtual enviroment for the application in the local path
ECHO **
ECHO ** create virtualenv in python3x for the application
CMD /C " "%PYTHON3x_HOME%/Scripts/virtualenv" -p "%PYTHON3x_HOME%/python" "%DB_HOME%/venv_win64/venv_win64_py3x" "

:: active the virtualenv and install the required packages
ECHO **
ECHO ** active the virtualenv and install the required packages for each enviroment
CMD /C " "%DB_HOME%/venv_win64/venv_win64_py3x/Scripts/activate.bat" && pip install biopython"

GOTO :EndProcess

:EndProcess1
    ECHO PYTHON3x_HOME env. variable is NOT defined
    GOTO :EndProcess

:: wait to Enter => Good installation
:EndProcess
    SET /P DUMMY=End of installation. Hit ENTER to continue...
