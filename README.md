# DecoyPYrat.v2 - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectromtery Analyses

Improvement on DecoyPYrat script from Sanger (https://www.sanger.ac.uk/science/tools/decoypyrat)

## Preliminar

1. Download and execute the latest Python 3.* installation package from [here](https://www.python.org/downloads/windows/).  
_While either 32-bit (x86) or 64-bit (x86-64) versions should work just fine_
* Verify a successful installation by opening a command prompt window and navigating to your Python installation directory (default is `C:\Python27`).  Type `python` from this location to launch the Python interpreter.
    ```
    Microsoft Windows [Version 6.2.9200]
    (c) 2012 Microsoft Corporation. All rights reserved.
    
    C:\Users\Username>cd C:\Python27
    
    C:\Python27>python
    Python 2.7.8 (default, Jun 30 2014, 16:03:49) [MSC v.1500 32 bit (Intel)] on win
    32
    Type "help", "copyright", "credits" or "license" for more information.
    >>>
    ```

2. Add the environment variable "PYTHON3x_HOME" with the installed python.
*_In Windows 7 and Windows 8, simply searching for "environment variables" will present the option to `Edit the system environment variables`. This will open the `System Properties / Advanced` tab_  
*_In Windows XP, right click on `My Computer->Properties` to open `System Properties` and click on the `Advanced` tab._  
 1. On the `System Properties / Advanced` tab, click `Environment Variables` to open `User Variables` and `System Variables`
 2. Create a new `User Variable` named Variable name: `PYTHON3x_HOME` and  Variable value: `whatever your installation path was`
 ![](venv_win64/add_env.png)
![](https://camo.githubusercontent.com/767e3e7294af750e7db47ffb119cdc1154e2c79f/68747470733a2f2f662e636c6f75642e6769746875622e636f6d2f6173736574732f323939363230352f313035383236332f38643062376334632d313138352d313165332d383532622d3863653063303263623464322e706e67)
 3. Find the system variable called `Path` and click `Edit`  
![](https://camo.githubusercontent.com/da06b60252e8293d278d2027544d23602daa853b/68747470733a2f2f662e636c6f75642e6769746875622e636f6d2f6173736574732f323939363230352f313035383239342f30643734343936382d313138362d313165332d383766302d6531326166323330353030612e706e67)
 4. Add the following text to the end of the Variable value:  `;%PYTHON_HOME%\;%PYTHON_HOME%\Scripts\`
![](https://camo.githubusercontent.com/fb28d689631f2f4012741f6cf599dd52ed720b92/68747470733a2f2f662e636c6f75642e6769746875622e636f6d2f6173736574732f323939363230352f313035383237362f63333566353334612d313138352d313165332d386631622d6439343033633836643939662e706e67)


2. Add the environment variable "PYTHON3x_HOME"



Creating a new global system variable is quite simple and is one of those features hiding in plain sight. Please note the screenshots are for Windows Server 2008, however the process for most versions of Windows is almost identical with only a few of the screens different.

In the Control Panel, open the System option (alternately, you can right-click on My Computer and select Properties). Select the “Advanced system settings” link.

![](https://www.howtogeek.com/wp-content/uploads/sg/2010/11/image_thumb21.png)

In the System Properties dialog, click “Environment Variables”.

![](https://www.howtogeek.com/wp-content/uploads/sg/2010/11/image_thumb22.png)

In the Environment Variables dialog, click the New button underneath the “System variables” section.

![](https://www.howtogeek.com/wp-content/uploads/sg/2010/11/image_thumb23.png)

Enter the name of your new variable as well the value and click OK.

![](https://www.howtogeek.com/wp-content/uploads/sg/2010/11/image_thumb24.png)

You should now see your new variable listed under the “System variables” section. Click OK to apply the changes.

![](https://www.howtogeek.com/wp-content/uploads/sg/2010/11/image_thumb25.png)

You can now access your new system environment variable like you would any other. You can use it from the command line or batch scripts without having to define it.

![](https://www.howtogeek.com/wp-content/uploads/sg/2010/11/image_thumb26.png)


## Execution
1. Execute the decoyPYrat.v2.bat script.
2. Add the Target FASTA file (with *.fasta or *.fa extension)
2. Add the Decoy FASTA file (with *.fasta or *.fa extension)


