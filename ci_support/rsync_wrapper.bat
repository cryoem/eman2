set "BASH_EXE=C:\Program Files\Git\bin\bash.exe"

"%BASH_EXE%" -c "rsync -avzh --stats /c/Users/EMAN/workspace/win-installers/eman2.win.exe %DEPLOY_DEST%/%1/eman2.win.%2.exe"
"%BASH_EXE%" -c "rsync -avzh --stats /c/Users/EMAN/workspace/win-installers/eman2_huge.win.exe %DEPLOY_DEST%/%1/eman2_huge.win.%2.exe"
