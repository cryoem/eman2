set "BASH_EXE=C:\Program Files\Git\bin\bash.exe"

"%BASH_EXE%" -c "rsync -avzh --stats /c/Users/EMAN/workspace/win-installers/eman2.win.exe %DEPLOY_DEST%/eman2.win.unstable.exe"
