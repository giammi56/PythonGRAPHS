@echo off
IF EXIST folder_list.txt del /q folder_list.txt
dir /ad /b /s *ipch *Debug* *Release* *Win32* .vs logs > folder_list.txt
For /f "delims=" %%i in (folder_list.txt) do rmdir /s /q "%%i" 
del /q /S folder_list.txt

del /S *.vs.db
del /S *.suo /A:H
del /S *.ilk
del /S LMF2root.exe
del /S LMF2root.lib
del /S *.ncb
del /S *.sdf
del /S *.user
del /S *.aps
del /S *.pdb
del /S *.manifest
del /S *.manifest.res
del /S *manifest.rc
del /S *.lastbuildstate
del /S *.exp
del /S *.tlog
del /S *.VS.db
del /S *.log
del /S *.bsc
del /S *.sbr
del /S *.root
del /S *.bak
del /S *.db
rmdir /s /q "LMF2root.tlog" 
rmdir /s /q ".vs" 
rmdir /s /q work
rmdir /s /q logs
