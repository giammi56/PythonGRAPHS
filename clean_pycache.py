import os
import glob


pyCacheFiles = glob.glob("ProjectLocation\**\*.pyc", recursive = True)
pyCacheFolders = glob.glob("ProjectLocation\**\__pycache__", recursive = True)


for file in pyCacheFiles:
    try:
        print(file)
        os.remove(file)
    except:
        print("Error removing .pyc files.")

for dir in pyCacheFolders:
    try:
        print(dir)
        os.removedirs(dir)
    except:
        print("Error removing _pycache__ directories.")