REM @pyinstaller csg_fileutil_dicoms_anonymizer.py
REM @pyinstaller --noconsole --onefile csg_fileutil_dicoms_anonymizer.spec
pyinstaller --onedir csg_dicoms_anonymizer.spec > pyinstaller-log.txt 2>&1 & type pyinstaller-log.txt
pyi-archive_viewer dist\csg_dicoms_anonymizer.exe -r -b > pyinstaller-dependencies.txt
pause
