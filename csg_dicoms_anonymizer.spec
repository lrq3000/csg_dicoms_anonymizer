# -*- mode: python -*-

block_cipher = None

import os, sys
cur_path = os.path.realpath('.')
sys.path.append(os.path.join(cur_path, 'csg_dicoms_anonymizer', 'csg_fileutil_libs'))  # for gooey spec file, because it does not support relative paths (yet?)

import gooey
gooey_root = os.path.dirname(gooey.__file__)
gooey_languages = Tree(os.path.join(gooey_root, 'languages'), prefix = 'gooey/languages')
gooey_images = Tree(os.path.join(gooey_root, 'images'), prefix = 'gooey/images')

a = Analysis([os.path.join('csg_dicoms_anonymizer', 'csg_dicoms_anonymizer.py')],
             pathex=[os.path.join(cur_path, 'csg_dicoms_anonymizer')],
             binaries=[],
             datas=[],
             hiddenimports=[os.path.join(cur_path, 'csg_dicoms_anonymizer', 'csg_fileutil_libs')],
             hookspath=[],
             runtime_hooks=[],
             excludes=['pandas', 'numpy', 'matplotlib', 'mpl-data', 'zmq', 'IPython', 'ipykernel', 'tcl', 'Tkinter', 'jupyter_client', 'ipywidgets', 'unittest', 'ipython', 'ipython_genutils', 'jupyter_core'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          gooey_languages, # Add them in to collected files
          gooey_images, # Same here.
          name='csg_dicoms_anonymizer',
          debug=False,
          strip=False,
          upx=True,
          windowed=True,
          console=True )
