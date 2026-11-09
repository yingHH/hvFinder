# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for the hvFinder Windows client.

Goal: be independent of the current working directory - paths are resolved
relative to this spec file
so the entry script and resources are found regardless of where pyinstaller is
invoked from.

Usage (from the repository root):
    pyinstaller --clean --noconfirm hvFinder.spec
Output:
    dist/hvFinder.exe
"""
import os
import sys

REPO = SPECPATH                      # directory of this spec (injected by PyInstaller)
ENTRY = os.path.join(REPO, "client", "gui.py")
ICON = os.path.join(REPO, "hvFinder.ico")

datas = []
for src, dst in ((os.path.join(REPO, "client", "settings.ini"), "."),
                 (os.path.join(REPO, "server", "configs"), "configs")):
    if os.path.exists(src):
        datas.append((src, dst))

hiddenimports = ["paramiko", "PySide6", "yaml", "cryptography"]

a = Analysis(
    [ENTRY],
    pathex=[REPO],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="hvFinder",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,               # PyInstaller skips UPX when it is not installed; the build still works
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    # .ico is a Windows resource; macOS wants .icns and Linux ignores icons,
    # so only Windows gets one (avoid a PyInstaller icon-format error).
    icon=[ICON] if (sys.platform == "win32" and os.path.exists(ICON)) else None,
)
