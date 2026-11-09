# Building the desktop client (Windows / macOS / Linux)

## Short version

**PyInstaller cannot cross-compile.** Each platform needs its own build:

| Target | Script | Output |
|---|---|---|
| Windows | `build\build_windows.ps1` | `dist\hvFinder.exe` |
| macOS | `bash build/build_unix.sh` | `dist/hvFinder` (Mach-O) |
| Linux | `bash build/build_unix.sh` | `dist/hvFinder` (ELF) |

There is no `.app` bundle on macOS: the same one-file executable is produced.

The analysis pipeline itself runs on the **server** and needs none of this — the
executable only wraps `client/gui.py`, so you can also simply run
`python client/gui.py` on any platform with PySide6 installed.

```powershell
powershell -ExecutionPolicy Bypass -File build\build_windows.ps1
```

## When a rebuild is needed

| Change | Rebuild? |
|---|---|
| `client/gui.py`, `hvFinder.ico` | **yes** |
| `client/settings.ini` defaults | **yes** (onefile bundles a copy) |
| `server/configs/config2.yaml` | **yes** (the GUI reads it and uploads it to the server) |
| `server/*.py` (pipeline code) | **no** — it runs on the server, not inside the exe |
| `server/configs/species_aliases.csv` | no separate step (it is bundled with `server/configs`) |

## Option A — a Windows machine (recommended for Windows releases)

Requirements: Windows 10/11 and Miniconda/Anaconda with `conda` on `PATH`.

```powershell
cd <repo>
powershell -ExecutionPolicy Bypass -File build\build_windows.ps1
```

The script

1. checks the platform (exits with an error elsewhere);
2. creates the `hvFinder` conda environment from `build/environment-win.yml`, or
   reuses it when it already exists (Python 3.9);
3. checks the build dependencies (PySide6 / paramiko / PyYAML / cryptography /
   PyInstaller);
4. runs `pyinstaller --clean --noconfirm hvFinder.spec`;
5. verifies `dist/hvFinder.exe` and writes its size, SHA256 and the Python /
   PyInstaller / PySide6 versions to `build/BUILD_INFO.txt`.

The equivalent manual steps are:

```bat
conda env create -f build\environment-win.yml
conda run -n hvFinder python -m PyInstaller --clean --noconfirm hvFinder.spec
```

## Option B — macOS / Linux

Requirements: macOS or Linux, and Miniconda/Anaconda with `conda` on `PATH`
(the same setup as Option A — both platforms share the environment name
`hvFinder` and the same PyInstaller spec).

```bash
bash build/build_unix.sh                  # create/reuse the conda env `hvFinder`
bash build/build_unix.sh --refresh-env    # conda env update --prune, then build
bash build/build_unix.sh --env-name hvf-414   # use a differently named env
PYTHON=/opt/miniconda3/envs/hvFinder/bin/python bash build/build_unix.sh  # skip conda
```

The script

1. detects the platform (`Darwin` → macOS, `Linux` → Linux) and refuses to run
   elsewhere;
2. locates `conda` (`$CONDA_EXE`, then `PATH`), creates the environment from
   `build/environment-unix.yml` (Python 3.9, PySide6, paramiko, PyYAML,
   PyInstaller) or reuses it when it already exists — pass `--refresh-env` to
   bring an existing environment up to date;
3. verifies PySide6 / paramiko / PyYAML / PyInstaller inside that environment
   (and prints the ready-to-run `conda env update … --prune` command when
   something is missing);
4. runs `pyinstaller --noconfirm hvFinder.spec` with
   `--workpath build/build --distpath dist` (`--no-clean` keeps the cached build
   directory);
5. writes `dist/hvFinder`, copies `settings.ini` next to it, and records size,
   SHA256 and tool versions (including the conda environment name) in
   `build/BUILD_INFO.txt`.

The equivalent manual steps are:

```bash
conda env create -f build/environment-unix.yml
conda run -n hvFinder python -m PyInstaller --noconfirm hvFinder.spec
```

To freeze the exact environment that produced a release, export it with
`conda list -n hvFinder --explicit > build/conda-explicit-$(uname -s).txt`.

Notes:

- the binary is unsigned. On macOS Gatekeeper may refuse it; clear the quarantine
  attribute (`xattr -d com.apple.quarantine dist/hvFinder`) and codesign +
  notarise it before distributing.
- on Linux the Qt platform plugin needs the usual desktop libraries
  (`libGL`, `libEGL`, `xcb`); a headless server without an X/Wayland session
  cannot run the GUI even though the build succeeds.
- the icon is applied on Windows only — macOS wants `.icns` and Linux ignores
  icons for executables.

## Option C — GitHub Actions

`build/build-windows.yml` is a ready-made workflow. To enable it:

```bash
mkdir -p .github/workflows
cp build/build-windows.yml .github/workflows/build-windows.yml
```

Then trigger `build-windows-client` manually; the artifact is `hvFinder-windows`.
The workflow assumes the repository root contains the package directories
(`client/`, `server/`, …). If you place this project in a subdirectory, adjust
`environment-file` and `working-directory` accordingly.

## Why the spec looks the way it does

Every path is derived from PyInstaller's injected `SPECPATH`: the entry point is
`client/gui.py`, the bundled data are `client/settings.ini` and `server/configs`,
and `pathex=[REPO]`. Missing resources are skipped instead of aborting the build,
and the icon is optional, so the build works from any working directory.

## File encodings

| File | Encoding | Note |
|---|---|---|
| `build/build_windows.ps1` | **UTF-8 with BOM** | Windows PowerShell 5.1 reads scripts with the system code page unless a BOM is present, which corrupts non-ASCII text and can break parsing |
| `client/settings.ini` | UTF-8 (BOM tolerated) | read with `utf-8-sig` |
| `.py` / `.md` / `.sh` / `.yaml` | UTF-8, no BOM | project-wide convention |

If you edit `build_windows.ps1` in VS Code, make sure the encoding shown in the
status bar stays **"UTF-8 with BOM"**.

## Post-build checklist

- [ ] `dist\hvFinder.exe` starts when launched **from the `dist\` directory**
- [ ] the three tabs render and switch
- [ ] Tab 1 connects over SSH
- [ ] `settings.ini` is written/updated on exit
- [ ] the SHA256 in `build\BUILD_INFO.txt` matches the distributed file
- [ ] one sample runs end to end (Tab 2 then Tab 3)

## Notes

1. `settings.ini` is resolved next to the executable (next to the script when run
   from source); the build copies it into `dist\`, and a onefile build releases a
   bundled copy on first start if none is present.
2. `pip install -e .` also works and exposes `hvfinder-gui`, `pathogen_detect`,
   `mNGS_pip` and `get_script_base`; the build does not depend on it.
3. A onefile PySide6 build is roughly 50 MB and takes 3-8 s to start the first time
   (it unpacks into a temporary directory).
4. Some antivirus products flag PyInstaller onefile binaries heuristically; sign
   the executable before distributing it widely.
5. `upx=True` is set but PyInstaller simply skips compression when UPX is absent.
