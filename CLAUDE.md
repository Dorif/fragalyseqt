# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FragalyseQt is a desktop application for DNA fragment analysis (MLPA, QF-PCR, STR, TP-PCR) data processing from capillary electrophoresis instruments. It reads FSA and HID (ABIF) file formats and provides peak detection, sizing, and CSV export. Licensed under GNU AGPL v3.

## Running the Application

```bash
pip install -e .
python -m fragalyseqt.main
```

Or after installation, use the entry point:
```bash
fragalyseqt
```

There is no automated test suite, no linter configuration, and no CI/CD pipeline.

## Architecture

The application is a PyQt GUI app built on top of PyQtGraph (which abstracts PyQt5/PyQt6/PySide6). All Qt imports go through `pyqtgraph.Qt` to maintain compatibility across Qt bindings.

### Package Structure

The project uses a PEP 517 src layout:

```
pyproject.toml              — Build configuration (setuptools backend)
src/fragalyseqt/
├── __init__.py             — Package init with __version__
├── main.py                 — Entry point. Creates QApplication, instantiates FragalyseApp, runs event loop.
├── fragalyseqt.py          — Core module. Ui_MainWindow class with all UI construction and application logic.
├── localize.py             — UI string localization (English, Russian, Ukrainian, Romanian, French, Bulgarian).
├── sizestandards.py        — Dictionary of 70+ capillary electrophoresis size standards.
├── setvar.py               — Helper functions for analysis parameters and equipment detection.
├── fillarray.py            — Custom binary parser for HID file format.
└── boxes.py                — Cross-compatible Qt message box wrapper.
```

Files at the root: `README.md`, `LICENSE`, `COPYING`, `requirements.txt`, `.gitignore`, `CLAUDE.md`, `FragalyseQt.png`, `docs/`, `contrib/`.

### Key Design Characteristics

- **Global state**: `fragalyseqt.py` uses module-level global variables for application state (channel data, peak results, plot state, file data). Functions modify these globals directly.
- **Single-class UI**: All UI elements and logic live in `Ui_MainWindow.setupUi()` and its methods — there are no separate widget classes or MVC separation.
- **No .ui files**: The entire interface is constructed in Python code within `fragalyseqt.py`.
- **Qt abstraction**: All Qt imports use `pyqtgraph.Qt` (e.g., `from pyqtgraph.Qt.QtWidgets import ...`) rather than importing PyQt5/PyQt6/PySide6 directly. This is intentional for cross-binding compatibility and must be preserved.
- **Relative imports**: All intra-package imports use relative imports (e.g., `from .boxes import msgbox`).

### Dependencies

- **pyqtgraph** — Qt plotting and Qt binding abstraction layer
- **biopython** — Reading ABIF/FSA files (`Bio.SeqIO.AbiIO`)
- **scipy** — Peak detection (`scipy.signal.find_peaks`) and spline interpolation (`scipy.interpolate`)
- **charset-normalizer** — Handling non-Latin character encodings in file metadata
- **pybaselines** — Baseline correction and signal denoising

### Test Data

Test FSA/HID files from various instruments are in `docs/TEST_FILES/`. These come from NCBI OSIRIS, MLPAinter, NIST Forensic DNA dataset, and author-provided samples.
