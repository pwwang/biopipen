"""Provide default variables

- BIOPIPEN_DIR: the root directory of the biopipen source
- REPORT_DIR: the root directory of the report
- SCRIPTS_DIR: the root directory of the scripts
"""
from pathlib import Path

BIOPIPEN_DIR = Path(__file__).parent.parent.resolve()
REPORT_DIR = BIOPIPEN_DIR / "reports"
SCRIPT_DIR = BIOPIPEN_DIR / "scripts"
