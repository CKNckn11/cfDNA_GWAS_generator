#!/usr/bin/env python3
"""Compatibility wrapper for the installable command.

Prefer:
    cfdna-gwas-generate --config config.yaml
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from cfdna_gwas_generator.generator import main


if __name__ == "__main__":
    main()
