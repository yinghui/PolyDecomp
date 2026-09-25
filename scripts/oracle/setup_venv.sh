#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Creates the Python environment for the GEOS oracle (shapely) in
# ${PD_ORACLE_VENV:-$HOME/.venvs/pd-oracle}.
set -euo pipefail
venv="${PD_ORACLE_VENV:-$HOME/.venvs/pd-oracle}"
python3 -m venv "${venv}"
"${venv}/bin/pip" install --quiet --upgrade pip
"${venv}/bin/pip" install --quiet -r "$(dirname "$0")/requirements.txt"
"${venv}/bin/python" -c "import shapely; print('shapely', shapely.__version__, 'GEOS', shapely.geos_version_string)"
