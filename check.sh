#!/bin/sh
# The regression net. Under a second, no file and no network access.
set -e
PYTHONPATH=. PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest tests -q "$@"
