#!/usr/bin/env bash
# ClinicalASL entrypoint. Default: run pipeline on /input /output /workingdir [--inference-method ssvb|vaby] [--config cfg.json].
set -uo pipefail

case "${1:-}" in
  /bin/bash|bash|/bin/sh|sh|ls|cat|pwd|which|env)
    exec "$@" ;;
  python|python3|qasl|hd-bet|dcm2niix)
    exec "$@" ;;
esac

if [ $# -lt 3 ]; then
  echo "Usage: docker run ... <image> <inputdir> <outputdir> <workingdir> [--mode cvr|single] [--inference-method ssvb|vaby] [--config cfg.json]" >&2
  echo "   or: docker run -it ... <image> bash" >&2
  exit 1
fi

cd /opt/ClinicalASL/python
python -u run_pipeline.py "$@"
rc=$?
echo "ClinicalASL exit code: ${rc}"
exit ${rc}
