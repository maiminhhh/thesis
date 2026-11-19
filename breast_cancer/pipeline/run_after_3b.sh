#!/usr/bin/env bash
# Starts 03 -> 04 -> 05 -> figures (run_pipeline.sh) once both species' section-3b step has finished cleanly.
set -u
cd "$(dirname "$0")/.."
for sp in human mouse; do
  until grep -qE "3b FINISHED|STOP" logs/02/${sp}_3b.log 2>/dev/null; do sleep 120; done
  if grep -q STOP logs/02/${sp}_3b.log; then echo "STOP: $sp 3b failed, pipeline not started ($(date))"; exit 1; fi
done
echo "both species done, starting run_pipeline.sh ($(date))"
SKIP_WAIT=1 bash run_pipeline.sh > logs/pipeline.log 2>&1
