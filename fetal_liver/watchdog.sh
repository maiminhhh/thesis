#!/usr/bin/env bash
# Server-side watchdog for the fetal_liver reruns: runs in screen, independent of any Claude
# session. Every 2 min it checks the job status files and the driver, and appends events to
# logs/watchdog.log:
#   FAIL     a 05 job exited non-zero (main or best-annotator run)
#   STAGE    the driver logged a new stage ([pcr] lines in logs/rerun_downstream.log)
#   DIED     the driver's screen (fetal_pcr) is gone without 'all done'
#   HOURLY   progress summary (jobs done / failed / running, server load); PROGRESS every 6 h
# If logs/watchdog_email_to holds an address, FAIL / STAGE / DIED events are also emailed.
# If logs/watchdog_ntfy_topic holds an ntfy.sh topic, FAIL / STAGE / DIED events and a progress
# summary every 6 hours are pushed to it (phone notifications via the ntfy app).
# Run in screen:  screen -dmS fetal_watchdog bash watchdog.sh
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/fetal_liver
LOG=logs/watchdog.log
log() {
  echo "$(date '+%F %T') $1 $2" >> $LOG
  if [ "$1" != HOURLY ] && [ -s logs/watchdog_ntfy_topic ]; then
    local prio=default tag=information_source
    case $1 in FAIL|DIED) prio=high; tag=warning ;; STAGE) tag=white_check_mark ;; esac
    curl -s -m 20 -H "Title: fetal_liver $1" -H "Priority: $prio" -H "Tags: $tag" -d "$2" \
      "https://ntfy.sh/$(head -1 logs/watchdog_ntfy_topic)" > /dev/null 2>> $LOG
  fi
  if [ "$1" != HOURLY ] && [ -s logs/watchdog_email_to ]; then
    echo "$2" | mail -s "fetal_liver: $1 $(echo "$2" | cut -c1-60)" "$(head -1 logs/watchdog_email_to)" 2>> $LOG
  fi
}
count() { if [ -f "$1" ]; then grep -c "$2" "$1"; else echo 0; fi; }
running() {
  local n=0
  for p in $(pgrep -f 'python.*05\.'); do [ "$(readlink /proc/$p/cwd 2>/dev/null)" == "$PWD" ] && n=$((n + 1)); done
  echo $n
}
seen_fail=$(cat logs/05/status.txt logs/05_best_annotator/status.txt logs/05_best_annotator_early/status.txt logs/05_best_annotator_scvi/status.txt 2>/dev/null | grep -c 'exit=[1-9]')
seen_stage=$(grep -c '\[\(pcr\|scvi\)\]' logs/rerun_downstream.log)
last_hour=$(date +%H)
log START "watchdog started; $seen_fail failed jobs and $seen_stage driver lines already seen"
while true; do
  sleep 120
  fails=$(cat logs/05/status.txt logs/05_best_annotator/status.txt logs/05_best_annotator_early/status.txt logs/05_best_annotator_scvi/status.txt 2>/dev/null | grep 'exit=[1-9]')
  n_fail=$(echo -n "$fails" | grep -c 'exit=')
  if [ "$n_fail" -gt "$seen_fail" ]; then
    echo "$fails" | tail -n $((n_fail - seen_fail)) | while read -r l; do log FAIL "$l"; done
    seen_fail=$n_fail
  fi
  n_stage=$(grep -c '\[\(pcr\|scvi\)\]' logs/rerun_downstream.log)
  if [ "$n_stage" -gt "$seen_stage" ]; then
    grep '\[\(pcr\|scvi\)\]' logs/rerun_downstream.log | tail -n $((n_stage - seen_stage)) | while read -r l; do log STAGE "$l"; done
    seen_stage=$n_stage
  fi
  if grep -q '\[scvi\] all done' logs/rerun_downstream.log; then log STAGE "pipeline finished; watchdog exiting"; exit 0; fi
  if ! screen -ls | grep -q '\.fetal_\(pcr\|scvi_chain\)'; then log DIED "driver screens fetal_pcr and fetal_scvi_chain are gone without 'all done' (see logs/rerun_downstream.log)"; exit 1; fi
  if [ "$(date +%H)" != "$last_hour" ]; then
    last_hour=$(date +%H)
    [ $((10#$last_hour % 6)) -eq 0 ] && log PROGRESS "05 main: $(count logs/05/status.txt 'exit=0') ok, $(count logs/05/status.txt 'exit=[1-9]') failed of 156; best-annotator: $(count logs/05_best_annotator_early/status.txt 'exit=0') ok, $(count logs/05_best_annotator_early/status.txt 'exit=[1-9]') failed of 156; running $(running); load $(cut -d' ' -f3 /proc/loadavg)"
    log HOURLY "05 main: $(count logs/05/status.txt 'exit=0') ok, $(count logs/05/status.txt 'exit=[1-9]') failed of 156; best-annotator: $(count logs/05_best_annotator_early/status.txt 'exit=0') ok, $(count logs/05_best_annotator_early/status.txt 'exit=[1-9]') failed of 156; running $(running); load $(cut -d' ' -f3 /proc/loadavg)"
  fi
done
