#!/usr/bin/env bash
# Colon pipeline watchdog: runs on the server in screen `colon_watchdog`, independent of any Claude session.
# Every 2 min it appends one line per event to logs/watchdog_events.txt:
#   STAGE   new line in logs/pipeline_status.txt
#   JOB     a job finished (any status file under logs/04, logs/05*) - FAIL if exit != 0
#   SCREEN  a colon screen appeared or disappeared
#   STALL   no new status line or finished job for STALL_H hours while colon jobs are running
#   IDLE    no colon screens and no colon jobs left
# Claude sessions follow this file (tail) and report; it also is the record to catch up from after a restart.
# Phone: events also go to ntfy.sh (topic in logs/ntfy_topic): every STAGE / FAIL / STALL / IDLE / screen
# ending at once, and finished-ok jobs as a digest at most every DIGEST_H hours.
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
EV=logs/watchdog_events.txt
STALL_H=${STALL_H:-6}
DIGEST_H=${DIGEST_H:-3}
TOPIC=$(cat logs/ntfy_topic)
push() {   # push <priority> <tags> <message>
  curl -s --max-time 20 -o /dev/null -H "Title: colon pipeline" -H "Priority: $1" -H "Tags: $2" -d "$3" "https://ntfy.sh/$TOPIC" || true
}
ev() {
  echo "$(date '+%F %T') $*" >> $EV
  case "${1%% *}" in
    STAGE)  push default arrow_forward "$*" ;;
    STALL)  push high warning "$*" ;;
    IDLE)   push default checkered_flag "$*" ;;
    SCREEN) if [[ "$*" == *ended* ]]; then push default information_source "$*"; fi ;;
    JOB)    if [[ "$*" == "JOB FAIL"* ]]; then push high x "$*"; else ok_count=$((ok_count + 1)); fi ;;
  esac
}
ok_count=0; last_digest=$(date +%s)
snap_jobs() { cat logs/04/status*.txt logs/05*/status*.txt 2>/dev/null | sort; }
screens() { screen -ls | grep -o '[0-9]*\.colon_[a-z0-9_]*' | sed 's/^[0-9]*\.//' | sort; }
n_status=$(wc -l < logs/pipeline_status.txt)
jobs_prev=$(snap_jobs); scr_prev=$(screens); last_progress=$(date +%s); stalled=0; idle=0
ev "WATCHDOG started (screens: $(echo $scr_prev | tr '\n' ' '))"
push low eyes "Watchdog (re)started; watching: $(echo $scr_prev | tr '\n' ' ')"
while true; do
  sleep 120
  n=$(wc -l < logs/pipeline_status.txt)
  if [ "$n" -gt "$n_status" ]; then
    tail -n +$((n_status + 1)) logs/pipeline_status.txt | sed 's/^[0-9-]* [0-9:]* //' | while read -r l; do ev "STAGE $l"; done
    n_status=$n; last_progress=$(date +%s); stalled=0
  fi
  jobs_now=$(snap_jobs)
  new=$(comm -13 <(echo "$jobs_prev") <(echo "$jobs_now"))
  if [ -n "$new" ]; then
    while read -r l; do
      [ -z "$l" ] && continue
      if echo "$l" | grep -q 'exit=0'; then ev JOB ok "$(echo $l | tr -s ' ')"; else ev JOB FAIL "$(echo $l | tr -s ' ')"; fi
    done <<< "$new"
    last_progress=$(date +%s); stalled=0
  fi
  jobs_prev=$jobs_now
  if [ $ok_count -gt 0 ] && [ $(( $(date +%s) - last_digest )) -gt $(( DIGEST_H * 3600 )) ]; then
    done05=$(grep -c 'exit=0' logs/05/status.txt logs/05/status_raw_rerun.txt 2>/dev/null | awk -F: '{s+=$2} END {print s}')
    push low hourglass "$ok_count jobs finished OK in the last ${DIGEST_H}h (05 scored so far: $done05 / 114)"
    ok_count=0; last_digest=$(date +%s)
  fi
  scr_now=$(screens)
  for s in $(comm -23 <(echo "$scr_prev") <(echo "$scr_now")); do ev "SCREEN ended: $s"; done
  for s in $(comm -13 <(echo "$scr_prev") <(echo "$scr_now")); do ev "SCREEN started: $s"; done
  scr_prev=$scr_now
  running=0
  for p in $(pgrep -u "$USER" -f 'python|Rscript|exec/R'); do
    [ "$(readlink /proc/$p/cwd 2>/dev/null)" = "$PWD" ] && { running=1; break; }
  done
  if [ $running = 1 ] && [ $stalled = 0 ] && [ $(( $(date +%s) - last_progress )) -gt $(( STALL_H * 3600 )) ]; then
    ev "STALL no stage line or finished job for ${STALL_H}h while colon jobs are running"; stalled=1
  fi
  if [ $running = 0 ] && [ -z "$(echo "$scr_now" | grep -v watchdog)" ]; then
    [ $idle = 0 ] && ev "IDLE no colon screens or jobs left"; idle=1
  else idle=0; fi
done
