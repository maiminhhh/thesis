#!/usr/bin/env bash
# Colon pipeline watchdog: runs on the server in screen `colon_watchdog`, independent of any Claude session.
# Every 2 min it appends one line per event to logs/watchdog_events.txt:
#   STAGE   new line in logs/pipeline_status.txt
#   JOB     a job finished (any status file under logs/04, logs/05*) - FAIL if exit != 0
#   SCREEN  a colon screen appeared or disappeared
#   STALL   no new status line or finished job for STALL_H hours while colon jobs are running
#   IDLE    no colon screens and no colon jobs left
# Claude sessions follow this file (tail) and report; it also is the record to catch up from after a restart.
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
EV=logs/watchdog_events.txt
STALL_H=${STALL_H:-6}
ev() { echo "$(date '+%F %T') $*" >> $EV; }
snap_jobs() { cat logs/04/status*.txt logs/05*/status*.txt 2>/dev/null | sort; }
screens() { screen -ls | grep -o '[0-9]*\.colon_[a-z0-9_]*' | sed 's/^[0-9]*\.//' | sort; }
n_status=$(wc -l < logs/pipeline_status.txt)
jobs_prev=$(snap_jobs); scr_prev=$(screens); last_progress=$(date +%s); stalled=0; idle=0
ev "WATCHDOG started (screens: $(echo $scr_prev | tr '\n' ' '))"
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
    echo "$new" | while read -r l; do
      [ -z "$l" ] && continue
      if echo "$l" | grep -q 'exit=0'; then ev "JOB ok   $(echo $l | tr -s ' ')"; else ev "JOB FAIL $(echo $l | tr -s ' ')"; fi
    done
    last_progress=$(date +%s); stalled=0
  fi
  jobs_prev=$jobs_now
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
