#!/usr/bin/env bash
# breast_cancer watchdog: runs in screen (bc_watchdog), independent of any Claude session.
# Every 5 minutes (important events are also pushed to the phone via ntfy, see pipeline/ntfy_send.sh):
#   logs/STATUS.txt  - snapshot: screens, pipeline stage, per-dataset 04/05 progress, failures, running jobs
#   logs/EVENTS.log  - append-only, one timestamped line per new event (stage change, job done/failed, STOP,
#                      repair exit, choice ready, pipeline finished, possible stall = a running job whose log has
#                      not changed for STALL_HOURS)
# Start:  screen -dmS bc_watchdog bash pipeline/watchdog.sh
# Do not edit while it runs (bash reads scripts as it goes); stop it with `screen -S bc_watchdog -X quit` first.
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/breast_cancer
EV=logs/EVENTS.log; ST=logs/STATUS.txt; SEEN=logs/.watchdog_seen
STALL_HOURS=${STALL_HOURS:-3}
DATASETS="full_mouse dissected_mouse genes2 genes4 genes6 genes8"
touch $SEEN
# events that also go to the phone (ntfy); failures, stops and stalls at high priority
NOTIFY_PAT='FAILED|STOP|POSSIBLE STALL|\] == |exit [0-9]|CHOICE READY|PIPELINE FINISHED|SKIP:|LABEL_COLS='
ev() {
  echo "$(date '+%F %T') $*" >> $EV
  if echo "$*" | grep -qE "$NOTIFY_PAT"; then
    prio=3; echo "$*" | grep -qE 'FAILED|STOP|STALL|exit [1-9]' && prio=5
    bash pipeline/ntfy_send.sh "breast_cancer" "$*" $prio > /dev/null 2>&1
  fi
}
ev "watchdog started (pid $$)"

while true; do
  # ---- events: lines worth knowing about in the pipeline / repair logs, each reported once ----
  for f in logs/pipeline.log logs/04/run_04.log logs/sensitivity_choice.log logs/after_genes6.log \
           logs/04/04.11_integrate_LIGER_UINMF_full_mouse.log logs/04/04.02_integrate_SAMap_full_mouse.log \
           logs/04/04.03_integrate_CCA_genes6.log logs/04/04.06_integrate_RPCA_genes6.log \
           logs/04/04.01_integrate_fastMNN_genes6.log logs/04/04.04_integrate_Harmony_genes6.log logs/05/FAILED; do
    [ -f "$f" ] || continue
    grep -aE "^==|STOP|^FAILED|^done|^exit [0-9]|CHOICE READY|PIPELINE FINISHED|SKIP:|LABEL_COLS=" "$f" 2>/dev/null |
      sed 's/nice -n 10 //; s#[^ ]*/envs/[a-z]*/bin/python ##; s/ > logs.*//' |
      while IFS= read -r line; do
        key="$f|$line"
        grep -qxF -- "$key" $SEEN || { echo "$key" >> $SEEN; ev "[${f#logs/}] $line"; }
      done
  done

  # ---- possible stalls: a running 04/05 job of this project whose own log (its stdout) has not changed for
  #      STALL_HOURS. The log is read from /proc/<pid>/fd/1, so reruns writing elsewhere (label sets, rescores)
  #      and other projects' jobs with the same names cannot be confused with it
  for pid in $(pgrep -u "$USER" -f '0[45]\.[0-9]+_[A-Za-z_]+\.(R|py) '); do
    [ "$(readlink /proc/$pid/cwd 2>/dev/null)" = "$PWD" ] || continue
    case "$(cat /proc/$pid/comm 2>/dev/null)" in python*|R|Rscript) ;; *) continue ;; esac   # not the bash wrappers
    log=$(readlink /proc/$pid/fd/1 2>/dev/null); [ -f "$log" ] || continue
    [ -n "$(find "$log" -mmin +$((STALL_HOURS * 60)) 2>/dev/null)" ] || continue
    key="stall|$log|$(stat -c %Y "$log")"
    grep -qxF -- "$key" $SEEN || { echo "$key" >> $SEEN; ev "POSSIBLE STALL: ${log#$PWD/} unchanged for >${STALL_HOURS}h while its job runs"; }
  done

  # ---- status snapshot ----
  {
    echo "breast_cancer status at $(date '+%F %T')   (events: logs/EVENTS.log)"
    uptime | sed 's/.*load/load/'
    echo; echo "screens:"; screen -ls | grep -E "bc_" | sed 's/^/  /'
    echo; echo "pipeline:"; grep -aE "^==|STOP" logs/pipeline.log 2>/dev/null | tail -3 | sed 's/^/  /'
    echo; echo "per dataset (04 outputs newer than 03 / 17, eval_cells, 05 scores / 19):"
    for g in $DATASETS; do
      ref=data/$g/03_concatenated_data.h5ad
      n4=$(find data/$g -maxdepth 1 -name '04.*.h5ad' -newer $ref 2>/dev/null | wc -l)
      e=$([ data/$g/eval_cells.txt -nt $ref ] && echo yes || echo no)
      n5=$(find data/$g -maxdepth 1 -name '05.*.csv' -newer $ref 2>/dev/null | wc -l)
      printf '  %-16s 04 %2s/17  eval_cells %-3s  05 %2s/19\n' $g $n4 $e $n5
    done
    echo; echo "04 jobs: done $(grep -cE '^done' logs/04/run_04.log 2>/dev/null), failed $(grep -cE '^FAILED' logs/04/run_04.log 2>/dev/null)"
    grep -E '^FAILED' logs/04/run_04.log 2>/dev/null | sed 's/nice -n 10 //; s#/dskh[^ ]*python ##; s/ > .*//; s/^/  /'
    echo; echo "last 10 events:"; tail -10 $EV | sed 's/^/  /'
  } > $ST.tmp && mv $ST.tmp $ST

  sleep 300
done
