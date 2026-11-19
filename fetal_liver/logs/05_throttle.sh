#!/usr/bin/env bash
# adjusts the running 05 xargs (-P) to the server load every 10 min: 15-min load > 95 -> one
# fewer job (min 6); < 65 -> one more (max 16). Exits when that xargs is gone.
pid=$1; P=$2
while kill -0 $pid 2>/dev/null; do
  sleep 600
  load=$(awk '{print int($3)}' /proc/loadavg)
  if [ $load -gt 95 ] && [ $P -gt 6 ]; then kill -USR2 $pid && P=$((P - 1)) && echo "$(date '+%F %T') load $load -> $P jobs"
  elif [ $load -lt 65 ] && [ $P -lt 16 ]; then kill -USR1 $pid && P=$((P + 1)) && echo "$(date '+%F %T') load $load -> $P jobs"; fi
done >> logs/05_throttle.log
