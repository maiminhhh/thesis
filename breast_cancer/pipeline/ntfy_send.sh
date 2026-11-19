#!/usr/bin/env bash
# send a phone notification via ntfy.sh:  bash pipeline/ntfy_send.sh "<title>" "<message>" [priority 1-5]
# topic in pipeline/ntfy_topic (random, private by obscurity; subscribe to it in the ntfy app)
T=$(cat "$(dirname "$0")/ntfy_topic")
curl -sS -m 20 -o /dev/null -w "%{http_code}" -H "Title: $1" -H "Priority: ${3:-3}" -H "Tags: dna" -d "$2" "https://ntfy.sh/$T"
