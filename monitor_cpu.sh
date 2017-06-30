#!/bin/bash

PROC=$1
N=600

ps aux | head -n 1
while true; do
  ps aux | grep "$PROC" | grep -v "grep" | grep -v "monitor_cpu.sh"
  sleep $N
done
