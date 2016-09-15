#!/bin/bash
# endless loop renicing cluster every 30 minutes
# needs to be run as either sudo or cluster

while [ 0 ]
do
renice +10 -u cluster
sleep 1800
done

