#!/bin/bash
JOBLIST=""
TASKLIST[0]=""
i=0

TASK="./run-one.sh issams r2/Q0"
echo $TASK
$TASK &
LASTJOB=$!
JOBLIST="$JOBLIST $LASTJOB"
TASKLIST[$i]=$TASK
i=$(($i+1))

TASK="./run-one.sh issams r2/Q1"
echo $TASK
$TASK &
LASTJOB=$!
JOBLIST="$JOBLIST $LASTJOB"
TASKLIST[$i]=$TASK
i=$(($i+1))

TASK="./run-one.sh issams-dev r2/Q0"
echo $TASK
$TASK &
LASTJOB=$!
JOBLIST="$JOBLIST $LASTJOB"
TASKLIST[$i]=$TASK
i=$(($i+1))

TASK="./run-one.sh issams-dev r2/Q1"
echo $TASK
$TASK &
LASTJOB=$!
JOBLIST="$JOBLIST $LASTJOB"
TASKLIST[$i]=$TASK
i=$(($i+1))

# now wait for all the processes to finish
i=0
for JOB in $JOBLIST
do
   wait $JOB
   echo “${TASKLIST[$i]} – Job $JOB exited with status $?”
   i=$(($i+1))
done

