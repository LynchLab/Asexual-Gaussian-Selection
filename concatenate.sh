#!/usr/bin/env bash

cd runs
batches=20*/

for batch in $batches
do
  batch=`echo $batch | sed s/\\\///g`
  echo > $batch.txt
  runs=${batch}/RUN*
  for run in $runs
  do
    awk 'ENDFILE { if (FNR == 0) printf("\n") >>FILENAME }' $run/dataout*.txt 2>/dev/null | sort -k1n >> $batch.txt
    cat $run/dataout*.txt 2>/dev/null | sort -k1n >> $batch.txt
    echo >> $batch.txt
  done

  echo >> $batch.txt
  echo >> $batch.txt
done