#!/bin/sh

while IFS='' read -r line || [[ -n "$line" ]]; do
  filesarray+=( "$line" )
done < histostoplot.txt

for i in "${filesarray[@]}"
do
  echo $i
  python plot.py --var $i --logScale 1
done
