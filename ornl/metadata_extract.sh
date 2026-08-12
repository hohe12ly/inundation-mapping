#!/bin/bash
# usage: cat file | metadata_extract.sh
# usage: unzip -l xxxxx.zip | metadata_extract.sh

[ ! -p /dev/stdin ] && echo "only take piped input" && exit 0 # only take piped input

source `dirname $0`/metadata.sh

while IFS= read aline
do
  #echo "IN: $aline"
  fn=`echo "$aline" | awk -F'/' '{print $NF}'`
  [ -z "$fn" ] && echo "$aline" && continue
  m=$(filemeta $fn)
  [ -z "$m" ] && echo "$aline" && continue
  newline=`echo "$aline" | sed -e "s/$fn/$fn : $m/"`
  echo "$newline"
done

#ddir=/gpfs/corral3/repl/utexas/NFIE-NWM-UT/nfiedata/pin2flood/2022-July

#h="$1"

#filelist=`grep -E "$h/[a-z][a-z]*" $ddir/${h}.txt | awk -F'/' '{print $NF}'
