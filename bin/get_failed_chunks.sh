cohort=$1
fileprefix=$2
chunk_file=$3

while read line
do
  chunk=`echo $line | sed 's:-:_:g;s/:/_/g'`
  if [ ! -f ${fileprefix}${chunk}.h5 ]; then
    echo $chunk
  fi
done < ${chunk_file} 


