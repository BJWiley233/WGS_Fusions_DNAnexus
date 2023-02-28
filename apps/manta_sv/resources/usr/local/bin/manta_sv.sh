set -o pipefail
set -o errexit

bam=$1
ref=$2
config=$3
non_wgs=$4
outputContig=$5

add_to_command=""
if [[ $non_wgs == "true" ]]; then
   add_to_command="$add_to_command --exome"
fi
if [[ $outputContig == "true" ]]; then
   add_to_command="$add_to_command --outputContig"
fi

/usr/bin/python /usr/bin/manta/bin/configManta.py \
   --config ${config} \
   --referenceFasta ${ref} \
   --tumorBam ${bam} \
   --runDir $(pwd) $add_to_command

/usr/bin/python runWorkflow.py -m local -g 20 -j8

