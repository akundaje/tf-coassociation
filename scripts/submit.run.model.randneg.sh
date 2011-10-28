#!/bin/bash
# Will batch submit several runs for all files of a particular type in a directory

if [[ "$#" -lt 4 ]]
    then
    echo $(basename $0) 1>&2
    echo "Computes percentage tags in peaks as a measure of IP enrichment" 1>&2
    echo "USAGE:" 1>&2
    echo "$(basename $0) [scriptName] [inputDir] [inputRegex] [outputDir]" 1>&2
    echo "   [scriptName]: script to be parallelized" 1>&2
    echo "   [inputDir]: input directory containing files" 1>&2
    echo "   [inputRegex]: regular expression used to match file names" 1>&2
    echo "   [outputDir]: output directory" 1>&2
    echo "   [niter]: (OPTIONAL) number of sampled models, DEFAULT:50" 1>&2
    echo "   [mem]: (OPTIONAL) memory in GB, Default: 8" 1>&2
    exit 1 
fi

SCRIPT_NAME=$1
if [[ ! -e ${SCRIPT_NAME} ]]
    then
    echo "Script ${SCRIPT_NAME} does not exist" 1>&2
    exit 1
fi

IDIR=$2
if [[ ! -d ${IDIR} ]]
    then
    echo "Input directory ${IDIR} does not exist" 1>&2
    exit 1
fi

IREGEX=$3

ODIR=$4
[[ ! -d ${ODIR} ]] && mkdir ${ODIR}

NITER=50
if [[ "$#" -gt 4 ]]
    then
    NITER=$5
fi

MEM=4
if [[ "$#" -gt 5 ]]
    then
    MEM=$6
fi
MEM=$(( MEM * 1024 ))

for iFile in $(find ${IDIR} -regextype posix-extended -type f -regex ${IREGEX})
  do
  baseFile=$(echo $(basename ${iFile}) | sed -r -e 's/^SigMtrx_//g' -e 's/\.mtrx.+$//g') #input file name prefix
  echo ${baseFile}
  for (( i=1; i<=${NITER}; i++ ))
    do
    OUTDIR="${ODIR}/${baseFile}" # Create a directory with input file prefix
    [[ ! -d ${OUTDIR} ]] && mkdir ${OUTDIR}
    OUTPREFIX="${OUTDIR}/${baseFile}.iter${i}"
    OUTFILE="${OUTPREFIX}.randneg.rfresults.Rdata"
    logFile="${OUTPREFIX}.out"
    errFile="${OUTPREFIX}.err"
    TMP_DIR="${TMP}/rulefit_${RANDOM}${RANDOM}"
    #bsub -M "${MEM}" -R "rusage[mem=${MEM}]" -o ${logFile} -e ${errFile} "Rscript ${SCRIPT_NAME} ${iFile} ${OUTFILE}"
    if [[ ! -e ${OUTFILE} ]]
	then
	#bsub -q research-rh6 -W 94:00 -M ${MEM} -R "rusage[mem=${MEM}]" -o ${logFile} -e ${errFile} "mkdir ${TMP_DIR} ; Rscript ${SCRIPT_NAME} ${iFile} ${OUTFILE} 0 ${TMP_DIR} ; rm -rf ${TMP_DIR}"
	bsub -q research-rh6 -W 94:00 -o ${logFile} -e ${errFile} "mkdir ${TMP_DIR} ; Rscript ${SCRIPT_NAME} ${iFile} ${OUTFILE} 0 ${TMP_DIR} ; rm -rf ${TMP_DIR}"
    fi
  done
done
