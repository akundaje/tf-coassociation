#!/bin/bash
# Will batch submit several runs for all files of a particular type in a directory

if [[ "$#" -lt 4 ]]
    then
    echo $(basename $0) 1>&2
    echo "Computes percentage tags in peaks as a measure of IP enrichment" 1>&2
    echo "USAGE:" 1>&2
    echo "$(basename $0) [scriptName] [inputDir] [inputRegex] [outputDir]" 1>&2
    echo "   [scriptName]: script to be parallelized" 1>&2
    echo "   [posDir]: input directory containing positive set files" 1>&2
    echo "   [negDir]: input directory containing negative set files" 1>&2
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

POSDIR=$2
if [[ ! -d ${POSDIR} ]]
    then
    echo "Positive set directory ${POSDIR} does not exist" 1>&2
    exit 1
fi

NEGDIR=$3
if [[ ! -d ${NEGDIR} ]]
    then
    echo "Negative set directory ${NEGDIR} does not exist" 1>&2
    exit 1
fi

IREGEX=$4

ODIR=$5
[[ ! -d ${ODIR} ]] && mkdir ${ODIR}

NITER=50
if [[ "$#" -gt 5 ]]
    then
    NITER=$6
fi

MEM=4
if [[ "$#" -gt 6 ]]
    then
    MEM=$7
fi
MEM=$(( MEM * 1024 ))

for posFile in $(find ${POSDIR} -regextype posix-extended -type f -regex ${IREGEX})
  do
  baseFile=$(echo $(basename ${posFile}) | sed -r -e 's/^SigMtrx_//g' -e 's/\.mtrx.+$//g') #input file name prefix
  echo ${baseFile}
  negFile=$( find ${NEGDIR} -type f -name "*${baseFile}*.Rdata" )
  if [[ ! -e ${negFile} ]]
      then
      echo "SKIPPING:Negative set File not found" 1>&2
      continue
  fi
  OUTDIR="${ODIR}/${baseFile}" # Create a directory with input file prefix
  [[ ! -d ${OUTDIR} ]] && mkdir ${OUTDIR}
  OUTPREFIX="${OUTDIR}/${baseFile}"
  OUTFILE="${OUTPREFIX}.posneg.rfresults.Rdata"
  logFile="${OUTPREFIX}.out"
  errFile="${OUTPREFIX}.err"
  #bsub -M "${MEM}" -R "rusage[mem=${MEM}]" -o ${logFile} -e ${errFile} "Rscript ${SCRIPT_NAME} ${iFile} ${OUTFILE}"
  if [[ ! -e ${OUTFILE} ]]
      then
      if [[ ${MEM} -gt $(( 4 * 1024)) ]]
	  then
	  bsub -W 94:00 -M "${MEM}" -R "rusage[mem=${MEM}]" -o ${logFile} -e ${errFile} "Rscript ${SCRIPT_NAME} ${posFile} ${negFile} ${OUTFILE}"
      else
	  bsub -W 94:00 -o ${logFile} -e ${errFile} "Rscript ${SCRIPT_NAME} ${posFile} ${negFile} ${OUTFILE}"
      fi
  fi
done

