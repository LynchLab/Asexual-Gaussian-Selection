#!/usr/bin/env bash

# Script to launch a batch of runs

PROGRAM_NAME=Gaussian

CURRENT_DATETIME=`date +%Y.%m.%d.%H_%M_%S`
RUN_DIRECTORY=runs/${CURRENT_DATETIME}
ROOT_DIRECTORY=`pwd`

# more legible to see TRUE/FALSE vs 1/0
TRUE=1
FALSE=0

# Load the Modules we need for compiling and running
source Modules.sh

# Create the ./runs directory
mkdir -p ${RUN_DIRECTORY}
cd ${RUN_DIRECTORY}

# Copy source code
RUN_PROGRAM_NAME=_${PROGRAM_NAME}_For_Run

# Remove Comments from Parameters
REFERENCE_PARAMETERS=Parameters_For_Reference
cp ${ROOT_DIRECTORY}/Parameters ./${REFERENCE_PARAMETERS}
awk '!/^#|^$/ || /^#!/' ${REFERENCE_PARAMETERS} > .${REFERENCE_PARAMETERS}

VARIABLE_NAMES=`awk '{ print $1 }' .${REFERENCE_PARAMETERS}`

declare -A VARIABLE_VALUES
declare -A RUN_INDICES
declare -A RUN_INDEX_COUNTS

for variable_name in ${VARIABLE_NAMES}
do
  VARIABLE_VALUES[$variable_name]=`grep ${variable_name} .${REFERENCE_PARAMETERS} | awk '{sub($1 FS,"")}7'`
  RUN_INDICES[$variable_name]=1 # annoyingly awk starts counting columns at 1 not 0
  RUN_INDEX_COUNTS[$variable_name]=`echo ${VARIABLE_VALUES[$variable_name]} | awk ' { print NF } '`
done

# iterate over all combinations of VARIABLE_VALUES and make substitutions

more_runs_to_prepare=TRUE

while [ $more_runs_to_prepare ]
  do

  for variable_name in ${VARIABLE_NAMES}
  do
    echo -n $variable_name:${RUN_INDICES[$variable_name]} " "
  done
  echo

  # create the directory for the run
  CURRENT_RUN_DIRECTORY="RUN"
  for variable_name in ${VARIABLE_NAMES}
  do
    CURRENT_RUN_DIRECTORY=${CURRENT_RUN_DIRECTORY}_${variable_name}_${RUN_INDICES[$variable_name]}
  done

  CURRENT_RUN_DIRECTORY=`echo $CURRENT_RUN_DIRECTORY|sed s/[[:space:]]/_/g`

  mkdir -p $CURRENT_RUN_DIRECTORY
  cd $CURRENT_RUN_DIRECTORY

  # copy the original program to the run
  cp ${ROOT_DIRECTORY}/${PROGRAM_NAME}.cpp ./${RUN_PROGRAM_NAME}.cpp
  cp ${ROOT_DIRECTORY}/Gaussian.sh ./Gaussian.sh

  for variable_name in ${VARIABLE_NAMES}
  do
    TEMPFILE=`mktemp --suffix .cpp`
    variable_value=`echo ${VARIABLE_VALUES[$variable_name]} | awk -v column=${RUN_INDICES[$variable_name]} '{ print $column }'`
    awk -v variable_name="$variable_name" -v variable_value="$variable_value" '{ $3 = (($1 == "#define" && $2 == variable_name) ? variable_value : $3) } " "' ${RUN_PROGRAM_NAME}.cpp > $TEMPFILE
    cp $TEMPFILE ${RUN_PROGRAM_NAME}.cpp
    rm $TEMPFILE
  done

  # Compile the program
  echo Compiling..
  icpx ${RUN_PROGRAM_NAME}.cpp -lm `pkg-config --cflags gsl --libs` -lgsl -lgslcblas -o ${PROGRAM_NAME} -O3

  echo Submitting Job to Queue
  sbatch --array=1-21 Gaussian.sh

  cd ..

  INCREMENT_VARIABLE=`echo ${VARIABLE_NAMES[@]} |awk '{ print $1 }'`
  RUN_INDICES[$INCREMENT_VARIABLE]=$((${RUN_INDICES[$INCREMENT_VARIABLE]}+1))

  carry=FALSE
  for variable_name in ${VARIABLE_NAMES}
  do
    if (($carry))
    then
      RUN_INDICES[$variable_name]=$((${RUN_INDICES[$variable_name]}+1))
      carry=FALSE
    fi

    if (( ${RUN_INDICES[$variable_name]} > ${RUN_INDEX_COUNTS[$variable_name]}))
    then
      RUN_INDICES[$variable_name]=1
      carry=TRUE
    fi
  done

  if (($carry))
  then
    echo "Job preparation complete"
    break
  fi
done