#!/bin/bash
#$ -S /bin/bash
set -e

#set -e -x

#######################################################################
#
#  Program:   ASHSTHK (Multi template thickness pipeline for ASHS)
#  Module:    $Id$
#  Language:  BASH Shell Script
#  Copyright (c) 2021 Long Xie, University of Pennsylvania
#
#  This file is the implementation of the multi-template thickness pipeline
#  to measure thickness of medial temporal lobe subregions. The software is 
#  distributed using Docker. Internet access is required.
#
#######################################################################

# some basic functions
function usage()
{
  cat <<-USAGETEXT

ashsthk_docker: multi-template thickness pipeline for ASHS
  usage:
    ashsthk_docker [options]

  required options:
    -n str            Subject ID
    -l str            Side (left or right)
    -i path           Path to the ASHS segmentation.
    -a path           Path to the multi-template thickness template
    -w path           Output directory

  optional:
    -T                Tidy mode. Cleans up files once they are unneeded.
    -R                Do not remove Docker image when finish. Use this argument if you 
                      need to run the script on more than one subjects to avoid loading
                      Docker image at every run. The user need to remove the Docker image
                      manually. 
    -d str            User can specify specific Docker image (Default: longxie/ashsthk:v3)
    -c number         Number of CPUs the docker container can use (Default: all CPUs)
    -u id             User can  specify a user ID for the docker image to avoid the output 
                      files being own by the root (Default: id -u)
    -g id             User can specify a group ID for the user ID given by the -u argument.
                      (Default: id -g). 
    -h                Print help
    -s integer        Run only one stage (see below); also accepts range (e.g. -s 1-3).
                      By default, Only steps 1 to 5 will be run (to variant template).
                      If the user pointwise correspondance is desired, steps 6-8 are needed.
                      Stages:
                        1: Perform affine and coarse deformable registration between
                           subject segmentation to all the atlases.
                        2: Determine group membership.
                        3: Perform deformable registration to the selected variant template.
                        4: Perform geodesic shooting to the variant template (VT)
                        5: Evaluate fit quality and measure thickness for VT.
                        6: Perform deformable registration to the unified template.
                        7: Perform geodesic shooting to the unified template (UT).
                        8: Evaluate fit quality and measure thickness for UT.

USAGETEXT
}

# Dereference a link - different calls on different systems
function dereflink ()
{
  if [[ $(uname) == "Darwin" ]]; then
    local SLTARG=$(readlink $1)
    if [[ $SLTARG ]]; then
      echo $SLTARG
    else
      echo $1
    fi
  else
    readlink -m $1
  fi
}

# Print usage by default
if [[ $# -lt 1 ]]; then
  echo "Try $0 -h for more information."
  exit 2
fi

# Read the options
unset DELETETMP DOCKERIMG RMDOCKERIMG USERID GROUPID
CPUS="DEFAULT"
while getopts "n:l:i:a:w:s:d:c:u:g:hTR" opt; do
  case $opt in

    n) id=$OPTARG;;
    l) SIDE=$OPTARG;;
    i) SEGORIG=$OPTARG;;
    a) GSTEMPDIR=$OPTARG;;
    w) OUTDIR=$OPTARG;;
    s) STAGE_SPEC=$OPTARG;;
    d) DOCKERIMG=$OPTARG;;
    u) USERID=$OPTARG;;
    g) GROUPID=$OPTARG;;
    c) CPUS=$OPTARG;;
    T) DELETETMP=1;;
    R) RMDOCKERIMG=0;;
    h) usage; exit 0;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;

  esac
done

##############################################
# Software PATH
if [[ $DELETETMP == "" ]]; then
  DELETETMP="0"
fi

if [[ $RMDOCKERIMG == "" ]]; then
  RMDOCKERIMG="1"
fi

if [[ $DOCKERIMG == "" ]]; then
  DOCKERIMG="longxie/ashsthk:v3"
fi

# Check if the required parameters were passed in
echo "id    : ${id?    "Subject id was not specified. See $0 -h"}"
echo "side  : ${SIDE?  "The side of the subject was not specified. See $0 -h"}"
echo "ASHS segmentation  : ${SEGORIG?  " The path to ASHS segmentation was not specified. See $0 -h"}"
echo "Tempalte  : ${GSTEMPDIR? "The path to the multi-template template was not specified. See $0 -h"}"
echo "OutputDir    : ${OUTDIR?    "CSV file with information of all the timepoints was not specified. See $0 -h"}"

# Check whether the variable Side is valide.
if [[ $SIDE != "left" && $SIDE != "right" ]]; then
  echo "The input to -l has to be left or right. See $0 -h."
  exit -2
fi

# Convert the work directory to absolute path
mkdir -p ${OUTDIR?}
OUTDIR=$(cd $OUTDIR; pwd)
if [[ ! -d $OUTDIR ]]; then
  echo "Work directory $OUTDIR cannot be created"
  exit -2
fi

# Check the CPU setting
re="^[0-9]+([.][0-9]+)?$"
CPUCMD=" "
if [[ $CPUS == "DEFAULT" ]]; then
  CPUCMD=" "
elif ! [[ $CPUS =~ $re ]]; then
  echo "Number of CPUs need to be a positive number"
  exit -2
else
  CPUCMD=" --cpus $CPUS "
fi

# Check the user and group ID option
if [[ $USERID == "" ]]; then
  USERID=$(id -u)
fi
if [[ $GROUPID == "" ]]; then
  GROUPID=$(id -g)
fi

# Set the start and end stages
if [[ $STAGE_SPEC ]]; then
  STAGE_START=$(echo $STAGE_SPEC | awk -F '-' '$0 ~ /^[0-9]+-*[0-9]*$/ {print $1}')
  STAGE_END=$(echo $STAGE_SPEC | awk -F '-' '$0 ~ /^[0-9]+-*[0-9]*$/ {print $NF}')
  if [[ $STAGE_START -lt 1 ]]; then
    STAGE_START=1
  fi
  if [[ $STAGE_END -gt 8 ]]; then
    STAGE_END=8
  fi
  if [[ $STAGE_END -lt $STAGE_START ]]; then
    STAGE_END=$STAGE_START
  fi
else
  STAGE_START=1
  STAGE_END=5
fi

#############################################
function main()
{
  # import docker image
  sudo docker pull $DOCKERIMG

  # copy segmentation into the output folder
  cp $SEGORIG $OUTDIR/${id}_${SIDE}_ASHSSeg.nii.gz

  # generate a docker container and run ashsthk
  CMDT=""
  if [[ $DELETETMP == "1" ]]; then
    CMDT=" -T "
  fi
  sudo docker run $CPUCMD \
    --name=ashsthk-${id}-${SIDE} \
    --mount type=bind,source="${GSTEMPDIR}",target=/app/template,readonly \
    --mount type=bind,source="${OUTDIR}",target=/app/output \
    -e OUT_UID=${USERID} \
    -e OUT_GID=${GROUPID} \
    $DOCKERIMG \
    /bin/bash /home/ashsthk/ashsthk_main.sh \
      -n $id \
      -i /app/output//${id}_${SIDE}_ASHSSeg.nii.gz \
      -l $SIDE \
      -a /app/template \
      -s ${STAGE_START}-${STAGE_END} \
      $CMDT \
      -w /app/output/
    

  # remove the docker container and image
  sudo docker rm ashsthk-${id}-${SIDE}
  if [[ $RMDOCKERIMG == "1" ]]; then
    sudo docker rmi $DOCKERIMG
  fi
}

main






