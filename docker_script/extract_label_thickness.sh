#!/bin/bash
#$ -S /bin/bash
set -e

#set -e -x

#######################################################################
#
#  Program:   extract_label_thickness
#  Module:    $Id$
#  Language:  BASH Shell Script
#  Copyright (c) 2021 Long Xie, University of Pennsylvania
#
#  This script will extract mean thickness of each label in the 
#  input segmentation. This method works best for sheet-like structures
#  (like the cortex) and may not work for tubular structures.
#
#######################################################################

# some basic functions
function usage()
{
  cat <<-USAGETEXT

extract_label_thickness: extract mean thickness of each label
  usage:
    extract_label_thickness [options]

  required options:
    -i path           Path to the segmentation.
    -w path           Output directory

  optional:
    -T                Tidy mode. Cleans up files once they are unneeded.
    -e N              Minimal number of mesh edges separating two generator
		      points of a VD face for it to be considered (default: 6)
    -p X.XX           Prune the mesh using factor X.XX (try 2.0). The
                      pruning algorithm deletes faces in the VD for
                      which the ratio of the geodesic distance between
                      the generating points and the euclidean distance
		      between these points is less than X.XX (default: 1.2)
    -R                Do not remove Docker image when finish. Use this argument if you
                      need to run the script on more than one subjects to avoid loading
                      Docker image at every run. The user need to remove the Docker image
                      manually.
    -d str            User can specify specific Docker image (Default: longxie/ashsthk:v3)
    -u id             User can  specify a user ID for the docker image to avoid the output
                      files being own by the root (Default: id -u)
    -g id             User can specify a group ID for the user ID given by the -u argument.
                      (Default: id -g).
    -c number         Number of CPUs the docker container can use (Default: all CPUs)
    -h                Print help

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
while getopts "i:w:d:c:e:p:u:g:hRT" opt; do
  case $opt in

    i) SEGORIG=$OPTARG;;
    w) OUTDIR=$OPTARG;;
    e) thick_e=$OPTARG;;
    p) thick_p=$OPTARG;;
    c) CPUS=$OPTARG;;
    d) DOCKERIMG=$OPTARG;;
    u) USERID=$OPTARG;;
    g) GROUPID=$OPTARG;;
    T) DELETETMP=1;;
    R) RMDOCKERIMG=0;;
    h) usage; exit 0;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;

  esac
done

##############################################
# Setup environment
if [[ $RMDOCKERIMG == "" ]]; then
  RMDOCKERIMG="1"
fi

if [[ $DOCKERIMG == "" ]]; then
  DOCKERIMG="longxie/ashsthk:v3"
fi

# Check if the required parameters were passed in
echo "Input segmentation  : ${SEGORIG?  " The path to input segmentation was not specified. See $0 -h"}"
echo "OutputDir    : ${OUTDIR?    "Output directory was not specified. See $0 -h"}"

filenameext=$(basename $SEGORIG)
filename=${SEGORIG%.*}
filename=${filename%.*}
filename=$(basename $filename)

if [[ $DELETETMP == "" ]]; then
  DELETETMP="0"
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

# Check value of the -p and -e command
if [[ $thick_p == "" ]]; then
  thick_p=1.2
fi
if [[ $thick_e == "" ]]; then
  thick_e=6
fi


#############################################
function main()
{
  # import docker image
  sudo docker pull $DOCKERIMG

  # copy segmentation into the output folder
  cp $SEGORIG $OUTDIR/${filenameext}

  # generate a docker container and run ashsthk
  CMDT=""
  if [[ $DELETETMP == "1" ]]; then
    CMDT=" -T "
  fi
  sudo docker run $CPUCMD \
    --name=extract-label-thickness-${filename} \
    --mount type=bind,source="${OUTDIR}",target=/app/output \
    -e OUT_UID=${USERID} \
    -e OUT_GID=${GROUPID} \
    $DOCKERIMG \
    /bin/bash /home/utils/extract_label_thickness.sh \
      -i /app/output/${filenameext} \
      $CMDT \
      -p $thick_p \
      -e $thick_e \
      -w /app/output/

  # remove the docker container and image
  sudo docker rm extract-label-thickness-${filename}
  if [[ $RMDOCKERIMG == "1" ]]; then
    sudo docker rmi $DOCKERIMG
  fi
}

main


