#!/bin/bash
#$ -S /bin/bash
#set -e

set -e -x

#######################################################################
#
#  Program:   ASHSTHK (Multi template thickness pipeline for ASHS)
#  Module:    $Id$
#  Language:  BASH Shell Script
#  Copyright (c) 2021 Long Xie, University of Pennsylvania
#
#  This file is the implementation of the multi-template thickness pipeline 
#  to measure thickness of medial temporal lobe subregions.
#
#######################################################################

# some config that is tailored for the docker
export LD_LIBRARY_PATH=/home/customlib/lib:/home/customlib/icclibs
export ASHSTHKROOT=/home/ashsthk
export MATLAB_RT=/home/pkg/MATLAB_2016a/v901
#NSLOTS=1



# some basic functions
function usage()
{
  cat <<-USAGETEXT

ashsthk_main: multi-template thickness pipeline for ASHS
  usage:
    ashsthk_main [options]

  required options:
    -n str            Subject ID
    -l str            Side (left or right)
    -i path           Path to the ASHS segmentation.
    -a path           Path to the multi-template thickness template
    -w path           Output directory

  optional:
    -T                Tidy mode. Cleans up files once they are unneeded.
    -h                Print help
    -s integer        Run only one stage (see below); also accepts range (e.g. -s 1-3).
		      By default, Only steps 1 to 5 will be run (to variant template).
                      If the user needs pointwise correspondance, steps 6 to 8 are needed. 
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
while getopts "n:l:i:a:w:s:hT" opt; do
  case $opt in

    n) id=$OPTARG;;
    l) SIDE=$OPTARG;;
    i) SEGORIG=$OPTARG;;
    a) GSTEMPDIR=$OPTARG;;
    w) OUTDIR=$OPTARG;;
    s) STAGE_SPEC=$OPTARG;;
    T) DELETETMP=1;;
    h) usage; exit 0;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;

  esac
done

##############################################
# Setup environment
#BASEDIR=$(dirname "$0")
#export ASHSTHKROOT=$BASEDIR

# Software PATH
CODEDIR=$ASHSTHKROOT
BINDIR=$ASHSTHKROOT/bin
C3DPATH=$ASHSTHKROOT/bin
FLIPLRMAT=$BINDIR/flip_LR.mat
FLIPLRITK=$BINDIR/flip_LR_itk.txt
WORKDIR=$OUTDIR/work_${SIDE}
DUMPDIR=$WORKDIR/dump
TMPDIR=$DUMPDIR/tmp
if [[ $DELETETMP == "" ]]; then
  DELETETMP="0"
fi

# GS template directories
ATLASLIST=$GSTEMPDIR/GSTemplate/MST/paths/IDSide.txt
GTGROUPS=$GSTEMPDIR/group_xval_2group_all.txt

# Check if the required parameters were passed in
echo "id    : ${id?    "Subject id was not specified. See $0 -h"}"
echo "side  : ${SIDE?  "The side of the subject was not specified. See $0 -h"}"
echo "ASHS segmentation  : ${SEGORIG?  " The path to ASHS segmentation was not specified. See $0 -h"}"
echo "Tempalte  : ${GSTEMPDIR? "The path to the multi-template template was not specified. See $0 -h"}"
echo "OutputDir    : ${OUTDIR?    "CSV file with information of all the timepoints was not specified. See $0 -h"}"

# Check the root dir
if [[ ! $ASHSTHKROOT ]]; then
  echo "Please set ASHSTHKROOT to the ASHS thicknes pipeline root directory before running $0"
  exit -2
elif [[ $ASHSTHKROOT != $(dereflink $ASHSTHKROOT) ]]; then
  echo "ASHSTHKROOT must point to an absolute path, not a relative path"
  exit -2
fi

# Check matlab directory
if [[ ! $MATLAB_RT ]]; then
  echo "Please set MATLAB_RT to the MATLAB rumtime directory before running $0"
  exit -2
fi

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
mkdir -p $WORKDIR
mkdir -p $DUMPDIR
mkdir -p $TMPDIR


# read parameters from the template config file.
if [[ ! -f $GSTEMPDIR/template_config.sh ]]; then
  echo "The configuration file of the multi-template template can not be found. Please check the completeness of the template."
  exit -2
fi
source $GSTEMPDIR/template_config.sh

# Redirect output/error to a log file in the dump directory
LOCAL_LOG=$(date +ashsthk_main.o%Y%m%d_%H%M%S)
mkdir -p $DUMPDIR
exec > >(tee -i $DUMPDIR/$LOCAL_LOG)
exec 2>&1

# Write into the log the arguments and environment
echo "ashsthk_main execution log"
echo "  timestamp:   $(date)"
echo "  invocation:  $0 $@"
echo "  directory:   $PWD"
echo "  environment:"

STAGE_NAMES=(\
  "Perform affine and coarse deformable registration between subject segmentation to all the atlases." \
  "Determine group membership." \
  "Perform deformable registration to the selected variant template." \
  "Perform geodesic shooting to the variant template (VT)." \
  "Evaluate quality of fit and measure thickness for VT." \
  "Perform deformable registration to the unified template." \
  "Perform geodesic shooting to the unified template (UT)." \
  "Evaluate quality of fit and measure thickness for UT.")

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
  # Run the various stages
  for ((STAGE=$STAGE_START; STAGE<=$STAGE_END; STAGE++)); do

    # The desription of the current stage
    STAGE_TEXT=${STAGE_NAMES[STAGE-1]}
    echo "****************************************"
    echo "Starting stage $STAGE: $STAGE_TEXT"
    echo "****************************************"

    case $STAGE in

      1) RegToAtlases $SIDE;;
      2) TempMembership $SIDE;;
      3) RegToInitTemp $SIDE;;
      4) GeoShooting $SIDE;;
      5) PostSteps $SIDE;;
      6) RegToUT $SIDE;;
      7) GeoShootingUT $SIDE;;
      8) PostSteps $SIDE;;

    esac

  done

  DeleteFile

  echo "Done!"
}


######################################################
function RegToAtlases()
{
  side=$1

  ########################################
  # resample the images GS template chunk
  ########################################
  SUBJDATADIR=$WORKDIR/data
  mkdir -p $SUBJDATADIR
  #SEGORIG=$SUBJASHSDIR/final/${id}_${side}_lfseg_heur.nii.gz
  SEG=$SUBJDATADIR/${id}_${side}_lfseg_heur_dividedCS.nii.gz
  if [[ ! -f $SEGORIG ]]; then
    echo "$SEGORIG does not exist."
    exit
  fi


  # divide the CS label
  if [ ! -f $SEG ]; then
  $C3DPATH/c3d $SEGORIG \
    -replace $ERCLABEL 999 $BA35LABEL 999 $BA36LABEL 999 \
    -thresh 999 999 1 0 -sdt -scale -1 -popas A \
    $SEGORIG -thresh $PHCLABEL $PHCLABEL 1 0 -sdt -scale -1 -popas P \
    -push P -push A -vote -popas AP \
    $SEGORIG -thresh $CSLABEL $CSLABEL 1 0 -push AP -multiply \
    $SEGORIG -add \
    -o $SEG
  fi

  if [ ! -f $SUBJDATADIR/${id}_${side}_seg.nii.gz ]; then
    # Generate a binary image for the label
    #RM="-rm $SUBJDIR/${PREFIX1}_${id}_T1w_trim_denoised_SR.nii.gz $SUBJDATADIR/${id}_${side}_tse.nii.gz "
    RM=""
    for ((i=0;i<${#LABEL_IDS_FIT[*]};i++)); do
      $C3DPATH/c3d $SEG -replace $(for k in ${LABEL_MRG_FIT[i]}; do echo $k 999; done) -thresh 999 999 1 0 \
        -o $SUBJDATADIR/${id}_${side}_${LABEL_IDS_FIT[i]}_orig.nii.gz
      RM="$RM -rm $SUBJDATADIR/${id}_${side}_${LABEL_IDS_FIT[i]}_orig.nii.gz $SUBJDATADIR/${id}_${side}_${LABEL_IDS_FIT[i]}.nii.gz "
    done
    $C3DPATH/c3d $(for sub in ${LABEL_IDS_FIT[*]}; do echo $SUBJDATADIR/${id}_${side}_${sub}_orig.nii.gz; done) \
    -vote -type ushort -o $SUBJDATADIR/${id}_${side}_seg_orig.nii.gz

  # mlaffine to align subject seg to template space
  if [ ! -f $SUBJDATADIR/${id}_${side}_to_MSTInitTemp_mlaffine.txt ]; then

    REGSEG=$SUBJDATADIR/${id}_${side}_seg_orig.nii.gz
    ADD_ON_TRANS=""
    if [[ $side == "right" ]]; then
      # moments alignment
      ADD_ON_TRANSMAT=$TMPDIR/${id}_${side}_seg_moments.mat
      ADD_ON_TRANS=$TMPDIR/${id}_${side}_seg_moments_itk.txt
      $BINDIR/greedy -d 3  \
        -i $GSTEMPDIR/GSUTemplate/gshoot/template_1/template/iter_2/template_1_gshoot_seg.nii.gz \
           $REGSEG \
        -moments 2 \
        -o $ADD_ON_TRANSMAT

      $BINDIR/greedy -d 3  \
        -rf $GSTEMPDIR/GSUTemplate/gshoot/template_1/template/iter_2/template_1_gshoot_seg.nii.gz \
        -ri LABEL 0.2vox \
        -rm $SUBJDATADIR/${id}_${side}_seg_orig.nii.gz \
            $TMPDIR/${id}_${side}_seg_orig_flipLR.nii.gz \
        -r $ADD_ON_TRANSMAT
      REGSEG=$TMPDIR/${id}_${side}_seg_orig_flipLR.nii.gz
      $C3DPATH/c3d_affine_tool $ADD_ON_TRANSMAT \
        -oitk $ADD_ON_TRANS
    fi

    $BINDIR/ml_affine \
      $GSTEMPDIR/GSUTemplate/gshoot/template_1/template/iter_2/template_1_gshoot_seg.nii.gz \
      $REGSEG \
      $SUBJDATADIR/${id}_${side}_to_MSTInitTemp_mlaffine.txt

    if [[ $side == "right" ]]; then

      $C3DPATH/c3d_affine_tool \
        $SUBJDATADIR/${id}_${side}_to_MSTInitTemp_mlaffine.txt \
        -oitk $TMPDIR/${id}_${side}_to_MSTInitTemp_mlaffine_itk.txt

      $BINDIR/ComposeMultiTransform 3 \
        $TMPDIR/${id}_${side}_to_MSTInitTemp_mlaffine_combined_itk.txt \
        -R $TMPDIR/${id}_${side}_to_MSTInitTemp_mlaffine_itk.txt \
        $TMPDIR/${id}_${side}_to_MSTInitTemp_mlaffine_itk.txt \
        $ADD_ON_TRANS

       $C3DPATH/c3d_affine_tool \
         -itk $TMPDIR/${id}_${side}_to_MSTInitTemp_mlaffine_combined_itk.txt \
         -o $SUBJDATADIR/${id}_${side}_to_MSTInitTemp_mlaffine.txt
    fi
  fi

  # transform the labels and the SRT1
  $BINDIR/greedy -d 3  \
    -rf $GSTEMPDIR/GSUTemplate/gshoot/template_1/template/iter_2/template_1_gshoot_seg.nii.gz \
    $RM \
    -r $SUBJDATADIR/${id}_${side}_to_MSTInitTemp_mlaffine.txt

  # generate segmentation in the template space
  $C3DPATH/c3d $(for sub in ${LABEL_IDS_FIT[*]}; do echo $SUBJDATADIR/${id}_${side}_${sub}.nii.gz; done) \
    -vote -type ushort -o $SUBJDATADIR/${id}_${side}_seg.nii.gz

  fi

  ######################################
  # register to all the atlases and compute similarity
  idside=${id}_${side}
  ATLASES=$(cat $ATLASLIST) 
  SUBJREGATLASDIR=$WORKDIR/RegToAtlases
 
  for idside_fix in $ATLASES; do
    OUTTMPDIR=$SUBJREGATLASDIR/$idside/${idside}_to_${idside_fix}
    mkdir -p $OUTTMPDIR

    # perform registration
    if [[ -f $OUTTMPDIR/${idside}_to_${idside_fix}_sim.txt ]]; then

      echo "Seg file exists."

    else

      # file names
      #TMPDIR=$OUTTMPDIR
      MAT_MOMENTS=$TMPDIR/${idside}_to_${idside_fix}_moment.mat
      MAT_AFFINE=$TMPDIR/${idside}_to_${idside_fix}_affine.mat
      WARP=$TMPDIR/${idside}_to_${idside_fix}_warp.nii.gz

      # greedy command
      CMD=""
      for sub in ${LABEL_IDS_FIT[*]}; do
        CMD="$CMD -w 1 -i $GSTEMPDIR/data//${idside_fix}_${sub}.nii.gz $SUBJDATADIR/${idside}_${sub}.nii.gz "
      done

      # Perform moments of intertia matching between the two masks
      #$BINDIR/greedy -d 3  \
      #  $CMD \
      #  -moments \
      #  -o $MAT_MOMENTS

      # Perform affine matching between the two masks
      $BINDIR/greedy -d 3  \
       $CMD \
       -a -ia-identity \
       -n 100x100 \
       -o $MAT_AFFINE

      #-a -ia $MAT_MOMENTS \

      # Run greedy between these two images
      $C3DPATH/c3d $GSTEMPDIR/data/${idside_fix}_seg.nii.gz -dup \
        $SUBJDATADIR/${idside}_seg.nii.gz \
        -int 0 -reslice-identity \
        -add -binarize -dilate 1 10x10x10vox \
        -o $TMPDIR/${idside}_to_${idside_fix}_mask.nii.gz
      MASK="-gm $TMPDIR/${idside}_to_${idside_fix}_mask.nii.gz"

     $BINDIR/greedy -d 3  \
        $CMD \
        -it $MAT_AFFINE \
        $MASK \
        -n 50x50x20x0 \
        -s 2.0mm 0.1mm -e 0.5 \
        -o $WARP

      # Reslice the segmentations from raw space
      RM=""
      for sub in ${LABEL_IDS_FIT[*]}; do
        RM="$RM -rm $SUBJDATADIR/${idside}_${sub}.nii.gz $TMPDIR/${idside}_to_${idside_fix}_reslice_${sub}.nii.gz"
      done

      $BINDIR/greedy -d 3  \
        -rf $GSTEMPDIR/data/${idside_fix}_seg.nii.gz \
        $RM \
        -r $WARP $MAT_AFFINE

      # Create seg
      $C3DPATH/c3d $(for sub in ${LABEL_IDS_FIT[*]}; do echo $TMPDIR/${idside}_to_${idside_fix}_reslice_${sub}.nii.gz; done) \
        -vote -type ushort \
        -o $TMPDIR/${idside}_to_${idside_fix}_reslice_seg.nii.gz

      # measure similarity
      OVL=$($C3DPATH/c3d $TMPDIR/${idside}_to_${idside_fix}_reslice_seg.nii.gz \
        $MSTEVALREPSTR \
        $GSTEMPDIR/data/${idside_fix}_seg.nii.gz \
        $MSTEVALREPSTR \
        -label-overlap \
        | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' )
      echo "$OVL" > \
        $OUTTMPDIR/${idside}_to_${idside_fix}_sim.txt

    fi

  done
}

######################################################
function TempMembership()
{
  side=$1
  idside=${id}_${side}
  SUBJDATADIR=$WORKDIR/data
  SUBJREGATLASDIR=$WORKDIR/RegToAtlases
  SUBJMEMBERSHIPDIR=$WORKDIR/membership
  ATLASES=$(cat $ATLASLIST)
  ATLASESGRPS=$(cat $GTGROUPS)
  mkdir -p $SUBJMEMBERSHIPDIR

  # get similarity matrix  
  PWSIM=""

  for idside_fix in $ATLASES; do

    PWREGDIR=$SUBJREGATLASDIR/${idside}/${idside}_to_${idside_fix}
    SIMFILE=$PWREGDIR/${idside}_to_${idside_fix}_sim.txt
    if [[ ! -f $SIMFILE ]]; then
      echo "error: $SIMFILE does not exist" >> $LOGTXT
      if [[ $PWSIM == "" ]]; then
        PWSIM="nan"
      else
        PWSIM="$PWSIM,nan"
      fi
    else
      if [[ $PWSIM == "" ]]; then
        PWSIM="$(cat $SIMFILE)"
      else
        PWSIM="$PWSIM,$(cat $SIMFILE)"
      fi
    fi

  done

  echo $PWSIM >> $SUBJMEMBERSHIPDIR/adj_${side}.csv
  echo "$id,$side,${id}_${side}" >>  $SUBJMEMBERSHIPDIR/all_info_${side}.txt

  # compute group membership
  source $CODEDIR/matlabcode/group_membership_MultiAndUnifiTemp/for_testing/run_group_membership_MultiAndUnifiTemp.sh \
     $MATLAB_RT \
     $SUBJMEMBERSHIPDIR/adj_${side}.csv \
     $GTGROUPS \
     6 \
     $SUBJMEMBERSHIPDIR/autogroup_${side}.txt \
     $SUBJMEMBERSHIPDIR/UTautogroup_${side}.txt
  unset LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=/home/customlib/lib/:/home/customlib/icclibs/
}

######################################################
function RegToInitTemp()
{
  side=$1
  SUBJMEMBERSHIPDIR=$WORKDIR/membership
  SUBJREGINITTEMPDIR=$WORKDIR/RegToInitTemp
  ATLASES=($(cat $ATLASLIST))
  mkdir -p $SUBJREGINITTEMPDIR

  ######################
  # step 1: register to the most similar subject
  Line=($(cat $SUBJMEMBERSHIPDIR/autogroup_${side}.txt | head -n 1 | tail -n 1))
  echo ${Line[*]}
  group=${Line[0]}
  idx_inter=${Line[1]}
  id_inter=${ATLASES[$((idx_inter-1))]}
  echo $id_inter
  echo $idx_inter

  # file names
  INITREGDIR=$SUBJREGINITTEMPDIR/init
  mkdir -p $INITREGDIR
  PREFIX=$INITREGDIR/${id}_${side}_to_${id_inter}
  MAT_MOMENTS=${PREFIX}_moment.mat
  MAT_AFFINE=${PREFIX}_affine.mat
  WARP=${PREFIX}_warp.nii.gz

  if [[ ! -f $MAT_AFFINE || ! -f $WARP ]]; then

    # greedy command
    CMD=""
    for sub in ${LABEL_IDS_FIT[*]}; do
      CMD="$CMD -w 1 -i $GSTEMPDIR/data/${id_inter}_${sub}.nii.gz $SUBJDATADIR/${idside}_${sub}.nii.gz "
    done

    # Perform moments of intertia matching between the two masks
    #$BINDIR/greedy -d 3  \
    #  $CMD \
    #  -moments \
    #  -o $MAT_MOMENTS

    $BINDIR/ml_affine \
      $GSTEMPDIR/data/${id_inter}_seg.nii.gz \
      $SUBJDATADIR/${idside}_seg.nii.gz \
      $MAT_MOMENTS

    # Perform affine matching between the two masks
    $BINDIR/greedy -d 3  \
      $CMD \
     -a -ia $MAT_MOMENTS \
     -n 100x100 \
     -o $MAT_AFFINE

    # Run greedy between these two images
    $C3DPATH/c3d $GSTEMPDIR/data/${id_inter}_seg.nii.gz -dup \
      $SUBJDATADIR/${idside}_seg.nii.gz \
      -int 0 -reslice-identity \
      -add -binarize -dilate 1 10x10x10vox \
      -o $TMPDIR/${idside}_to_${id_inter}_mask.nii.gz

    $BINDIR/greedy -d 3 \
      $CMD \
      -it $MAT_AFFINE \
      -gm $TMPDIR/${idside}_to_${id_inter}_mask.nii.gz \
      -n 50x40x20 -float \
      -s 2vox 1vox -e 0.5 \
      -o $WARP

  fi

  # load and update the warp chian
  TMPWARPCHAIN=($(cat $GSTEMPDIR/GSTemplate/MST/registration/template_${group}/$id_inter/final/chain_unwarp_to_final.txt))
  WARPCHAIN=$(for ((i=0;i<${#TMPWARPCHAIN[*]};i++)); do echo "$GSTEMPDIR/${TMPWARPCHAIN[i]}"; done)
  WARPCHAIN="$WARPCHAIN $WARP $MAT_AFFINE"
  echo $WARPCHAIN > $INITREGDIR/chain_unwarp_to_final_${side}.txt
  rm -rf $INITREGDIR/chain_unwarp_to_final.txt

  # Reslice the segmentations from raw space
  RM=""
  for sub in ${KINDS_FIT[*]}; do
    RM="$RM -rm $SUBJDATADIR/${idside}_${sub}.nii.gz $INITREGDIR/${idside}_to_MSTRoot_${group}_reslice_${sub}.nii.gz"
  done

  $BINDIR/greedy -d 3  \
    -rf $GSTEMPDIR/GSTemplate/MST/template/template_${group}/template_${group}_seg.nii.gz \
    $RM \
    -r $WARPCHAIN

  # Create seg
  $C3DPATH/c3d $(for sub in ${LABEL_IDS_FIT[*]}; do echo $INITREGDIR/${idside}_to_MSTRoot_${group}_reslice_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $INITREGDIR/${idside}_to_MSTRoot_${group}_reslice_seg.nii.gz

  ######################
  # step 2: additional registration to initial template
  INITTEMPREGDIR=$SUBJREGINITTEMPDIR/inittemp
  mkdir -p $INITTEMPREGDIR
  INITWORKDIR=$GSTEMPDIR/GSTemplate/InitTemp/template_${group}/work/

  CMD=""
  for sub in ${LABEL_IDS_FIT[*]}; do
    CMD="$CMD -w 1 -i $INITWORKDIR/template_${group}_${sub}.nii.gz $INITREGDIR/${idside}_to_MSTRoot_${group}_reslice_${sub}.nii.gz "
  done

  # Run greedy between these two images
  WARP=$INITTEMPREGDIR/${idside}_totempWarp.nii.gz

  # Run greedy between these two images
  if [[ ! -f $WARP ]]; then

    $C3DPATH/c3d $INITWORKDIR/template_${group}_seg.nii.gz -dup \
      $INITREGDIR/${idside}_to_MSTRoot_${group}_reslice_seg.nii.gz \
      -int 0 -reslice-identity \
      -add -binarize -dilate 1 10x10x10vox \
      -o $TMPDIR/${idside}_totemp_mask.nii.gz

    $BINDIR/greedy -d 3  \
      $CMD \
      -n 120x120x40 \
      -gm $TMPDIR/${idside}_totemp_mask.nii.gz \
      -s 0.6mm 0.1mm \
      -e 0.5 \
      -o $WARP

      # 120x120x40

  fi

  # Reslice the segmentations from raw space
  RM=""
  for sub in ${KINDS_FIT[*]}; do
    RM="$RM -rm $SUBJDATADIR/${idside}_${sub}.nii.gz $INITTEMPREGDIR/${idside}_totemp_${group}_reslice_${sub}.nii.gz"
  done

  $BINDIR/greedy -d 3 \
    -rf $INITWORKDIR/template_${group}_BKG.nii.gz \
    $RM \
    -r $WARP $WARPCHAIN

  $C3DPATH/c3d $(for sub in ${LABEL_IDS_FIT[*]}; do echo $INITTEMPREGDIR/${idside}_totemp_${group}_reslice_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $INITTEMPREGDIR/${idside}_totemp_${group}_reslice_seg.nii.gz

  WARPCHAIN="$WARP $WARPCHAIN"
  echo $WARPCHAIN > $INITTEMPREGDIR/chain_unwarp_to_final_${side}.txt
}

######################################################
function GeoShooting()
{
  side=$1
  idside=${id}_${side}
  SUBJDATADIR=$WORKDIR/data
  SUBJREGATLASDIR=$WORKDIR/RegToAtlases
  SUBJMEMBERSHIPDIR=$WORKDIR/membership
  SUBJREGINITTEMPDIR=$WORKDIR/RegToInitTemp
  SUBJGEOSHOOTDIR=$WORKDIR/GeoShoot/$side
  mkdir -p $SUBJGEOSHOOTDIR
  Line=($(cat $SUBJMEMBERSHIPDIR/autogroup_${side}.txt | head -n 1 | tail -n 1))
  group=${Line[0]}

  # The path to the landmarks
  INITLANDMARKS=$GSTEMPDIR/GSTemplate/InitTemp/template_${group}/work/iter_04/template_${group}_MRGcombined_sampled.vtk
  LANDMARKS=$GSTEMPDIR/GSTemplate/gshoot/template_${group}/shape_avg/iter_1/shavg_landmarks.vtk

  # Reference space (root node in cm-rep space)
  REFSPACE=$SUBJGEOSHOOTDIR/refspace.nii.gz
  if [ ! -f $REFSPACE ]; then
    PERC=25
    PAD=60
    $C3DPATH/c3d $SUBJDATADIR/${idside}_seg_orig.nii.gz \
      -resample ${PERC}x${PERC}x${PERC}% \
      -pad ${PAD}x${PAD}x${PAD}vox ${PAD}x${PAD}x${PAD}vox 0 \
      -o $TMPDIR/${idside}_seg_orig.nii.gz

    $BINDIR/AverageImages 3 \
      $REFSPACE 0 \
      $GSTEMPDIR/GSTemplate/InitTemp/template_${group}/work/iter_04/template_${group}_seg.nii.gz \
      $GSTEMPDIR/GSTemplate/gshoot/template_${group}/refspace_${group}.nii.gz \
      $TMPDIR/${idside}_seg_orig.nii.gz

    $C3DPATH/c3d $REFSPACE \
      -thresh 0.0001 inf 1 0 \
      -trim 20vox \
      -resample-mm 0.4x0.4x0.4mm \
      -o $REFSPACE
  fi

  # Result meshes
  TARGET=$SUBJGEOSHOOTDIR/shooting_target_native.vtk
  LM_PROCRUSTES=$SUBJGEOSHOOTDIR/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$SUBJGEOSHOOTDIR/target_to_root_procrustes.mat
  SHOOTING_WARP=$SUBJGEOSHOOTDIR/shooting_warp.nii.gz

  # Mesh containing the momenta
  MOMENTA=$SUBJGEOSHOOTDIR/shooting_momenta.vtk

  # Target-related stuff in the WORK directory that is only done in the first iter
  # Get the warp chain from file
  WARPCHAIN=$(cat $SUBJREGINITTEMPDIR/inittemp/chain_unwarp_to_final_${side}.txt)

  # Apply the warp chain to the landmark mesh in template space, creating
  # the target locations for the geodesic shooting registration
  # -- this code works when the WARPCHAIN is empty (MST root)
  if [ ! -f $TARGET ]; then
    $BINDIR/greedy -d 3 \
      -rf $REFSPACE \
      -rs $INITLANDMARKS $TARGET \
      -r $WARPCHAIN \
         $SUBJDATADIR/${id}_${side}_to_MSTInitTemp_mlaffine.txt
  fi

  # Landmarks in reference space
  rm -f $SUBJGEOSHOOTDIR/landmarks.vtk
  ln -sf $LANDMARKS $SUBJGEOSHOOTDIR/landmarks.vtk

  # Bring the target mesh back near the root mesh using procrustes alignment
  if [ ! -f $LM_PROCRUSTES_MAT ]; then
    $BINDIR/vtkprocrustes $TARGET $LANDMARKS $LM_PROCRUSTES_MAT
  fi

  # Apply procrustes to the landmarks.
  #warpmesh $TARGET $LM_PROCRUSTES $LM_PROCRUSTES_MAT
  if [ ! -f $LM_PROCRUSTES ]; then
    $BINDIR/greedy -d 3 \
      -rf $REFSPACE \
      -rs $TARGET $LM_PROCRUSTES \
      -r $LM_PROCRUSTES_MAT
  fi

  # Perform geodesic shooting between the procrustes landmarks and the
  # warped landmarks - this is going to allow us to interpolate the correspondence
  # found by the MST to the rest of the images
  if [[ ! -f $MOMENTA ]]; then

    $BINDIR/lmshoot -d 3 \
      -m $LANDMARKS $LM_PROCRUSTES \
      -s $GSSIGMA -l $GSWEIGHT -n $GSTP -i $GSITER 0 \
      -o $MOMENTA

  fi
  
  #Warp template to native space image
  VOLITERDIR=$GSTEMPDIR/GSTemplate/gshoot/template_${group}/template/iter_2
  if [ ! -f $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_seg.nii.gz ]; then

  TRANS=""
  TRANS1=""
  for ((i=0;i<${#MESHES[*]};i++)); do
    TRANS="$TRANS -M $VOLITERDIR/template_${group}_gshoot_${MESHES[i]}.vtk $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_${MESHES[i]}.vtk"
    TRANS1="$TRANS1 -rs $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_${MESHES[i]}.vtk $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_${MESHES[i]}.vtk"
  done
  $BINDIR/lmtowarp -d 3 -n $GSTP -s $GSSIGMA \
    -m $MOMENTA \
    $TRANS

  $BINDIR/greedy -d 3  \
    -rf $REFSPACE \
    $TRANS1 \
    -r $LM_PROCRUSTES_MAT,-1

  for ((i=0;i<${#MESHES[*]};i++)); do

    $BINDIR/mesh2img -f -vtk $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_${MESHES[i]}.vtk \
      -a 0.3 0.3 0.3 4 \
      $TMPDIR/template_${group}_to_${idside}_GSShoot_${MESHES[i]}.nii.gz

    $C3DPATH/c3d $SUBJDATADIR/${idside}_lfseg_heur_dividedCS.nii.gz \
      $TMPDIR/template_${group}_to_${idside}_GSShoot_${MESHES[i]}.nii.gz \
      -int 0 -reslice-identity \
      -o $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_${MESHES[i]}.nii.gz

  done

  $C3DPATH/c3d $(for ((i=1;i<${#LABEL_IDS_MESH[*]};i++)); do echo $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_${LABEL_IDS_MESH[i]}.nii.gz; done) \
    -foreach -sdt -scale -1 -endfor \
    -vote -shift 1 \
    $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_NOBKG.nii.gz \
    -multiply -type ushort \
    -o $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_seg.nii.gz
  fi
}

######################################################
function PostSteps()
{
  side=$1
  EvalFit $side
  MeasureThickness $side
  CleanUp $side
}

######################################################
function RegToUT()
{
  side=$1
  idside=${id}_${side}
  SUBJDATADIR=$WORKDIR/data
  SUBJREGATLASDIR=$WORKDIR/RegToAtlases
  SUBJMEMBERSHIPDIR=$WORKDIR/membership
  SUBJREGUTTEMPDIR=$WORKDIR/RegToUT
  ATLASES=($(cat $ATLASLIST))
  mkdir -p $SUBJREGUTTEMPDIR

  ######################
  # step 1: register to the most similar subject
  Line=($(cat $SUBJMEMBERSHIPDIR/UTautogroup_${side}.txt | head -n 1 | tail -n 1))
  echo ${Line[*]}
  group=${Line[0]}
  idx_inter=${Line[1]}
  id_inter=${ATLASES[$((idx_inter-1))]}
  echo $id_inter
  echo $idx_inter

  # file names
  INITREGDIR=$SUBJREGUTTEMPDIR/init
  mkdir -p $INITREGDIR
  PREFIX=$INITREGDIR/${id}_${side}_to_${id_inter}
  #MAT_INIT=$GSTEMPDIR/GSUTemplate/gshoot/template_1/$id_inter/iter_final/target_to_root_procrustes.mat
  MAT_MOMENTS=${PREFIX}_moment.mat
  MAT_AFFINE=${PREFIX}_affine.mat
  WARP=${PREFIX}_warp.nii.gz

  if [[ ! -f $MAT_AFFINE || ! -f $WARP ]]; then

    # get the segmentation smaller
    $C3DPATH/c3d $GSTEMPDIR/data/${id_inter}_seg_orig.nii.gz \
      -trim 10x10x10vox \
      -o $TMPDIR/${id_inter}_seg_orig.nii.gz

    # greedy command
    CMD=""
    for sub in ${LABEL_IDS_FIT[*]}; do
      $C3DPATH/c3d $TMPDIR/${id_inter}_seg_orig.nii.gz \
        $GSTEMPDIR/data/${id_inter}_${sub}_orig.nii.gz \
        -reslice-identity \
        -o $TMPDIR/${id_inter}_${sub}_orig.nii.gz
      CMD="$CMD -w 1 -i $TMPDIR/${id_inter}_${sub}_orig.nii.gz $SUBJDATADIR/${idside}_${sub}.nii.gz "
    done

    # Perform moments of intertia matching between the two masks
    #greedy -d 3  \
    #  $CMD \
    #  -moments \
    #  -o $MAT_MOMENTS

    $BINDIR/ml_affine \
      $TMPDIR/${id_inter}_seg_orig.nii.gz \
      $SUBJDATADIR/${idside}_seg.nii.gz \
      $MAT_MOMENTS

    # Perform affine matching between the two masks
    $BINDIR/greedy -d 3  \
      $CMD \
     -a -ia $MAT_MOMENTS \
     -n 100x100 \
     -o $MAT_AFFINE

    $BINDIR/greedy -d 3  \
      $CMD \
      -it $MAT_AFFINE \
      -n 50x40x20 -float \
      -s 2vox 1vox -e 0.5 \
      -o $WARP

  fi

  # load and update the warp chian
  #TMPWARPCHAIN=($(cat $GSTEMPDIR/GSUTemplate/MST/registration/template_1/$id_inter/final/chain_unwarp_to_final.txt))
  #WARPCHAIN=$(for ((i=0;i<${#TMPWARPCHAIN[*]};i++)); do echo "$GSTEMPDIR/${TMPWARPCHAIN[i]}"; done)
  WARPCHAIN="$GSTEMPDIR/GSUTemplate/gshoot/template_1/$id_inter/iter_final/shooting_warp.nii.gz $GSTEMPDIR/GSUTemplate/gshoot/template_1/$id_inter/iter_final/target_to_root_procrustes.mat,-1"
  WARPCHAIN="$WARPCHAIN $WARP $MAT_AFFINE"
  echo $WARPCHAIN > $INITREGDIR/chain_unwarp_to_final_${side}.txt
  #rm -rf $INITREGDIR/chain_unwarp_to_final.txt

  # Reslice the segmentations from raw space
  RM=""
  for sub in ${KINDS_FIT[*]}; do
    RM="$RM -rm $SUBJDATADIR/${idside}_${sub}.nii.gz $INITREGDIR/${idside}_to_MSTRoot_reslice_${sub}.nii.gz"
  done

  $BINDIR/greedy -d 3   \
    -rf $GSTEMPDIR/GSUTemplate/gshoot/template_1/template/iter_2/template_1_gshoot_seg.nii.gz \
    $RM \
    -r $WARPCHAIN

  # Create seg
  $C3DPATH/c3d $(for sub in ${LABEL_IDS_FIT[*]}; do echo $INITREGDIR/${idside}_to_MSTRoot_reslice_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $INITREGDIR/${idside}_to_MSTRoot_reslice_seg.nii.gz

  ######################
  # step 2: additional registration to initial template
  INITTEMPREGDIR=$SUBJREGUTTEMPDIR/UTemp
  mkdir -p $INITTEMPREGDIR
  INITWORKDIR=$GSTEMPDIR/GSUTemplate/gshoot/template_1/template/iter_2/

  CMD=""
  for sub in ${LABEL_IDS_FIT[*]}; do
    CMD="$CMD -w 1 -i $INITWORKDIR/template_1_gshoot_${sub}.nii.gz $INITREGDIR/${idside}_to_MSTRoot_reslice_${sub}.nii.gz "
  done

  # Run greedy between these two images
  WARP=$INITTEMPREGDIR/${idside}_totempWarp.nii.gz

  # Run greedy between these two images
  if [[ ! -f $WARP ]]; then

    $C3DPATH/c3d $INITWORKDIR/template_1_gshoot_seg.nii.gz -dup \
      $INITREGDIR/${idside}_to_MSTRoot_reslice_seg.nii.gz \
      -int 0 -reslice-identity \
      -add -binarize -dilate 1 10x10x10vox \
      -o $TMPDIR/${idside}_totemp_mask.nii.gz

    $BINDIR/greedy -d 3  \
      $CMD \
      -n 120x120x40 \
      -gm $TMPDIR/${idside}_totemp_mask.nii.gz \
      -s 0.6mm 0.3mm \
      -e 0.5 \
      -o $WARP

      # 120x120x40

  fi

  # Reslice the segmentations from raw space
  RM=""
  for sub in ${KINDS_FIT[*]}; do
    RM="$RM -rm $SUBJDATADIR/${idside}_${sub}.nii.gz $INITTEMPREGDIR/${idside}_totemp_reslice_${sub}.nii.gz"
  done

  $BINDIR/greedy -d 3  \
    -rf $INITWORKDIR/template_1_gshoot_BKG.nii.gz \
    $RM \
    -r $WARP $WARPCHAIN

  $C3DPATH/c3d $(for sub in ${LABEL_IDS_FIT[*]}; do echo $INITTEMPREGDIR/${idside}_totemp_reslice_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $INITTEMPREGDIR/${idside}_totemp_reslice_seg.nii.gz

  WARPCHAIN="$WARP $WARPCHAIN"
  echo $WARPCHAIN > $INITTEMPREGDIR/chain_unwarp_to_final_${side}.txt
}

######################################################
function GeoShootingUT()
{
  side=$1
  idside=${id}_${side}
  SUBJDATADIR=$WORKDIR/data
  SUBJREGATLASDIR=$WORKDIR/RegToAtlases
  SUBJMEMBERSHIPDIR=$WORKDIR/membership
  SUBJREGUTTEMPDIR=$WORKDIR/RegToUT
  SUBJGEOSHOOTDIR=$WORKDIR/GeoShootUT/$side
  mkdir -p $SUBJGEOSHOOTDIR
  group=1
  
  # The path to the landmarks
  LANDMARKS=$GSTEMPDIR/GSUTemplate/gshoot/template_${group}/template/iter_2/template_${group}_gshoot_MRGcombined_sampled.vtk

  # Reference space (root node in cm-rep space)
  REFSPACE=$SUBJGEOSHOOTDIR/refspace.nii.gz
  if [ ! -f $REFSPACE ]; then
    PERC=25
    PAD=60
    $C3DPATH/c3d $SUBJDATADIR/${idside}_seg_orig.nii.gz \
      -resample ${PERC}x${PERC}x${PERC}% \
      -pad ${PAD}x${PAD}x${PAD}vox ${PAD}x${PAD}x${PAD}vox 0 \
      -o $TMPDIR/${idside}_seg_orig.nii.gz

    $BINDIR/AverageImages 3 \
      $REFSPACE 0 \
      $GSTEMPDIR/GSUTemplate/InitTemp/template_${group}/work/iter_04/template_${group}_seg.nii.gz \
      $GSTEMPDIR/GSUTemplate/gshoot/template_${group}/refspace_${group}.nii.gz \
      $TMPDIR/${idside}_seg_orig.nii.gz

    $C3DPATH/c3d $REFSPACE \
      -thresh 0.0001 inf 1 0 \
      -trim 20vox \
      -resample-mm 0.4x0.4x0.4mm \
      -o $REFSPACE
  fi

  # Result meshes
  TARGET=$SUBJGEOSHOOTDIR/shooting_target_native.vtk
  LM_PROCRUSTES=$SUBJGEOSHOOTDIR/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$SUBJGEOSHOOTDIR/target_to_root_procrustes.mat
  SHOOTING_WARP=$SUBJGEOSHOOTDIR/shooting_warp.nii.gz

  # Mesh containing the momenta
  MOMENTA=$SUBJGEOSHOOTDIR/shooting_momenta.vtk

  # Target-related stuff in the WORK directory that is only done in the first iter
  # Get the warp chain from file
  WARPCHAIN=$(cat $SUBJREGUTTEMPDIR/UTemp/chain_unwarp_to_final_${side}.txt)

  # Apply the warp chain to the landmark mesh in template space, creating
  # the target locations for the geodesic shooting registration
  # -- this code works when the WARPCHAIN is empty (MST root)
  if [ ! -f $TARGET ]; then
    $BINDIR/greedy -d 3  \
      -rf $REFSPACE \
      -rs $LANDMARKS $TARGET \
      -r $WARPCHAIN \
         $SUBJDATADIR/${id}_${side}_to_MSTInitTemp_mlaffine.txt
  fi

  # Landmarks in reference space
  rm -f $SUBJGEOSHOOTDIR/landmarks.vtk
  ln -sf $LANDMARKS $SUBJGEOSHOOTDIR/landmarks.vtk

  # Bring the target mesh back near the root mesh using procrustes alignment
  if [ ! -f $LM_PROCRUSTES_MAT ]; then
    $BINDIR/vtkprocrustes $TARGET $LANDMARKS $LM_PROCRUSTES_MAT
  fi

  # Apply procrustes to the landmarks.
  #warpmesh $TARGET $LM_PROCRUSTES $LM_PROCRUSTES_MAT
  if [ ! -f $LM_PROCRUSTES ]; then
    $BINDIR/greedy -d 3  \
      -rf $REFSPACE \
      -rs $TARGET $LM_PROCRUSTES \
      -r $LM_PROCRUSTES_MAT
  fi

  # Perform geodesic shooting between the procrustes landmarks and the
  # warped landmarks - this is going to allow us to interpolate the correspondence
  # found by the MST to the rest of the images
  if [[ ! -f $MOMENTA ]]; then

    $BINDIR/lmshoot -d 3 \
      -m $LANDMARKS $LM_PROCRUSTES \
      -s $GSSIGMA -l $GSWEIGHT -n $GSTP -i $GSITER 0 \
      -o $MOMENTA

  fi

  #Warp template to native space image
  VOLITERDIR=$GSTEMPDIR/GSUTemplate/gshoot/template_${group}/template/iter_2
  if [ ! -f $SUBJGEOSHOOTDIR/template_to_${idside}_GSShoot_seg.nii.gz ]; then

  TRANS=""
  TRANS1=""
  for ((i=0;i<${#MESHES[*]};i++)); do
    TRANS="$TRANS -M $VOLITERDIR/template_${group}_gshoot_${MESHES[i]}.vtk $SUBJGEOSHOOTDIR/template_to_${idside}_GSShoot_${MESHES[i]}.vtk"
    TRANS1="$TRANS1 -rs $SUBJGEOSHOOTDIR/template_to_${idside}_GSShoot_${MESHES[i]}.vtk $SUBJGEOSHOOTDIR/template_to_${idside}_GSShoot_${MESHES[i]}.vtk"
  done
  $BINDIR/lmtowarp -d 3 -n $GSTP -s $GSSIGMA \
    -m $MOMENTA \
    $TRANS

  $BINDIR/greedy -d 3  \
    -rf $REFSPACE \
    $TRANS1 \
    -r $LM_PROCRUSTES_MAT,-1

  for ((i=0;i<${#MESHES[*]};i++)); do

    $BINDIR/mesh2img -f -vtk $SUBJGEOSHOOTDIR/template_to_${idside}_GSShoot_${MESHES[i]}.vtk \
      -a 0.3 0.3 0.3 4 \
      $TMPDIR/template_to_${idside}_GSShoot_${MESHES[i]}.nii.gz

    $C3DPATH/c3d $SUBJDATADIR/${idside}_lfseg_heur_dividedCS.nii.gz \
      $TMPDIR/template_to_${idside}_GSShoot_${MESHES[i]}.nii.gz \
      -int 0 -reslice-identity \
      -o $SUBJGEOSHOOTDIR/template_to_${idside}_GSShoot_${MESHES[i]}.nii.gz

  done

  $C3DPATH/c3d $(for ((i=1;i<${#LABEL_IDS_MESH[*]};i++)); do echo $SUBJGEOSHOOTDIR/template_to_${idside}_GSShoot_${LABEL_IDS_MESH[i]}.nii.gz; done) \
    -foreach -sdt -scale -1 -endfor \
    -vote -shift 1 \
    $SUBJGEOSHOOTDIR/template_to_${idside}_GSShoot_NOBKG.nii.gz \
    -multiply -type ushort \
    -o $SUBJGEOSHOOTDIR/template_to_${idside}_GSShoot_seg.nii.gz
  fi
}

######################################################
function EvalFit()
{
  side=$1
  idside=${id}_${side}
  SUBJDATADIR=$WORKDIR/data
  SUBJREGATLASDIR=$WORKDIR/RegToAtlases
  SUBJMEMBERSHIPDIR=$WORKDIR/membership
  SUBJREGINITTEMPDIR=$WORKDIR/RegToInitTemp
  SUBJGEOSHOOTDIR=$WORKDIR/GeoShoot/$side
  SUBJUTGEOSHOOTDIR=$WORKDIR/GeoShootUT/$side
  SUBJEVALDIR=$WORKDIR/evaluation
  mkdir -p $SUBJEVALDIR
  Line=($(cat $SUBJMEMBERSHIPDIR/autogroup_${side}.txt | head -n 1 | tail -n 1))
  group=${Line[0]}

  ASHSSEG=$SUBJDATADIR/${idside}_seg_orig.nii.gz
  FITTED=$SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_seg.nii.gz
  if [[ -f $FITTED && ! -f $SUBJEVALDIR/${idside}_GSShootASHS_overlap.csv ]]; then
  do_pair $ASHSSEG $FITTED
  echo "$id,$side,$group$FULLOVL" > \
    $SUBJEVALDIR/${idside}_GSShootASHS_overlap.csv
  fi

  UTFITTED=$SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_seg.nii.gz
  if [[ -f $UTFITTED && ! -f $SUBJEVALDIR/${idside}_UTGSShootASHS_overlap.csv ]]; then
  do_pair $ASHSSEG $UTFITTED
  echo "$id,$side,$group$FULLOVL" > \
    $SUBJEVALDIR/${idside}_UTGSShootASHS_overlap.csv
  fi
}

function do_pair()
{
  # Get a pair of segmentations
  seg_a=$1
  seg_b=$2

  # out dice file
  #out_dice_file=$3

  # Iterate over all relevant labels
  FULLOVL=""
  for ((i=0; i<${#EVALLABELS[*]}; i++)); do

    # Do the analysis on full-size meshes
    REPRULE=$(for lab in ${RANGES[i]}; do echo $lab 99; done)

    # Extract the binary images and compute overlap
    $C3DPATH/c3d \
      $seg_a -dup $seg_b -int 0 -reslice-identity \
      -foreach -replace $REPRULE -thresh 99 99 1 0 -endfor \
      -overlap 1 | tee $TMPDIR/ovl.txt

    # Get the full-extent overlap
    OVL=$(cat $TMPDIR/ovl.txt | grep OVL | awk -F '[ ,]+' '{print $6}')

    #echo $id ${LABELS[i]} full $OVL $DIST >> $out_file
    FULLOVL="${FULLOVL},${OVL}"

  done
}

######################################################
function MeasureThickness()
{
  side=$1
  idside=${id}_${side}
  SUBJDATADIR=$WORKDIR/data
  SUBJREGATLASDIR=$WORKDIR/RegToAtlases
  SUBJMEMBERSHIPDIR=$WORKDIR/membership
  SUBJREGINITTEMPDIR=$WORKDIR/RegToInitTemp
  SUBJGEOSHOOTDIR=$WORKDIR/GeoShoot/$side
  SUBJUTGEOSHOOTDIR=$WORKDIR/GeoShootUT/$side
  SUBJEVALDIR=$WORKDIR/evaluation
  Line=($(cat $SUBJMEMBERSHIPDIR/autogroup_${side}.txt | head -n 1 | tail -n 1))
  group=${Line[0]}

  if [[ ! -f $SUBJGEOSHOOTDIR/${idside}_mean_thickness.vtk && ! -f $SUBJGEOSHOOTDIR/${idside}_median_thickness.vtk ]]; then

  # Extract the thickness of the subfield
  if [[ ! -f $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap.vtk && -f $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG.vtk ]]; then
    $BINDIR/cmrep_vskel -Q $BINDIR/qvoronoi \
      -T $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap.vtk \
      -p $thick_p -e $thick_e \
      $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG.vtk \
      $TMPDIR/template_${group}_to_${idside}_GSShoot_MRG_skel.vtk
  fi

  # sample thickness map on probability maps
  if [[ ! -f $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap_withlabel.vtk && -f $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG.vtk ]]; then
    TMPMESHES=""
    LABELDEF=(${MESHESDEF[${MESHMRGIDX}]})
    for ((i=0;i<${#LABEL_IDS_FIT[*]};i++)); do
      if [[ ${LABELDEF[$i]} == "1" ]]; then
        sub=${LABEL_IDS_FIT[$i]}

        # binarize seg
        $C3DPATH/c3d $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_seg.nii.gz \
          -thresh $i $i 1 0 \
          -smooth 1vox \
          -o $TMPDIR/template_${group}_to_${idside}_GSShoot_${sub}_smooth.nii.gz

        $BINDIR/mesh_image_sample \
          $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap.vtk \
          $TMPDIR/template_${group}_to_${idside}_GSShoot_${sub}_smooth.nii.gz \
          $TMPDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap_${sub}.vtk \
          PROB
        TMPMESHES="$TMPMESHES $TMPDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap_${sub}.vtk"
      fi
    done

    # merge prob arrays
    $BINDIR/mesh_merge_arrays \
      -r $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap.vtk \
      $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap_withlabel.vtk \
      PROB $TMPMESHES
  fi

  # run matlab script to generate label and mean thickness for each label
  if [[ ! -f $SUBJGEOSHOOTDIR/${idside}_mean_thickness.vtk && ! -f $SUBJGEOSHOOTDIR/${idside}_median_thickness.vtk && -f $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG.vtk ]]; then
  source $CODEDIR/matlabcode/compute_meshlabel_meanthickness/for_testing/run_compute_meshlabel_meanthickness.sh \
      $MATLAB_RT \
      $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap_withlabel.vtk \
      $id \
      $side \
      $SUBJGEOSHOOTDIR/${idside}_mean_thickness.vtk \
      $SUBJGEOSHOOTDIR/${idside}_median_thickness.vtk
  unset LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=/home/customlib/lib:/home/customlib/icclibs
  fi

  fi

  # unified template
   
  if [[ ! -f $SUBJUTGEOSHOOTDIR/${idside}_mean_thickness.vtk && ! -f $SUBJUTGEOSHOOTDIR/${idside}_median_thickness.vtk ]]; then
	  
  # Extract the thickness of the subfield
  if [[ ! -f $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG_thickmap.vtk && -f $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG.vtk  ]]; then
    $BINDIR/cmrep_vskel -Q $BINDIR/qvoronoi \
      -T $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG_thickmap.vtk \
      -p $thick_p -e $thick_e \
      $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG.vtk \
      $TMPDIR/template_to_${idside}_GSShoot_MRG_skel.vtk
  fi

  # sample thickness map on probability maps
  if [[ ! -f $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG_thickmap_withlabel.vtk && -f $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG.vtk ]]; then
    TMPMESHES=""
    LABELDEF=(${MESHESDEF[${MESHMRGIDX}]})
    for ((i=0;i<${#LABEL_IDS_FIT[*]};i++)); do
      if [[ ${LABELDEF[$i]} == "1" ]]; then
        sub=${LABEL_IDS_FIT[$i]}

        # binarize seg
        $C3DPATH/c3d $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_seg.nii.gz \
          -thresh $i $i 1 0 \
          -smooth 1vox \
          -o $TMPDIR/template_to_${idside}_GSShoot_${sub}_smooth.nii.gz

        $BINDIR/mesh_image_sample \
          $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG_thickmap.vtk \
          $TMPDIR/template_to_${idside}_GSShoot_${sub}_smooth.nii.gz \
          $TMPDIR/template_to_${idside}_GSShoot_MRG_thickmap_${sub}.vtk \
          PROB
        TMPMESHES="$TMPMESHES $TMPDIR/template_to_${idside}_GSShoot_MRG_thickmap_${sub}.vtk"
      fi
    done

    # merge prob arrays
    $BINDIR/mesh_merge_arrays \
      -r $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG_thickmap.vtk \
      $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG_thickmap_withlabel.vtk \
      PROB $TMPMESHES
  fi

  # run matlab script to generate label and mean thickness for each label
  if [[ ! -f $SUBJUTGEOSHOOTDIR/${idside}_mean_thickness.vtk && ! -f $SUBJUTGEOSHOOTDIR/${idside}_median_thickness.vtk && -f $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG.vtk ]]; then
  source $CODEDIR/matlabcode/compute_meshlabel_meanthickness/for_testing/run_compute_meshlabel_meanthickness.sh \
    $MATLAB_RT \
    $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG_thickmap_withlabel.vtk \
    $id \
    $side \
    $SUBJUTGEOSHOOTDIR/${idside}_mean_thickness.vtk \
    $SUBJUTGEOSHOOTDIR/${idside}_median_thickness.vtk
  unset LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=/home/customlib/lib/:/home/customlib/icclibs/
  fi

  fi
}

######################################################
function CleanUp()
{
  side=$1
  idside=${id}_${side}

  # remove output
  OUTCSV=$OUTDIR/${idside}_thickness.csv
  rm -f $OUTCSV

  # header
  header="ID,SIDE,TempType,GROUP"
  for type in MeanThk MedianThk FitQuality; do
  for sub in ${MESH_LABEL_FIT[*]}; do
    header="$header,${sub}_${type}"
  done
  done
  header="$header,All_FitQuality"
  echo $header > $OUTCSV
  
  SUBJDATADIR=$WORKDIR/data
  SUBJREGATLASDIR=$WORKDIR/RegToAtlases
  SUBJMEMBERSHIPDIR=$WORKDIR/membership
  SUBJREGINITTEMPDIR=$WORKDIR/RegToInitTemp
  SUBJGEOSHOOTDIR=$WORKDIR/GeoShoot/$side
  SUBJUTGEOSHOOTDIR=$WORKDIR/GeoShootUT/$side
  SUBJEVALDIR=$WORKDIR/evaluation

  # get group 
  Line=($(cat $SUBJMEMBERSHIPDIR/autogroup_${side}.txt | head -n 1 | tail -n 1))
  group=${Line[0]}
  outstr="$id,$side,MultiTemp,$group"

  set +e
  # extract mean thickness
  if [[ -f $SUBJGEOSHOOTDIR/${idside}_mean_thickness.vtk && -f $SUBJGEOSHOOTDIR/${idside}_median_thickness.vtk ]]; then
  begin=3
  end=$(($begin+${#MESH_LABEL_FIT[*]}))
  end=$((end-1))
  outstr="$outstr,$(cat -A $SUBJGEOSHOOTDIR/${idside}_mean_thickness.vtk | cut -d , -f ${begin}-${end} | sed -e 's/\$//g')"
    
  # extract median thickness
  outstr="$outstr,$(cat -A $SUBJGEOSHOOTDIR/${idside}_median_thickness.vtk | cut -d , -f ${begin}-${end} | sed -e 's/\$//g')"

  outstr="$outstr,$(cat -A $SUBJEVALDIR/${idside}_GSShootASHS_overlap.csv | cut -d , -f ${EVALINX} | sed -e 's/\$//g')"
    
  echo $outstr >> $OUTCSV
  fi

  # unified template
  outstr="$id,$side,UnifiedTemp,1"
    
  # extract mean thickness
  if [[ -f $SUBJUTGEOSHOOTDIR/${idside}_mean_thickness.vtk && -f $SUBJUTGEOSHOOTDIR/${idside}_median_thickness.vtk ]]; then
  outstr="$outstr,$(cat -A $SUBJUTGEOSHOOTDIR/${idside}_mean_thickness.vtk | cut -d , -f ${begin}-${end} | sed -e 's/\$//g')"

  # extract median thickness
  outstr="$outstr,$(cat -A $SUBJUTGEOSHOOTDIR/${idside}_median_thickness.vtk | cut -d , -f ${begin}-${end} | sed -e 's/\$//g')"

  outstr="$outstr,$(cat -A $SUBJEVALDIR/${idside}_UTGSShootASHS_overlap.csv | cut -d , -f ${EVALINX} | sed -e 's/\$//g')"

  echo $outstr >> $OUTCSV
  fi
  set -e

  # keep the fitted mesh and momenta
  if [[ -f $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap_withlabel.vtk ]]; then
    cp $SUBJGEOSHOOTDIR/template_${group}_to_${idside}_GSShoot_MRG_thickmap_withlabel.vtk \
       $OUTDIR/${idside}_template_${group}_fitted_mesh.vtk
  fi
  if [[ -f $SUBJGEOSHOOTDIR/shooting_momenta.vtk ]]; then
    cp $SUBJGEOSHOOTDIR/shooting_momenta.vtk \
       $OUTDIR/${idside}_template_${group}_momenta.vtk
  fi
  
  # Unified template
  if [[ -f $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG_thickmap_withlabel.vtk ]]; then 
    cp $SUBJUTGEOSHOOTDIR/template_to_${idside}_GSShoot_MRG_thickmap_withlabel.vtk \
       $OUTDIR/${idside}_UT_template_fitted_mesh.vtk
  fi
  if [[ -f $SUBJUTGEOSHOOTDIR/shooting_momenta.vtk ]]; then
    cp $SUBJUTGEOSHOOTDIR/shooting_momenta.vtk \
       $OUTDIR/${idside}_UT_template_momenta.vtk
  fi

}

#################################################
function DeleteFile()
{
  # delete intermediate file if specify
  if [[ $DELETETMP == "1" ]]; then
    rm -rf $WORKDIR
  fi
}

main
