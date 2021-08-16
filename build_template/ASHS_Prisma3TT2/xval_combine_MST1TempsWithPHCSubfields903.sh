#!/bin/bash
#$ -S /bin/bash
set -x -e

##############################################################################
# INIT
#ANTSPATH=/home/avants/bin/ants
ANTSPATH=/data/picsl/longxie/pkg/bin/antsbin/bin/
#ANTSPATH=/data/picsl/pauly/bin/ants/
NEWANTSDIR=/share/apps/ANTs/2014-06-23/build/bin/
C3DPATH=/data/picsl/longxie/pkg/c3d_tool/bin
PKGDIR=/data/picsl/longxie/pkg/
MATLAB_BIN=/share/apps/matlab/R2017a/bin/matlab
GREEDY_BIN=/data/picsl/pauly/bin/
export PATH=$ANTSPATH:$C3DPATH:$PATH
RCODEDIR=/home/longxie/ASHS_T1/scripts

# TMPDIR
#if [[ ! $TMPDIR ]]; then
#  TMPDIR=/tmp
#fi

TMPDIR=/data/picsl/longxie/tmp/$(mktemp -u  | awk -F "/" '{print $NF}')
rm -rf $TMPDIR
mkdir -p $TMPDIR

# Directories
ROOT=/home/longxie/ASHS_PrismaT2/
WORKDIR=$ROOT/thickness
ASHSRUNDIR=/data/picsl/pauly/wolk/abc_prisma/exp03/atlas
ASHSTEMPDIR=/data/picsl/pauly/wolk/abc_prisma/exp03/final/template
#MATLABCODEDIR=$WORKDIR/matlabcode/
GTGROUPDIR=$ROOT/analysis_input/
GTSEGDIR=$ROOT/Prisma3T_atlas
SUBJ_TXT=$ROOT/analysis_input/subj_xval_vertical.txt
IDENTITYFN=$ROOT/analysis_input/identity_itk.txt
LMTOWARPDIR=/home/longxie/LMTemplateMatching/PointSetUtilities/build_release

##############################################################################
# Parameters needs to be specify
# Experiment number
expid=903
expdir=$WORKDIR/exp/exp${expid}
DUMPDIR=${expdir}/dump
mkdir -p ${expdir}/dump

########################################
# 0. copy data
CLEANUPDIR=$expdir/cleanup

LABEL_IDS_ALL=(BKG        CA      DG  TAIL SUB ERC  BA35 BA36 PHC  aCS  pCS)
LABEL_MRG_ALL=("0 6 7 13" "1 2 4" "3" "5"  "8" "9"   "10" "11" "12" "15" "14")
# Labels to get segmentation
LABEL_IDS=(BKG        CA      DG  TAIL SUB ERC  BA35 BA36 PHC  aCS  pCS)
LABEL_MRG=("0 6 7 13" "1 2 4" "3" "5"  "8" "9"   "10" "11" "12" "15" "14")
LABEL_NEW=(0          1       2   3    4   5    6    7    8    9    10)
MESH_LABEL=(CA DG TAIL SUB ERC  BA35 BA36 PHC)
KINDS="tse mprage ${LABEL_IDS[*]}"
DATADIR=${expdir}/data

########################################
# 1 pairwise registration and clustering
# 1.1 ANTs parameters
PWDIR=${expdir}/pairwise/

# similarity type
SIM_TYPE="PRCCS_seg_dice"
SIM_DIR=$expdir/PWsim/sim_${SIM_TYPE}

# 1.2 path directory
groups=(1)
GSTEMPDIR=$expdir/GSTemplate
ALLPATHDIR=$GSTEMPDIR/MST/paths

########################################
# 2. MST registration
MSTREGDIR=$GSTEMPDIR/MST/registration
MSTTEMPDIR=$GSTEMPDIR/MST/template

# 2.2 generate template
LABEL_FG=(CA DG TAIL SUB ERC BA35 BA36 PHC CS)
MESHES=(CA DG TAIL SUB ERC BA35 BA36 PHC CS MRG NOBKG)
#           BKG HIPPO ERC BA35 BA36 PHC aCS pCS
#           BKG CA DG TAIL SUB ERC BA35 BA36 PHC aCS pCS
MESHESDEF=("-1   1 -1 -1   -1  -1  -1   -1   -1  -1  -1" \
           "-1  -1  1 -1   -1  -1  -1   -1   -1  -1  -1" \
           "-1  -1 -1  1   -1  -1  -1   -1   -1  -1  -1" \
           "-1  -1 -1 -1    1  -1  -1   -1   -1  -1  -1" \
           "-1  -1 -1 -1   -1   1  -1   -1   -1  -1  -1" \
           "-1  -1 -1 -1   -1  -1   1   -1   -1  -1  -1" \
           "-1  -1 -1 -1   -1  -1  -1    1   -1  -1  -1" \
           "-1  -1 -1 -1   -1  -1  -1   -1    1  -1  -1" \
           "-1  -1 -1 -1   -1  -1  -1   -1   -1   1   1" \
           "-1   1 -1  1    1   1   1    1    1  -1  -1" \
           "-1   1  1  1    1   1   1    1    1   1   1")
#            CA  DG  TAIL SUB ERC BA35 BA36 PHC CS  MRG  NOBKG
SAMPLEPOINT=(500 300 200  400 250 300  400  350 350 3050 10)
COMBINEMESHES=(CA TAIL SUB ERC BA35 BA36 PHC)

# 2.3 additional registration
INITITER=5
INITTEMPDIR=$GSTEMPDIR/InitTemp

# 2.3. shooting
SHOOTDIR=$GSTEMPDIR/gshoot
GSHOOT_NITER=3

# 2.5
LABELWARPDIR=$expdir/labelwarp

# 2.6 evaluation
EVALDIR=$expdir/evaluation
EVALLABELS=(CA DG TAL SUB ERC BA35 BA36 PRC   PHC aCS pCS HIPP        EXPHIPP   ALL)
RANGES=(    1  2  3   4   5   6    7    "6 7" 8   9   10   "1 2 3 4"  "5 6 7 8" "1 2 3 4 5 6 7 8")

##############################################################################
function main()
{
  #reset_dir

  #######################################
  # 0. preparation
  #PreprocessLabels
  #copy_data

  #######################################
  # 1 perform registration
  # 1.1 pairwise registration
  #pairwise

  # 1.2 find paths between all subjects and the root for both template 1 and template 2
  #PrintAdjacency

  #######################################
  # 2. build template using MST for template 1 and template 2
  # 2.1 registration through MST
  #RegisterMST

  # 2.2 generate initial template and generate landmarks
  #GenerateInitTemp

  # 2.3 generate a better initial template by performing iterative averaging
  #RegisterInitTemp

  # 2.3. refine template using geodesic shooting
  # perform shape correction
  #ShootingCorrection
  #ShootCorrectedMesh
  #QC_shooting

  # 2.4. Perform shooting one more time using the new template
  #FinalShooting

  # 2.5. compute mean and mode
  ShapeStats
  MakeMovies

  # 2.6 warp templates back to subject space
  WarpTempToSubj

  # 2.7 evaluate overlap
  EvalGSTempFit
  #CopyShootMeshes
}

##########################################################
function PreprocessLabels()
{
  mkdir -p $CLEANUPDIR

  # Load id number
  IDs=($(cat $SUBJ_TXT))

  # Submit job to copy data
  PREFIX=PL${expid}
  for side in left right; do
    for ((i=0;i<${#IDs[*]};i++)); do

      id=${IDs[i]}

      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${id}_${side}" \
           $0 PreprocessLabels_sub $id $side
      sleep 0.1

     done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function PreprocessLabels_sub()
{
  id=$1
  side=$2

  in=$GTSEGDIR/orig_${id}/${id}_${side}_finalLW.nii.gz
  out5=$TMPDIR/${id}_${side}_sf_full_cleanup.nii.gz
  out=$CLEANUPDIR/${id}_${side}_dividedCS.nii.gz

  $C3DPATH/c3d $in -as S \
    -thresh 9 11 1 0 -sdt -scale -1 -popas A \
    -push S -thresh 12 12 1 0 -sdt -scale -1 -popas P \
    -push P -push A -vote -popas AP \
    -push S -thresh 14 14 1 0 -push AP -multiply \
    -push S -add \
    -o $out
  ln -sf $in $CLEANUPDIR/
}

###########################################################
# Copy data
function copy_data()
{
  # perform initial affine registration between left and flip right temps
  if [ ! -f $DATADIR/tempreg/right_to_left_init_affine.txt ]; then

    mkdir -p $DATADIR/tempreg
    ln -sf $ASHSTEMPDIR/refspace_meanseg_left.nii.gz \
       $DATADIR/tempreg/refspace_meanseg_left.nii.gz
    ln -sf $ASHSTEMPDIR/refspace_meanseg_right.nii.gz \
       $DATADIR/tempreg/refspace_meanseg_right.nii.gz
    greedy -d 3 \
      -i $DATADIR/tempreg/refspace_meanseg_left.nii.gz \
         $DATADIR/tempreg/refspace_meanseg_right.nii.gz \
      -moments 2 \
      -o $DATADIR/tempreg/right_to_left_init_affine.txt
    c3d_affine_tool $DATADIR/tempreg/right_to_left_init_affine.txt \
      -oitk $DATADIR/tempreg/right_to_left_init_affine.txt

  fi

  # Get the data
  mkdir -p $DATADIR

  # Load id number
  IDs=($(cat $SUBJ_TXT))

  # Submit job to copy data
  PREFIX=CP${expid}
  for side in left right; do
    for ((i=0;i<${#IDs[*]};i++)); do

      id=${IDs[i]}
      fn=$id #$(ls $ASHSRUNDIR | grep $id)

      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${id}_${side}" \
           $0 copy_subject $id $side
      sleep 0.1

     done
  done

    # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function copy_subject()
{
  id=$1
  side=$2
  fn=$id #$(ls $ASHSRUNDIR | grep $id)
  SUBJASHSRUNDIR=$ASHSRUNDIR/orig_${id}

  # ASHS segmentation
  SEG=$CLEANUPDIR/${id}_${side}_dividedCS.nii.gz
  if [[ ! -f $SEG ]]; then
    echo "$SEG does not exist."
    exit
  fi

  # Link the subfield images
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
    ln -sf $SUBJASHSRUNDIR/tse_to_chunktemp_${side}.nii.gz \
           $DATADIR/${id}_${side}_tse.nii.gz
    ln -sf $SUBJASHSRUNDIR/tse_to_chunktemp_${side}.nii.gz \
           $DATADIR/${id}_${side}_mprage.nii.gz
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
    WarpImageMultiTransform 3  \
      $SUBJASHSRUNDIR/tse_to_chunktemp_${side}.nii.gz \
      $DATADIR/${id}_${side}_tse.nii.gz \
      -R $SUBJASHSRUNDIR/tse_to_chunktemp_left.nii.gz \
      $ADD_ON_TRANS
    WarpImageMultiTransform 3  \
      $SUBJASHSRUNDIR/tse_to_chunktemp_${side}.nii.gz \
      $DATADIR/${id}_${side}_mprage.nii.gz \
      -R $SUBJASHSRUNDIR/tse_to_chunktemp_left.nii.gz \
      $ADD_ON_TRANS
  fi

  #if [ ! -f $SUBJASHSRUNDIR/affine_t1_to_template/t1_to_template_affine_itk.txt ]; then
  #  c3d_affine_tool $SUBJASHSRUNDIR/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
  #    -oitk $SUBJASHSRUNDIR/ants_t1_to_temp/t1_to_template_affine_itk.txt
    #c3d_affine_tool $SUBJASHSRUNDIR/affine_t1_to_template/t1_to_template_affine.mat \
    #  -oitk $SUBJASHSRUNDIR/affine_t1_to_template/t1_to_template_affine_itk.txt
  #fi

  #if [ ! -f $SUBJASHSRUNDIR/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt ]; then
  #  c3d_affine_tool $SUBJASHSRUNDIR/flirt_t2_to_t1/flirt_t2_to_t1.mat \
  #    -oitk $SUBJASHSRUNDIR/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt
  #fi

  # Generate a binary image for the label
  for ((i=0;i<${#LABEL_IDS_ALL[*]};i++)); do

    c3d $SEG -replace $(for k in ${LABEL_MRG_ALL[i]}; do echo $k 999; done) -thresh 999 999 1 0 \
      -o $DATADIR/${id}_${side}_${LABEL_IDS_ALL[i]}_orig.nii.gz \

    greedy -d 3 \
      -rf $DATADIR/${id}_${side}_tse.nii.gz \
      -rm $DATADIR/${id}_${side}_${LABEL_IDS_ALL[i]}_orig.nii.gz \
          $DATADIR/${id}_${side}_${LABEL_IDS_ALL[i]}.nii.gz \
      -r $ADD_ON_TRANS \
         $SUBJASHSRUNDIR/ants_t1_to_temp/greedy_t1_to_template_warp.nii.gz \
         $SUBJASHSRUNDIR/affine_t1_to_template/t1_to_template_affine.mat \
         $SUBJASHSRUNDIR/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat

    #WarpImageMultiTransform 3  \
    #  $DATADIR/${id}_${side}_${LABEL_IDS_ALL[i]}_orig.nii.gz \
    #  $DATADIR/${id}_${side}_${LABEL_IDS_ALL[i]}.nii.gz \
    #  -R $DATADIR/${id}_${side}_tse.nii.gz \
    #  $ADD_ON_TRANS \
    #  $SUBJASHSRUNDIR/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz \
    #  $SUBJASHSRUNDIR/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
    #  $SUBJASHSRUNDIR/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

    #if [[ $side == "right" ]]; then
    #  c3d $DATADIR/${id}_${side}_${LABEL_IDS_ALL[i]}_orig.nii.gz \
    #    -flip x  \
    #    -o $DATADIR/${id}_${side}_${LABEL_IDS_ALL[i]}_orig.nii.gz
    #fi

  done

  # Vote in the ASHS template subject space to create segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $DATADIR/${id}_${side}_${sub}.nii.gz; done) \
    -vote -type ushort -o $DATADIR/${id}_${side}_seg.nii.gz

  c3d $(for sub in ${LABEL_IDS[*]}; do echo $DATADIR/${id}_${side}_${sub}_orig.nii.gz; done) \
    -vote -type ushort -o $DATADIR/${id}_${side}_seg_orig.nii.gz
}

############################################################
# Pairwise registration
function pairwise()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=PW${expid}
  for side in left right; do

    for id_fix in $IDS; do

      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${id_fix}_${side}" $0 \
           pairwise_sub $id_fix $side
      sleep 0.1

    done

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function pairwise_sub()
{
  id_fix=$1
  side_fix=$2
  idside_fix=${id_fix}_${side_fix}
  IDS=$(cat $SUBJ_TXT)
  fn_fix=$id_fix #$(ls $ASHSRUNDIR | grep $id_fix)
  fnside_fix=${fn_fix}_${side_fix}

  # preprocess for similarity measure
  fn_seg_fix=$TMPDIR/${fnside_fix}_seg_tmp.nii.gz
  c3d $DATADIR/${fnside_fix}_seg.nii.gz \
    -replace 1 0 2 0 3 0 4 0 5 0 8 0 10 0 \
    -o $fn_seg_fix

  for id_mov in $IDS; do
  for side_mov in left right; do

    idside_mov=${id_mov}_${side_mov}
     if [[ $idside_mov != $idside_fix ]]; then

      fn_mov=$id_mov #$(ls $ASHSRUNDIR | grep $id_mov)
      fnside_mov=${fn_mov}_${side_mov}
      OUTDIR=$PWDIR/${fnside_fix}/${fnside_mov}_to_${fnside_fix}
      #TMPDIR=$OUTDIR
      mkdir -p $PWDIR/${fnside_fix}

      # perform registration
      if [[ -f $OUTDIR/${fnside_mov}_to_${fnside_fix}_sim.txt ]]; then

        echo "Seg file exists."

      else

        mkdir -p $OUTDIR

        # file names
        MAT_MOMENTS=$TMPDIR/${fnside_mov}_to_${fnside_fix}_moment.mat
        MAT_AFFINE=$TMPDIR/${fnside_mov}_to_${fnside_fix}_affine.mat
        WARP=$TMPDIR/${fnside_mov}_to_${fnside_fix}_warp.nii.gz

        # greedy command
        CMD=""
        for sub in ${LABEL_IDS[*]}; do
          CMD="$CMD -w 1 -i $DATADIR/${fnside_fix}_${sub}.nii.gz $DATADIR/${fnside_mov}_${sub}.nii.gz "
        done

        # Perform moments of intertia matching between the two masks
        greedy -d 3 -threads $NSLOTS \
          $CMD \
          -moments \
          -o $MAT_MOMENTS

        # Perform affine matching between the two masks
        greedy -d 3 -threads $NSLOTS \
          $CMD \
         -a -ia $MAT_MOMENTS \
         -n 100x100 \
         -o $MAT_AFFINE

        # Run greedy between these two images
        c3d $DATADIR/${fnside_fix}_seg.nii.gz -dup \
          $DATADIR/${fnside_mov}_seg.nii.gz \
          -int 0 -reslice-identity \
          -add -binarize -dilate 1 10x10x10vox \
          -o $TMPDIR/${fnside_mov}_to_${fnside_fix}_mask.nii.gz
        MASK="-gm $TMPDIR/${fnside_mov}_to_${fnside_fix}_mask.nii.gz"

        greedy -d 3 -threads $NSLOTS \
          $CMD \
          -it $MAT_AFFINE \
          $MASK \
          -n 50x50x20x0 \
          -s 2.0mm 0.1mm -e 0.5 \
          -o $WARP

        # Reslice the segmentations from raw space
        RM=""
        for sub in ${KINDS[*]}; do
          RM="$RM -rm $DATADIR/${fnside_mov}_${sub}.nii.gz $TMPDIR/${fnside_mov}_to_${fnside_fix}_reslice_${sub}.nii.gz"
        done

        greedy -d 3 \
          -rf $DATADIR/${fnside_fix}_seg.nii.gz \
          $RM \
          -r $WARP $MAT_AFFINE

        # Create seg
        c3d $(for sub in ${LABEL_IDS[*]}; do echo $TMPDIR/${fnside_mov}_to_${fnside_fix}_reslice_${sub}.nii.gz; done) \
          -vote -type ushort \
          -o $TMPDIR/${fnside_mov}_to_${fnside_fix}_reslice_seg.nii.gz

        # measure similarity
        OVL=$(c3d $TMPDIR/${fnside_mov}_to_${fnside_fix}_reslice_seg.nii.gz \
          -replace 1 0 2 0 3 0 4 0 5 0 8 0 10 0 \
          $fn_seg_fix -label-overlap \
          | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' )
        echo "$OVL" > \
          $OUTDIR/${fnside_mov}_to_${fnside_fix}_sim.txt

      fi
    fi

  done
  done
}

##########################################################
function PrintAdjacency()
{
  IDS=($(cat $SUBJ_TXT))
  GRPLEFT=($(cat $GTGROUPDIR/group_xval_2group_left.txt))
  GRPRIGHT=($(cat $GTGROUPDIR/group_xval_2group_right.txt))
  rm -rf $ALLPATHDIR
  mkdir -p $ALLPATHDIR

  # create adjacency matrix
  if [ ! -f $ALLPATHDIR/adj.txt ]; then
  idx_fix=0
  idx_fix_1=0
  idx_fix_2=0
  for side_fix in left right; do
  for ((i_fix=0;i_fix<${#IDS[*]};i_fix++)); do
  #for id_fix in $IDS; do

    id_fix=${IDS[$i_fix]}
    idside_fix=${id_fix}_${side_fix}
    fn_fix=$id_fix #fn_fix=$(ls $ASHSRUNDIR | grep $id_fix)
    fnside_fix=${fn_fix}_${side_fix}
    if [[ $side_fix == "left" ]]; then
      grp_fix=${GRPLEFT[$i_fix]}
    else
      grp_fix=${GRPRIGHT[$i_fix]}
    fi

    idx_mov=0
    idx_mov_1=0
    for side_mov in left right; do
    for ((i_mov=0;i_mov<${#IDS[*]};i_mov++)); do
    #for id_mov in $IDS; do

      id_mov=${IDS[$i_mov]}
      idside_mov=${id_mov}_${side_mov}
      fn_mov=$id_mov #fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)
      fnside_mov=${fn_mov}_${side_mov}
      if [[ $side_mov == "left" ]]; then
        grp_mov=${GRPLEFT[$i_mov]}
      else
        grp_mov=${GRPRIGHT[$i_mov]}
      fi

      PWREGDIR=$PWDIR/${fnside_fix}/${fnside_mov}_to_${fnside_fix}
      SIMFILE=$PWREGDIR/${fnside_mov}_to_${fnside_fix}_sim.txt
      if [[ $idside_fix != $idside_mov ]]; then

        if [[ ! -f $SIMFILE ]]; then
          echo "error: $SIMFILE does not exist"
          exit
        fi

        echo "$((idx_fix+1)) $((idx_mov+1)) $(cat $SIMFILE)" \
          >> $ALLPATHDIR/adj.txt

        if [[ $grp_fix == $grp_mov ]]; then
          if [[ $grp_fix == "1" ]]; then
            echo "$((idx_fix_1+1)) $((idx_mov_1+1)) $(cat $SIMFILE)" \
              >> $ALLPATHDIR/adj_${grp_mov}.txt
          else
            echo "$((idx_fix_2+1)) $((idx_mov_1+1)) $(cat $SIMFILE)" \
              >> $ALLPATHDIR/adj_${grp_mov}.txt
          fi
        fi
      fi

      idx_mov=$(($idx_mov+1))
      if [[ $grp_fix == $grp_mov ]]; then
        idx_mov_1=$(($idx_mov_1+1))
      fi

    done
    done

    idx_fix=$(($idx_fix+1))
    if [[ $grp_fix == "1" ]]; then
      idx_fix_1=$(($idx_fix_1+1))
    else
      idx_fix_2=$(($idx_fix_2+1))
    fi

    echo "$fnside_fix" >> $ALLPATHDIR/IDSide.txt
    echo "$fnside_fix" >> $ALLPATHDIR/IDSide_${grp_fix}.txt
    echo "$id_fix $side_fix" >>  $ALLPATHDIR/ID_and_Side_${grp_fix}.txt
  done
  done
  fi

  /share/apps/R/R-3.1.1/bin/Rscript $RCODEDIR/compute_mst_onegraph.R $ALLPATHDIR/adj.txt \
    > $ALLPATHDIR/mst_paths.txt
  /share/apps/R/R-3.1.1/bin/Rscript $RCODEDIR/compute_mst.R $ALLPATHDIR/adj_1.txt \
    > $ALLPATHDIR/mst_paths_1.txt
  /share/apps/R/R-3.1.1/bin/Rscript $RCODEDIR/compute_mst.R $ALLPATHDIR/adj_2.txt \
    > $ALLPATHDIR/mst_paths_2.txt
}

############################################################
function RegisterMST()
{
  # id
  #IDS=$(cat $SUBJ_TXT)
  #IDS="DW124"

  # submit job for each subject
  PREFIX=RMST${expid}
  for group in ${groups[*]}; do
    IDS=($(cat $ALLPATHDIR/IDSide.txt))
    for ((idx=1;idx<=${#IDS[*]};idx++)); do

      id=${IDS[$((idx-1))]}
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=6.1G,s_vmem=6G \
           -N "${PREFIX}_${group}_${idx}_${id}" $0 \
           RegisterMST_sub $group $idx $id
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function RegisterMST_sub()
{
  group=$1
  idxmov=$2
  #IDS=$(cat $SUBJ_TXT)
  #fnside_mov=$(cat $ALLPATHDIR/IDSide.txt | head -n $idxmov | tail -n 1)
  fnside_mov=$3

  # get path
  GPATH=($(cat $ALLPATHDIR/mst_paths.txt | head -n $idxmov | tail -n 1))

  # We start with the warp chain as an empty list
  WARPCHAIN=""

  # The id in whose space the moving image is - at the beginning this is the
  # same as the subject's ID
  IDMOV=$fnside_mov

  # Loop. The first element of GPATH is the image itself
  for ((j=1;j<${#GPATH[*]};++j)); do

    # The next image in the path
    idxref=${GPATH[j]}
    IDREF=$(cat $ALLPATHDIR/IDSide.txt | head -n $idxref | tail -n 1)

    # Create directory for this registration
    WORK=$MSTREGDIR/template_${group}/$fnside_mov/step_${j}_${IDMOV}_to_${IDREF}
    mkdir -p $WORK

    # file names
    PREFIX=$WORK/${IDMOV}_to_${IDREF}
    MAT_MOMENTS=${PREFIX}_moment.mat
    MAT_AFFINE=${PREFIX}_affine.mat
    WARP=${PREFIX}_warp.nii.gz

    if [[ ! -f $MAT_AFFINE || ! -f $WARP ]]; then

    # greedy command
    CMD=""
    for sub in ${LABEL_IDS[*]}; do
      w=1
      #if [[ $sub == "PHC" ]]; then
      #  w="0.1"
      #fi
      CMD="$CMD -w ${w} -i $DATADIR/${IDREF}_${sub}.nii.gz $DATADIR/${IDMOV}_${sub}.nii.gz "
    done

    # Perform moments of intertia matching between the two masks
    #greedy -d 3 -threads $NSLOTS \
    #  $CMD \
    #  -moments \
    #  -o $MAT_MOMENTS

    # Perform affine matching between the two masks
    #greedy -d 3 -threads $NSLOTS \
    #  $CMD \
    # -a -ia $MAT_MOMENTS \
    # -n 100x100 \
    # -o $MAT_AFFINE

    # Run greedy between these two images
    c3d $DATADIR/${IDREF}_seg.nii.gz -dup \
      $DATADIR/${IDMOV}_seg.nii.gz \
      -int 0 -reslice-identity \
      -add -binarize -dilate 1 10x10x10vox \
      -o $TMPDIR/${IDMOV}_to_${IDREF}_mask.nii.gz

    greedy -d 3 -threads $NSLOTS \
      $CMD \
      -ia-identity \
      -gm $TMPDIR/${IDMOV}_to_${IDREF}_mask.nii.gz \
      -n 50x40x20 -float \
      -s 2vox 1vox -e 0.5 \
      -o $WARP

      #-it $MAT_AFFINE \

      #-s 0.6mm 0.1mm

    fi

    # Update the warp chain
    #WARPCHAIN="$WARP $MAT_AFFINE $WARPCHAIN"
    WARPCHAIN="$WARP $WARPCHAIN"
    IDMOV=${IDREF}

    # link the reference image
    for sub in ${KINDS[*]} seg; do
      ln -sf $DATADIR/${IDREF}_${sub}.nii.gz \
        $WORK/${IDREF}_${sub}.nii.gz
    done

    # Reslice the segmentations from raw space
    RM=""
    for sub in ${KINDS[*]}; do
      RM="$RM -rm $DATADIR/${fnside_mov}_${sub}.nii.gz $WORK/${fnside_mov}_to_${IDREF}_reslice_${sub}.nii.gz"
    done

    greedy -d 3 \
      -rf $DATADIR/${IDREF}_seg.nii.gz \
      $RM \
      -r $WARPCHAIN

    # Create seg
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $WORK/${fnside_mov}_to_${IDREF}_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $WORK/${fnside_mov}_to_${IDREF}_reslice_seg.nii.gz

  done

  # make link to the final directory
  FINALWORK=$MSTREGDIR/template_${group}/$fnside_mov/final
  mkdir -p $FINALWORK
  if [[ ${#GPATH[*]} == 1 ]]; then
    # for root subject, just link original data
    for sub in ${KINDS[*]} seg; do
      ln -sf $DATADIR/${fnside_mov}_${sub}.nii.gz \
        $FINALWORK/${fnside_mov}_to_MSTRoot_reslice_${sub}.nii.gz
    done
  else
    for sub in ${KINDS[*]} seg; do
      ln -sf $WORK/${fnside_mov}_to_${IDREF}_reslice_${sub}.nii.gz \
        $FINALWORK/${fnside_mov}_to_MSTRoot_reslice_${sub}.nii.gz
    done
  fi
  echo $WARPCHAIN > $FINALWORK/chain_unwarp_to_final.txt
}

##########################################################
function GenerateInitTemp()
{
  # create mean image for different labels
  PREFIX=GIT${expid}
  for group in ${groups[*]}; do
    mkdir -p $MSTTEMPDIR/template_${group}
    for sub in ${LABEL_IDS[*]}; do

      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${group}_${sub}" $0 \
           GenerateTemp_sub $group $sub INIT
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # generate seg
  for group in ${groups[*]}; do

  c3d $(for sub in ${LABEL_IDS[*]}; do echo $MSTTEMPDIR/template_${group}/template_${group}_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $MSTTEMPDIR/template_${group}/template_${group}_seg.nii.gz

  # generate meshes and subsample
  for ((i=0;i<${#MESHES[*]};i++)); do

    LABELDEF=(${MESHESDEF[i]})
    echo $LABELDEF

    # generate meshes
    c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$MSTTEMPDIR/template_${group}/template_${group}_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
      -accum -add -endaccum \
      -o $MSTTEMPDIR/template_${group}/template_${group}_${MESHES[i]}_vtk.nii.gz

    vtklevelset $MSTTEMPDIR/template_${group}/template_${group}_${MESHES[i]}_vtk.nii.gz \
      $MSTTEMPDIR/template_${group}/template_${group}_${MESHES[i]}.stl 0.0

    # perform subsampling
    SAM=$MSTTEMPDIR/template_${group}/template_${group}_${MESHES[i]}_sampled.ply
    /home/longxie/pkg/mesh_poisson_sample \
      $MSTTEMPDIR/template_${group}/template_${group}_${MESHES[i]}.stl \
      $SAM ${SAMPLEPOINT[i]}

    #/data/picsl-build/pauly/vcg/gcc64rel/mesh_poisson_sample

  done

  # combine and convert into vtk mesh
  TEMPLATE=$MSTTEMPDIR/template_${group}/template_${group}_MRGcombined_sampled.vtk
  NV=0
  for ((i=0;i<${#COMBINEMESHES[*]};i++)); do
    SAM=$MSTTEMPDIR/template_${group}/template_${group}_${COMBINEMESHES[i]}_sampled.ply
    NV=$(($NV+$(cat $SAM | grep 'element vertex' | awk '{print $3}')))
  done

  # Write the header of the VTK file
  echo "# vtk DataFile Version 4.0" > $TEMPLATE
  echo "vtk output" >> $TEMPLATE
  echo "ASCII" >> $TEMPLATE
  echo "DATASET POLYDATA" >> $TEMPLATE
  echo "POINTS ${NV} float" >> $TEMPLATE

  # Write the point coordinates
  for ((i=0;i<${#COMBINEMESHES[*]};i++)); do
    SAM=$MSTTEMPDIR/template_${group}/template_${group}_${COMBINEMESHES[i]}_sampled.ply
    NV=$(cat $SAM | grep 'element vertex' | awk '{print $3}')
    grep -A $NV end_header $SAM | tail -n $NV >> $TEMPLATE
  done

  done
}

function GenerateTemp_sub()
{
  group=$1
  sub=$2
  MODE=$3

  if [[ $MODE == "INIT" ]]; then
    c3d $MSTREGDIR/template_${group}/*/final/*reslice_${sub}.nii.gz \
      -mean -o $MSTTEMPDIR/template_${group}/template_${group}_${sub}.nii.gz
  elif [[ $MODE == "INITREG" ]]; then
    c3d $INITTEMPDIR/template_${group}/work/*_totemp_reslice_${sub}.nii.gz \
      -mean -o $INITTEMPDIR/template_${group}/work/template_${group}_${sub}.nii.gz
  else
    c3d $SHOOTDIR/template_${group}/*/iter_$((GSHOOT_NITER-1))/reslice*${sub}.nii.gz \
      -mean -o $SHOOTDIR/template_${group}/template/template_${group}_gshoot_${sub}.nii.gz
  fi
}

##########################################################
function RegisterInitTemp()
{
  # choose the best template to be the super template
  PREFIX=RST${expid}
  for ((iter=0;iter<$INITITER;iter++)); do
  #for ((iter=0;iter<0;iter++)); do
    for group in ${groups[*]}; do

      # prepare for the first iteration
      IDS=($(cat $ALLPATHDIR/IDSide.txt))
      INITWORKDIR=$INITTEMPDIR/template_${group}/work/
      mkdir -p $INITWORKDIR
      if [[ $iter -eq 0 ]]; then
        for sub in ${LABEL_IDS[*]} seg; do
          cp $MSTTEMPDIR/template_${group}/template_${group}_${sub}.nii.gz \
             $INITWORKDIR/template_${group}_${sub}.nii.gz
        done
      fi

      # Back up template
      ITDIR=$INITWORKDIR/$(printf iter_%02d $iter)
      mkdir -p $ITDIR
      cp -a $INITWORKDIR/template_${group}_*.nii.gz $ITDIR/

      # run registration
      for id in ${IDS[*]}; do

        # submit job
        qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${group}_${id}" $0 \
           RegisterInitTemp_sub $group $id
        sleep 0.1

      done
    done

    # Wait for completion
    qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

    # If this is the last iteration, we don't want to recompute the template
    if [[ $iter -lt $((INITITER-1)) ]]; then
    for group in ${groups[*]}; do
      for sub in ${KINDS[*]}; do

        qsub -cwd -o $DUMPDIR -j y \
             -q all.q,basic.q \
             -l h_vmem=4.1G,s_vmem=4G \
             -N "${PREFIX}_${group}_${sub}" $0 \
             GenerateTemp_sub $group $sub INITREG
        sleep 0.1

      done
    done

    # Wait for completion
    qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

    # Create the segmentation for the template
    for group in ${groups[*]}; do
      INITWORKDIR=$INITTEMPDIR/template_${group}/work/
      c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITWORKDIR/template_${group}_${sub}.nii.gz; done) \
        -vote -type ushort \
        -o $INITWORKDIR/template_${group}_seg.nii.gz
    done

    fi

  done

  # generate meshes and subsample
  for group in ${groups[*]}; do
  INITWORKDIR=$INITTEMPDIR/template_${group}/work/
  ITERDIR=$INITWORKDIR/$(printf iter_%02d $((INITITER-1)))

  # modify the background prob map
  cp $ITERDIR/template_${group}_BKG.nii.gz \
     $ITERDIR/template_${group}_BKG_backup.nii.gz
  c3d $ITERDIR/template_${group}_seg.nii.gz \
    -binarize -dilate 1 10x10x10vox \
    -scale -1 -shift 1 \
    $ITERDIR/template_${group}_BKG.nii.gz \
    -add -clip 0 1 \
    -o $ITERDIR/template_${group}_BKG.nii.gz

  for ((i=0;i<${#MESHES[*]};i++)); do

    LABELDEF=(${MESHESDEF[i]})
    echo $LABELDEF

    # generate meshes
    c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$ITERDIR/template_${group}_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
      -accum -add -endaccum \
      -o $ITERDIR/template_${group}_${MESHES[i]}_vtk.nii.gz

    vtklevelset $ITERDIR/template_${group}_${MESHES[i]}_vtk.nii.gz \
      $ITERDIR/template_${group}_${MESHES[i]}.stl 0.0

    vtklevelset $ITERDIR/template_${group}_${MESHES[i]}_vtk.nii.gz \
      $ITERDIR/template_${group}_${MESHES[i]}.vtk 0.0

    # perform subsampling
    SAM=$ITERDIR/template_${group}_${MESHES[i]}_sampled.ply
    /home/longxie/pkg/mesh_poisson_sample \
      $ITERDIR/template_${group}_${MESHES[i]}.stl \
      $SAM ${SAMPLEPOINT[i]}

  done

  # combine and convert into vtk mesh
  TEMPLATE=$ITERDIR/template_${group}_MRGcombined_sampled.vtk
  NV=0
  for ((i=0;i<${#COMBINEMESHES[*]};i++)); do
    SAM=$ITERDIR/template_${group}_${COMBINEMESHES[i]}_sampled.ply
    NV=$(($NV+$(cat $SAM | grep 'element vertex' | awk '{print $3}')))
  done

  # Write the header of the VTK file
  echo "# vtk DataFile Version 4.0" > $TEMPLATE
  echo "vtk output" >> $TEMPLATE
  echo "ASCII" >> $TEMPLATE
  echo "DATASET POLYDATA" >> $TEMPLATE
  echo "POINTS ${NV} float" >> $TEMPLATE

  # Write the point coordinates
  for ((i=0;i<${#COMBINEMESHES[*]};i++)); do
    SAM=$ITERDIR/template_${group}_${COMBINEMESHES[i]}_sampled.ply
    NV=$(cat $SAM | grep 'element vertex' | awk '{print $3}')
    grep -A $NV end_header $SAM | tail -n $NV >> $TEMPLATE
  done

  done
}

function RegisterInitTemp_sub()
{
  group=$1
  id=$2
  INITWORKDIR=$INITTEMPDIR/template_${group}/work/

  CMD=""
  for sub in ${LABEL_IDS[*]}; do
    CMD="$CMD -w 1 -i $INITWORKDIR/template_${group}_${sub}.nii.gz $MSTREGDIR/template_${group}/${id}/final/${id}_to_MSTRoot_reslice_${sub}.nii.gz "
  done

  # Run greedy between these two images
  WARP=$INITWORKDIR/${id}_totempWarp.nii.gz

  if [[ ! -f $WARP ]]; then

  # Run greedy between these two images
  c3d $INITWORKDIR/template_${group}_seg.nii.gz -dup \
      $MSTREGDIR/template_${group}/${id}/final/${id}_to_MSTRoot_reslice_seg.nii.gz \
      -int 0 -reslice-identity \
      -add -binarize -dilate 1 10x10x10vox \
      -o $TMPDIR/${id}_totemp_mask.nii.gz

  greedy -d 3 -threads $NSLOTS \
    $CMD \
    -n 120x120x40 \
    -gm $TMPDIR/${id}_totemp_mask.nii.gz \
    -s 1.5vox 1vox \
    -e 0.5 \
    -o $WARP

  fi

  # Reslice the segmentations from raw space
  RM=""
  for sub in ${KINDS[*]}; do
    RM="$RM -rm $DATADIR/${id}_${sub}.nii.gz $INITWORKDIR/${id}_totemp_reslice_${sub}.nii.gz"
  done

  WARPCHAIN=$(cat $MSTREGDIR/template_${group}/$id/final/chain_unwarp_to_final.txt)

  greedy -d 3 \
    -rf $INITWORKDIR/template_${group}_BKG.nii.gz \
    $RM \
    -r $WARP $WARPCHAIN

  c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITWORKDIR/${id}_totemp_reslice_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $INITWORKDIR/${id}_totemp_reslice_seg.nii.gz
}

##########################################################
function QC_shooting()
{
  #for group in ${groups[*]}; do
  for group in 2; do

      ITERDIR=$SHOOTDIR/template_${group}/shape_avg/iter_$((GSHOOT_NITER-1))
      VOLITERDIR=$SHOOTDIR/template_${group}/template/iter_$((GSHOOT_NITER-1))

      itksnap \
          -g $VOLITERDIR/template_${group}_gshoot_seg.nii.gz \
          -s $VOLITERDIR/template_${group}_gshoot_seg.nii.gz &

      IDS=($(cat $ALLPATHDIR/IDSide_${group}.txt))
      for ((idx=1;idx<=${#IDS[*]};idx++)); do

        id=${IDS[$((idx-1))]}
        SHOOTSUBJDIR=$SHOOTDIR/template_${group}/$id/iter_$((GSHOOT_NITER-1))

        itksnap \
          -g $SHOOTSUBJDIR/reslice_${id}_shooting_to_template_seg.nii.gz \
          -s $SHOOTSUBJDIR/reslice_${id}_shooting_to_template_seg.nii.gz

    done
  done
}

function ShootingCorrection()
{
  START_ITER=0
  END_ITER=$GSHOOT_NITER
  #END_ITER=1
  mkdir -p $SHOOTDIR

  # Loop
  for ((iter=$START_ITER;iter<$END_ITER;iter++)); do

    # loop through all subjects
    # optional: speed up shooting using GPU
    if [[ 1 == 1 ]]; then
    for group in ${groups[*]}; do

      IDS=($(cat $ALLPATHDIR/IDSide.txt))
      mkdir -p $SHOOTDIR/template_${group}/

      if [[ $iter -eq 0 ]]; then
        #LANDMARKS=$MSTTEMPDIR/template_${group}/template_${group}_MRGcombined_sampled.vtk
        INITWORKDIR=$INITTEMPDIR/template_${group}/work/
        ITERDIR=$INITWORKDIR/$(printf iter_%02d $((INITITER-1)))
        LANDMARKS=$ITERDIR/template_${group}_MRGcombined_sampled.vtk
        # generate ref space for the template
        if [ ! -f  $SHOOTDIR/template_${group}/refspace_${group}.nii.gz ]; then
        SEGS=($(ls $DATADIR/*seg_orig.nii.gz))
        c3d ${SEGS[0]} -resample 50x50x100% -pad 50x50x25vox 50x50x25vox 0 -dup  -popas REF \
          $(for seg in ${SEGS[*]}; do echo " -push REF $seg -int 0 -reslice-identity -add "; done) \
          -thresh 0.1 inf 1 0 \
          -trim 10vox \
          -resample-mm 0.4x0.4x0.4mm \
          -o $SHOOTDIR/template_${group}/refspace_${group}.nii.gz

        #/data/picsl/longxie/pkg/bin/antsbin/bin/AverageImages 3 $SHOOTDIR/template_${group}/refspace_${group}.nii.gz 0 \
        #  $(for id in ${IDS[*]}; do echo $DATADIR/${id}_seg_orig.nii.gz; done)
        #c3d $SHOOTDIR/template_${group}/refspace_${group}.nii.gz \
        #  -thresh 0.0001 inf 1 0 \
        #  -trim 20vox \
        #  -resample-mm 0.4x0.4x0.4mm \
        #  -o $SHOOTDIR/template_${group}/refspace_${group}.nii.gz
        fi
      else
        LANDMARKS=$SHOOTDIR/template_${group}/shape_avg/iter_$((iter-1))/shavg_landmarks.vtk
      fi

      for ((idx=1;idx<=${#IDS[*]};idx++)); do

        id=${IDS[$((idx-1))]}

        # optional: speed up shooting using GPU
        PREFIX=SCSP${expid}
        qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=6.1G,s_vmem=6G \
           -N "${PREFIX}_${group}_${idx}_${iter}" $0 \
           ShootingCorrectionSubjPre_sub $group $idx $id $iter $LANDMARKS

      done
    done
    fi

    # Wait for completion
     qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

    # perform shooting
    for group in ${groups[*]}; do

      IDS=($(cat $ALLPATHDIR/IDSide.txt))
      mkdir -p $SHOOTDIR/template_${group}/

      if [[ $iter -eq 0 ]]; then
        #LANDMARKS=$MSTTEMPDIR/template_${group}/template_${group}_MRGcombined_sampled.vtk
        INITWORKDIR=$INITTEMPDIR/template_${group}/work/
        ITERDIR=$INITWORKDIR/$(printf iter_%02d $((INITITER-1)))
        LANDMARKS=$ITERDIR/template_${group}_MRGcombined_sampled.vtk
        # generate ref space for the template
        if [ ! -f  $SHOOTDIR/template_${group}/refspace_${group}.nii.gz ]; then
        /data/picsl/longxie/pkg/bin/antsbin/bin/AverageImages 3 $SHOOTDIR/template_${group}/refspace_${group}.nii.gz 0 \
          $(for id in ${IDS[*]}; do echo $DATADIR/${id}_seg_orig.nii.gz; done)
        c3d $SHOOTDIR/template_${group}/refspace_${group}.nii.gz \
          -thresh 0.0001 inf 1 0 \
          -trim 20vox \
          -resample-mm 0.4x0.4x0.4mm \
          -o $SHOOTDIR/template_${group}/refspace_${group}.nii.gz
        fi
      else
        LANDMARKS=$SHOOTDIR/template_${group}/shape_avg/iter_$((iter-1))/shavg_landmarks.vtk
      fi

      for ((idx=1;idx<=${#IDS[*]};idx++)); do

        id=${IDS[$((idx-1))]}

        #PREFIX=GPUSCS${expid}
        #qsub -cwd -o $DUMPDIR -j y \
        #   -q gpu.q \
        #   -l h_vmem=6.1G,s_vmem=6G \
        #   -N "${PREFIX}_${group}_${idx}_${iter}" $0 \
        #   GPUShootingCorrectionSubj_sub $group $idx $id $iter $LANDMARKS

        #set +e
        #status=$(qstat | grep GPUS)
        #while [[ $status != "" ]]; do
        #  sleep 10
        #  status=$(qstat | grep GPUS)
        #  echo "wait till previous job is finished"
        #done
        #set -e


        # submit jobs
        PREFIX=SCS${expid}
        qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=6.1G,s_vmem=6G \
           -N "${PREFIX}_${group}_${idx}_${iter}" $0 \
           ShootingCorrectionSubj_sub $group $idx $id $iter $LANDMARKS
        sleep 0.1

      done
    done

    # Wait for completion
    qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

    # perform averaging
    PREFIX=SCA${expid}
    for group in ${groups[*]}; do

    if [[ $iter -eq 0 ]]; then
      #LANDMARKS=$MSTTEMPDIR/template_${group}/template_${group}_MRGcombined_sampled.vtk
      INITWORKDIR=$INITTEMPDIR/template_${group}/work/
      ITERDIR=$INITWORKDIR/$(printf iter_%02d $((INITITER-1)))
      LANDMARKS=$ITERDIR/template_${group}_MRGcombined_sampled.vtk
      LANDMARKSDIR=$ITERDIR/
    else
      LANDMARKS=$SHOOTDIR/template_${group}/shape_avg/iter_$((iter-1))/shavg_landmarks.vtk
      LANDMARKSDIR=$SHOOTDIR/template_${group}/shape_avg/iter_$((iter-1))/
    fi

    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=8.1G,s_vmem=8G \
         -N "${PREFIX}_shape_${group}_${iter}" $0 \
         ShootingCorrectionShapeAvg_sub $group $iter $LANDMARKS $LANDMARKSDIR

    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=15.1G,s_vmem=15G \
         -N "${PREFIX}_volume_${group}_${iter}" $0 \
         ShootingCorrectionVolumeAvg_sub $group $iter

    done

    # Wait for completion
    qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  done
}

function ShootingCorrectionSubjPre_sub()
{
  # The index of the subject being registered using the MST
  group=${1?}
  idx=${2?}
  ID=${3?}
  Line=($(cat $ALLPATHDIR/IDSide.txt | head -n $idx | tail -n 1))
  side=$(echo $Line | rev | cut -d _ -f 1 | rev)
  idraw=$(echo $Line | sed -e "s/_${side}//g")
  fnraw=$(ls $ASHSRUNDIR | grep $idraw | grep orig)

  # add on transforms
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS_1=$IDENTITYFN
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
  fi

  # The iteration - 0 is the initial iteration, involves extra work
  iter=${4?}

  # The path to the landmarks
  LANDMARKS=${5?}

  # directory
  WORK=$SHOOTDIR/template_${group}/$ID
  WITER=$WORK/iter_${iter}
  mkdir -p $WITER

  # Reference space (root node in cm-rep space)
  REFSPACE=$SHOOTDIR/template_${group}/refspace_${group}.nii.gz

  # Result meshes
  TARGET=$WORK/shooting_target_native.vtk
  LM_PROCRUSTES=$WITER/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$WITER/target_to_root_procrustes.mat
  SHOOTING_WARP=$WITER/shooting_warp.nii.gz

  # Mesh containing the momenta
  MOMENTA=$WITER/shooting_momenta.vtk

  # Target-related stuff in the WORK directory that is only done in the first iter
  if [[ $iter -eq 0 ]]; then

    # Get the warp chain from file
    WARPCHAIN=$(cat $MSTREGDIR/template_${group}/${ID}/final/chain_unwarp_to_final.txt)

    # Apply the warp chain to the landmark mesh in template space, creating
    # the target locations for the geodesic shooting registration
    # -- this code works when the WARPCHAIN is empty (MST root)
    if [ ! -f $TARGET ]; then
    SUBJASHSRUNDIR=$ASHSRUNDIR/$fnraw
    greedy -d 3 \
      -rf $REFSPACE \
      -rs $LANDMARKS $TARGET \
      -r $INITTEMPDIR/template_${group}/work/${ID}_totempWarp.nii.gz \
         $WARPCHAIN \
         $ADD_ON_TRANS \
         $SUBJASHSRUNDIR/ants_t1_to_temp/greedy_t1_to_template_warp.nii.gz \
         $SUBJASHSRUNDIR/affine_t1_to_template/t1_to_template_affine.mat \
         $SUBJASHSRUNDIR/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat
    fi

  fi

  # Landmarks in reference space
  ln -sf $LANDMARKS $WITER/landmarks.vtk

  # Bring the target mesh back near the root mesh using procrustes alignment
  if [ ! -f $LM_PROCRUSTES_MAT ]; then
    vtkprocrustes $TARGET $LANDMARKS $LM_PROCRUSTES_MAT
  fi

  # Apply procrustes to the landmarks.
  #warpmesh $TARGET $LM_PROCRUSTES $LM_PROCRUSTES_MAT
  if [ ! -f $LM_PROCRUSTES ]; then
  greedy -d 3 \
    -rf $REFSPACE \
    -rs $TARGET $LM_PROCRUSTES \
    -r $LM_PROCRUSTES_MAT
  fi

}

function ShootingCorrectionSubj_sub()
{
  # The index of the subject being registered using the MST
  group=${1?}
  idx=${2?}
  ID=${3?}
  Line=($(cat $ALLPATHDIR/IDSide.txt | head -n $idx | tail -n 1))
  side=$(echo $Line | rev | cut -d _ -f 1 | rev)
  idraw=$(echo $Line | sed -e "s/_${side}//g")
  fnraw=$(ls $ASHSRUNDIR | grep $idraw | grep orig)

  # add on transforms
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
  fi

  # The iteration - 0 is the initial iteration, involves extra work
  iter=${4?}

  # The path to the landmarks
  LANDMARKS=${5?}

  # directory
  WORK=$SHOOTDIR/template_${group}/$ID
  WITER=$WORK/iter_${iter}
  mkdir -p $WITER

  # Reference space (root node in cm-rep space)
  REFSPACE=$SHOOTDIR/template_${group}/refspace_${group}.nii.gz

  # Result meshes
  TARGET=$WORK/shooting_target_native.vtk
  LM_PROCRUSTES=$WITER/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$WITER/target_to_root_procrustes.mat
  SHOOTING_WARP=$WITER/shooting_warp.nii.gz

  # Mesh containing the momenta
  MOMENTA=$WITER/shooting_momenta.vtk

  # Target-related stuff in the WORK directory that is only done in the first iter
  if [[ $iter -eq 0 ]]; then

    # Get the warp chain from file
    WARPCHAIN=$(cat $MSTREGDIR/template_${group}/${ID}/final/chain_unwarp_to_final.txt)

    # Apply the warp chain to the landmark mesh in template space, creating
    # the target locations for the geodesic shooting registration
    # -- this code works when the WARPCHAIN is empty (MST root)
    if [ ! -f $TARGET ]; then
    SUBJASHSRUNDIR=$ASHSRUNDIR/$fnraw
    greedy -d 3 \
      -rf $REFSPACE \
      -rs $LANDMARKS $TARGET \
      -r $INITTEMPDIR/template_${group}/work/${ID}_totempWarp.nii.gz \
         $WARPCHAIN \
         $ADD_ON_TRANS \
         $SUBJASHSRUNDIR/ants_t1_to_temp/greedy_t1_to_template_warp.nii.gz \
         $SUBJASHSRUNDIR/affine_t1_to_template/t1_to_template_affine.mat \
         $SUBJASHSRUNDIR/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat
    fi

  fi

  # Landmarks in reference space
  ln -sf $LANDMARKS $WITER/landmarks.vtk

  # Bring the target mesh back near the root mesh using procrustes alignment
  if [ ! -f $LM_PROCRUSTES_MAT ]; then
    vtkprocrustes $TARGET $LANDMARKS $LM_PROCRUSTES_MAT
  fi

  # Apply procrustes to the landmarks.
  #warpmesh $TARGET $LM_PROCRUSTES $LM_PROCRUSTES_MAT
  if [ ! -f $LM_PROCRUSTES ]; then
  greedy -d 3 \
    -rf $REFSPACE \
    -rs $TARGET $LM_PROCRUSTES \
    -r $LM_PROCRUSTES_MAT
  fi

  # Perform geodesic shooting between the procrustes landmarks and the
  # warped landmarks - this is going to allow us to interpolate the correspondence
  # found by the MST to the rest of the images
  if [[ ! -f $MOMENTA ]]; then

    time lmshoot -d 3 \
      -m $LANDMARKS $LM_PROCRUSTES \
      -s 2.0 -l 5000 -n 40 -i 300 0 -f \
      -o $MOMENTA \

  fi

  if [[ ! -f $SHOOTING_WARP ]]; then

    # Convert the shooting result into a warp
    lmtowarp -d 3 -n 40 -r $REFSPACE \
      -m $MOMENTA -o $SHOOTING_WARP \
      -s 2.0

  fi

  # Warp the native space image into the template
  if [ ! -f $WITER/reslice_${ID}_shooting_to_template_seg.nii.gz ]; then

  RM=""
  for sub in ${LABEL_IDS[*]}; do
    RM="$RM -rm $DATADIR/${ID}_${sub}_orig.nii.gz $WITER/reslice_${ID}_shooting_to_template_${sub}.nii.gz"
  done

  greedy -d 3 \
    -rf $REFSPACE \
    $RM \
    -r $SHOOTING_WARP \
       $LM_PROCRUSTES_MAT,-1

  # Create seg
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $WITER/reslice_${ID}_shooting_to_template_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $WITER/reslice_${ID}_shooting_to_template_seg.nii.gz

  fi
}

function GPUShootingCorrectionSubj_sub()
{
  # The index of the subject being registered using the MST
  group=${1?}
  idx=${2?}
  ID=${3?}
  Line=($(cat $ALLPATHDIR/IDSide.txt | head -n $idx | tail -n 1))
  side=$(echo $Line | rev | cut -d _ -f 1 | rev)
  idraw=$(echo $Line | sed -e "s/_${side}//g")
  fnraw=$(ls $ASHSRUNDIR | grep $idraw | grep orig)

  # add on transforms
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
  fi

  # The iteration - 0 is the initial iteration, involves extra work
  iter=${4?}

  # The path to the landmarks
  LANDMARKS=${5?}

  # directory
  WORK=$SHOOTDIR/template_${group}/$ID
  WITER=$WORK/iter_${iter}
  mkdir -p $WITER

  # Reference space (root node in cm-rep space)
  REFSPACE=$SHOOTDIR/template_${group}/refspace_${group}.nii.gz

  # Result meshes
  TARGET=$WORK/shooting_target_native.vtk
  LM_PROCRUSTES=$WITER/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$WITER/target_to_root_procrustes.mat
  SHOOTING_WARP=$WITER/shooting_warp.nii.gz

  # Mesh containing the momenta
  MOMENTA=$WITER/shooting_momenta.vtk

  # Target-related stuff in the WORK directory that is only done in the first iter
  if [[ $iter -eq 0 ]]; then

    # Get the warp chain from file
    WARPCHAIN=$(cat $MSTREGDIR/template_${group}/${ID}/final/chain_unwarp_to_final.txt)

    # Apply the warp chain to the landmark mesh in template space, creating
    # the target locations for the geodesic shooting registration
    # -- this code works when the WARPCHAIN is empty (MST root)
    if [ ! -f $TARGET ]; then
    SUBJASHSRUNDIR=$ASHSRUNDIR/$fnraw
    greedy -d 3 \
      -rf $REFSPACE \
      -rs $LANDMARKS $TARGET \
      -r $INITTEMPDIR/template_${group}/work/${ID}_totempWarp.nii.gz \
         $WARPCHAIN \
         $ADD_ON_TRANS \
         $SUBJASHSRUNDIR/ants_t1_to_temp/greedy_t1_to_template_warp.nii.gz \
         $SUBJASHSRUNDIR/affine_t1_to_template/t1_to_template_affine.mat \
         $SUBJASHSRUNDIR/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat
     fi

  fi

  # Landmarks in reference space
  ln -sf $LANDMARKS $WITER/landmarks.vtk

  # Bring the target mesh back near the root mesh using procrustes alignment
  if [ ! -f $LM_PROCRUSTES_MAT ]; then
    vtkprocrustes $TARGET $LANDMARKS $LM_PROCRUSTES_MAT
  fi

  # Apply procrustes to the landmarks.
  #warpmesh $TARGET $LM_PROCRUSTES $LM_PROCRUSTES_MAT
  if [ ! -f $LM_PROCRUSTES ]; then
  #RS=""
  #for sub in ${MESHES[*]}; do
  #  RS="$RS -rs $WORK/template_${group}_to_${ID}_${sub}_native.vtk $WITER/template_${group}_to_${ID}_${sub}_procrustes.vtk"
  #done

  greedy -d 3 \
    -rf $REFSPACE \
    -rs $TARGET $LM_PROCRUSTES \
    -r $LM_PROCRUSTES_MAT
  fi

  # Perform geodesic shooting between the procrustes landmarks and the
  # warped landmarks - this is going to allow us to interpolate the correspondence
  # found by the MST to the rest of the images
  if [[ ! -f $MOMENTA ]]; then

    time /home/cwang/pcmrep/PointSetGeodesicShooting_CUDA_c00/lmshoot_cuda -d 3 \
      -m $LANDMARKS $LM_PROCRUSTES \
      -s 2.0 -l 5000 -n 40 -i 300 0 \
      -o $MOMENTA \

  fi
}

function ShootingCorrectionShapeAvg_sub()
{
  # The iteration - 0 is the initial iteration, involves extra worki
  group=${1?}
  iter=${2?}

  # The path to the landmarks
  SRC_LANDMARKS=${3?}
  SRC_LANDMARKSDIR=${4?}

  # directory
  WORK=$SHOOTDIR/template_${group}/shape_avg/iter_${iter}
  mkdir -p $WORK

  # Reference space (root node in cm-rep space)
  REFSPACE=$MSTTEMPDIR/template_${group}/template_${group}_seg.nii.gz

  # The result landmarks - after shape averating
  SHAVG_LANDMARKS_NOPROC=$WORK/shavg_landmarks_noprocrustes.vtk
  SHAVG_LANDMARKS=$WORK/shavg_landmarks.vtk

  # Average the momentum maps from the previous iteration
  avgmesharr \
    $SHOOTDIR/template_${group}/*/iter_${iter}/shooting_momenta.vtk \
    InitialMomentum $SRC_LANDMARKS $WORK/average_momenta.vtk

  # Perform the shooting and generate warp
  if [ ! -f $WORK/average_momenta_warp.nii.gz ]; then
    lmtowarp \
      -d 3 -n 40 -r $REFSPACE \
      -m $WORK/average_momenta.vtk \
      -o $WORK/average_momenta_warp.nii.gz -s 2.0
  fi

  # Apply the warp to the landmarks to bring them to shape-averaging position
  greedy \
    -d 3 -rs $SRC_LANDMARKS $SHAVG_LANDMARKS_NOPROC \
    -rf $REFSPACE -r $WORK/average_momenta_warp.nii.gz

  # This transformation of the landmarks can cause shrinkage of the template. This
  # is not at all what we want in the template, we actually want the template to
  # keep its size during this iterative process. The way to correct this is to perform
  # procrustes between the source lanmarks and the new shape average
  vtkprocrustes $SRC_LANDMARKS $SHAVG_LANDMARKS_NOPROC $WORK/residual_procrustes.mat \
    | grep RMS_ | tee $WORK/procrustes_metric.txt

  # Applying the inverse of this procrustes to the SHAVG_LANDMARKS_NOPROC gives the new
  # template landmarks which are shape averaged but still the same size as the original
  # template.
  RS=""
  for sub in ${MESHES[*]}; do
    RS="$RS -rs $SRC_LANDMARKSDIR/template_${group}_${sub}.vtk $WORK/template_${group}_${sub}.vtk "
  done

  greedy \
    -d 3 -rs $SRC_LANDMARKS $SHAVG_LANDMARKS \
    $RS \
    -rf $REFSPACE \
    -r $WORK/average_momenta_warp.nii.gz $WORK/residual_procrustes.mat,-1
}

function ShootingCorrectionVolumeAvg_sub()
{
  # The iteration - 0 is the initial iteration, involves extra work
  group=${1?}
  iter=${2?}
  ITERDIR=$SHOOTDIR/template_${group}/template/iter_${iter}
  mkdir -p $ITERDIR

  # create mean image for different labels
  for sub in ${LABEL_IDS[*]}; do

    IMGS=($(ls $SHOOTDIR/template_${group}/*/iter_${iter}/reslice*${sub}.nii.gz))
    #c3d $(for file in $(ls $SHOOTDIR/template_${group}/*/iter_${iter}/reslice*${sub}.nii.gz); do 
    c3d ${IMGS[0]} $(for ((i=1;i<${#IMGS[*]};i++)); do echo " ${IMGS[i]} -add "; done) -scale $(echo 1/${#IMGS[*]} | bc -l) -o $ITERDIR/template_${group}_gshoot_${sub}.nii.gz

    #c3d $SHOOTDIR/template_${group}/*/iter_${iter}/reslice*${sub}.nii.gz \
    #  -mean -o $ITERDIR/template_${group}_gshoot_${sub}.nii.gz

  done

  # generate seg
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $ITERDIR/template_${group}_gshoot_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $ITERDIR/template_${group}_gshoot_seg.nii.gz

  # generate meshes
  for ((i=0;i<${#MESHES[*]};i++)); do

    LABELDEF=(${MESHESDEF[i]})
    echo $LABELDEF

    c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$ITERDIR/template_${group}_gshoot_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
      -accum -add -endaccum \
      -o $ITERDIR/template_${group}_gshoot_${MESHES[i]}_vtk.nii.gz

    vtklevelset \
      $ITERDIR/template_${group}_gshoot_${MESHES[i]}_vtk.nii.gz \
      $ITERDIR/template_${group}_gshoot_${MESHES[i]}.vtk 0.0

  done
}

##########################################################
function ShootCorrectedMesh()
{
  for group in ${groups[*]}; do

    IDS=($(cat $ALLPATHDIR/IDSide.txt))

    for ((idx=1;idx<=${#IDS[*]};idx++)); do

      id=${IDS[$((idx-1))]}

      # optional: speed up shooting using GPU
      PREFIX=SCM${expid}
      qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX}_${group}_${idx}" $0 \
         ShootCorrectedMesh_sub $group $idx $id
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
     -q all.q,basic.q \
     -hold_jid "${PREFIX}_*" -sync y -b y \
     sleep 1
}

function ShootCorrectedMesh_sub()
{
  # The index of the subject being registered using the MST
  group=${1?}
  idx=${2?}
  ID=${3?}
  Line=($(cat $ALLPATHDIR/IDSide.txt | head -n $idx | tail -n 1))
  side=$(echo $Line | rev | cut -d _ -f 1 | rev)
  idraw=$(echo $Line | sed -e "s/_${side}//g")
  fnraw=$(ls $ASHSRUNDIR | grep $idraw | grep orig)

  # directory
  WORK=$SHOOTDIR/template_${group}/$ID
  WITER=$WORK/iter_$((GSHOOT_NITER-1))

  # Reference space (root node in cm-rep space)
  REFSPACE=$SHOOTDIR/template_${group}/refspace_${group}.nii.gz

  # Result meshes
  TARGET=$WORK/shooting_target_native.vtk
  LM_PROCRUSTES=$WITER/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$WITER/target_to_root_procrustes.mat
  SHOOTING_WARP=$WITER/shooting_warp.nii.gz

  # Mesh containing the momenta
  MOMENTA=$WITER/shooting_momenta.vtk

  # template mesh
  TEMPITERDIR=$SHOOTDIR/template_${group}/template/iter_$((GSHOOT_NITER-1))
  TEMPMESH=$TEMPITERDIR/template_${group}_gshoot_MRG.vtk

  # warp it back to subject space
  #greedy -d 3 \
  #  -rf $REFSPACE \
  #  -rs $TEMPMESH $WITER/${ID}_MRG_tempfit.vtk \
  #  -r $SHOOTING_WARP \
  #     $LM_PROCRUSTES_MAT,-1

  if [[ ! -f $WITER/${ID}_MRG_tempfit.vtk ]]; then
  $LMTOWARPDIR/lmtowarp -d 3 -n 40 -s 2.0 \
    -m $MOMENTA \
    -M $TEMPMESH $WITER/${ID}_MRG_tempfit.vtk

  greedy -d 3 \
    -rf $REFSPACE \
    -rs $WITER/${ID}_MRG_tempfit.vtk $WITER/${ID}_MRG_tempfit.vtk \
    -r $LM_PROCRUSTES_MAT,-1
  fi

  # get a mesh for the labels
  c3d $DATADIR/${ID}_seg_orig.nii.gz \
    -replace 2 0 8 0 9 0 \
    $TMPDIR/${ID}_seg_orig.nii.gz
  vtklevelset \
    $TMPDIR/${ID}_seg_orig.nii.gz \
    $WITER/${ID}_seg_orig.vtk 0.5

  # sample the label
  LABELIDS=(1 3 4 5 6 7)
  MESHLABELS=('CA' 'SUB' 'ERC' 'BA35' 'BA36' 'PHC')
  MESHES=""
  for ((i=0;i<${#MESHLABELS[*]};i++)); do
    c3d $TMPDIR/${ID}_seg_orig.nii.gz \
      -thresh ${LABELIDS[i]} ${LABELIDS[i]} 1 0 \
      -smooth 1.5vox \
      -o $TMPDIR/${ID}_seg_orig_${MESHLABELS[i]}.nii.gz
    mesh_image_sample \
      $WITER/${ID}_seg_orig.vtk \
      $TMPDIR/${ID}_seg_orig_${MESHLABELS[i]}.nii.gz \
      $TMPDIR/${ID}_seg_orig_${MESHLABELS[i]}.vtk \
      PROB
    MESHES="$MESHES $TMPDIR/${ID}_seg_orig_${MESHLABELS[i]}.vtk"
  done

  mesh_merge_arrays \
    -r $WITER/${ID}_seg_orig.vtk \
    $WITER/${ID}_seg_orig.vtk \
    PROB $MESHES

  # run matlab script to generate label
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/home/longxie/ASHS_PHC/thickness_newlabel/matlabcode');
    compute_mesh_label('$WITER/${ID}_seg_orig.vtk');
MATCODE

}

##########################################################
function FinalShooting()
{
  #################################
  # perform sampling
  group="1"
  IDS=($(cat $ALLPATHDIR/IDSide.txt))
  VOLTEMPDIR=$SHOOTDIR/template_${group}/template/iter_$((GSHOOT_NITER-1))/
  TEMPLATE=$VOLTEMPDIR/template_${group}_gshoot_MRGcombined_sampled.vtk
  iter="final"

  if [[ ! -f $TEMPLATE ]]; then

  # modify the background prob map
  cp $VOLTEMPDIR/template_${group}_gshoot_BKG.nii.gz \
     $VOLTEMPDIR/template_${group}_gshoot_BKG_backup.nii.gz
  c3d $VOLTEMPDIR/template_${group}_gshoot_seg.nii.gz \
    -binarize -dilate 1 10x10x10vox \
    -scale -1 -shift 1 \
    $VOLTEMPDIR/template_${group}_gshoot_BKG.nii.gz \
    -add -clip 0 1 \
    -o $VOLTEMPDIR/template_${group}_gshoot_BKG.nii.gz

  for ((i=0;i<${#MESHES[*]};i++)); do

    LABELDEF=(${MESHESDEF[i]})
    echo $LABELDEF

    c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$VOLTEMPDIR/template_${group}_gshoot_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
      -accum -add -endaccum \
      -o $VOLTEMPDIR/template_${group}_gshoot_${MESHES[i]}_vtk.nii.gz

    vtklevelset $VOLTEMPDIR/template_${group}_gshoot_${MESHES[i]}_vtk.nii.gz \
      $VOLTEMPDIR/template_${group}_gshoot_${MESHES[i]}.stl 0.0

    vtklevelset $VOLTEMPDIR/template_${group}_gshoot_${MESHES[i]}_vtk.nii.gz \
      $VOLTEMPDIR/template_${group}_gshoot_${MESHES[i]}.vtk 0.0

    # perform subsampling
    SAM=$VOLTEMPDIR/template_${group}_gshoot_${MESHES[i]}_sampled.ply
    /home/longxie/pkg/mesh_poisson_sample \
      $VOLTEMPDIR/template_${group}_gshoot_${MESHES[i]}.stl \
      $SAM ${SAMPLEPOINT[i]}

  done

  # combine and convert into vtk mesh
  rm -rf $TEMPLATE
  NV=0
  for ((i=0;i<${#COMBINEMESHES[*]};i++)); do
    SAM=$VOLTEMPDIR/template_${group}_gshoot_${COMBINEMESHES[i]}_sampled.ply
    NV=$(($NV+$(cat $SAM | grep 'element vertex' | awk '{print $3}')))
    grep -A $NV end_header $SAM | tail -n $NV >> $TEMPLATE
  done

  # Write the header of the VTK file
  echo "# vtk DataFile Version 4.0" > $TEMPLATE
  echo "vtk output" >> $TEMPLATE
  echo "ASCII" >> $TEMPLATE
  echo "DATASET POLYDATA" >> $TEMPLATE
  echo "POINTS ${NV} float" >> $TEMPLATE

  # Write the point coordinates
  for ((i=0;i<${#COMBINEMESHES[*]};i++)); do
    SAM=$VOLTEMPDIR/template_${group}_gshoot_${COMBINEMESHES[i]}_sampled.ply
    NV=$(cat $SAM | grep 'element vertex' | awk '{print $3}')
    grep -A $NV end_header $SAM | tail -n $NV >> $TEMPLATE
  done

  fi

  #################################
  # loop through all subjects
  # optional: speed up shooting using GPU
  if [[ 1 == 1 ]]; then
  REFSPACE=$SHOOTDIR/template_${group}/refspace_${group}.nii.gz
  LANDMARKS=$TEMPLATE

  for ((idx=1;idx<=${#IDS[*]};idx++)); do

    id=${IDS[$((idx-1))]}

    # optional: speed up shooting using GPU
    PREFIX=FSSP${expid}
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=6.1G,s_vmem=6G \
         -N "${PREFIX}_${group}_${idx}_${iter}" $0 \
         FinalShootingSubjPre_sub $group $idx $id final $LANDMARKS

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # perform shooting
  for ((idx=1;idx<=${#IDS[*]};idx++)); do

    id=${IDS[$((idx-1))]}
    #PREFIX=GPUFSS${expid}
    #qsub -cwd -o $DUMPDIR -j y \
    #   -q gpu.q \
    #   -l h_vmem=6.1G,s_vmem=6G \
    #   -N "${PREFIX}_${group}_${idx}_${iter}" $0 \
    #   GPUFinalShootingSubj_sub $group $idx $id final $LANDMARKS

    #set +e
    #status=$(qstat | grep GPUF)
    #while [[ $status != "" ]]; do
    #  sleep 10
    #  status=$(qstat | grep GPUF)
    #  echo "wait till previous job is finished"
    #done
    #set -e

    # submit jobs
    PREFIX=FSS${expid}
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=6.1G,s_vmem=6G \
         -N "${PREFIX}_${group}_${idx}_${iter}" $0 \
         FinalShootingSubj_sub $group $idx $id final $LANDMARKS
    sleep 0.1

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
  fi
}

function FinalShootingSubjPre_sub()
{
  # The index of the subject being registered using the MST
  group=${1?}
  idx=${2?}
  ID=${3?}
  Line=($(cat $ALLPATHDIR/IDSide.txt | head -n $idx | tail -n 1))
  side=$(echo $Line | rev | cut -d _ -f 1 | rev)
  idraw=$(echo $Line | sed -e "s/_${side}//g")
  fnraw=$(ls $ASHSRUNDIR | grep $idraw | grep orig)
  iter=${4?}
  LANDMARKS=${5?}

  # directory
  WORK=$SHOOTDIR/template_${group}/$ID
  PITER=$WORK/iter_$((GSHOOT_NITER-1))
  WITER=$WORK/iter_${iter}
  mkdir -p $WITER

  # Reference space (root node in cm-rep space)
  REFSPACE=$SHOOTDIR/template_${group}/refspace_${group}.nii.gz

  # Result meshes
  TARGET=$WITER/shooting_target_native.vtk
  LM_PROCRUSTES=$WITER/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$WITER/target_to_root_procrustes.mat
  SHOOTING_WARP=$WITER/shooting_warp.nii.gz

  # Mesh containing the momenta
  MOMENTA=$WITER/shooting_momenta.vtk

  # warp landmarks to subjects space
  PMOMENTA=$PITER/shooting_momenta.vtk
  PLM_PROCRUSTES_MAT=$PITER/target_to_root_procrustes.mat
  if [[ ! -f $TARGET ]]; then
    $LMTOWARPDIR/lmtowarp -d 3 -n 40 -s 2.0 \
      -m $PMOMENTA \
      -M $LANDMARKS $TARGET

    greedy -d 3 \
      -rf $REFSPACE \
      -rs $TARGET $TARGET \
      -r $PLM_PROCRUSTES_MAT,-1
  fi

  # Landmarks in reference space
  ln -sf $LANDMARKS $WITER/landmarks.vtk

  # Bring the target mesh back near the root mesh using procrustes alignment
  #if [ ! -f $LM_PROCRUSTES_MAT ]; then
    vtkprocrustes $TARGET $LANDMARKS $LM_PROCRUSTES_MAT
  #fi

  # Apply procrustes to the landmarks.
  #warpmesh $TARGET $LM_PROCRUSTES $LM_PROCRUSTES_MAT
  #if [ ! -f $LM_PROCRUSTES ]; then
  greedy -d 3 \
    -rf $REFSPACE \
    -rs $TARGET $LM_PROCRUSTES \
    -r $LM_PROCRUSTES_MAT
  #fi
}

function GPUFinalShootingSubj_sub()
{
  # The index of the subject being registered using the MST
  group=${1?}
  idx=${2?}
  ID=${3?}
  Line=($(cat $ALLPATHDIR/IDSide.txt | head -n $idx | tail -n 1))
  side=$(echo $Line | rev | cut -d _ -f 1 | rev)
  idraw=$(echo $Line | sed -e "s/_${side}//g")
  fnraw=$(ls $ASHSRUNDIR | grep $idraw | grep orig)
  iter=${4?}
  LANDMARKS=${5?}

  # directory
  WORK=$SHOOTDIR/template_${group}/$ID
  WITER=$WORK/iter_${iter}
  mkdir -p $WITER

  # Reference space (root node in cm-rep space)
  REFSPACE=$SHOOTDIR/template_${group}/refspace_${group}.nii.gz

  # Result meshes
  TARGET=$WITER/shooting_target_native.vtk
  LM_PROCRUSTES=$WITER/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$WITER/target_to_root_procrustes.mat
  SHOOTING_WARP=$WITER/shooting_warp.nii.gz

  # Mesh containing the momenta
  MOMENTA=$WITER/shooting_momenta.vtk

  # Perform geodesic shooting between the procrustes landmarks and the
  # warped landmarks - this is going to allow us to interpolate the correspondence
  # found by the MST to the rest of the images
  if [[ ! -f $MOMENTA ]]; then

    time /home/cwang/pcmrep/PointSetGeodesicShooting_CUDA_c00/lmshoot_cuda -d 3 \
      -m $LANDMARKS $LM_PROCRUSTES \
      -s 2.0 -l 5000 -n 40 -i 300 0 \
      -o $MOMENTA \

  fi
}

function FinalShootingSubj_sub()
{
  # The index of the subject being registered using the MST
  group=${1?}
  idx=${2?}
  ID=${3?}
  Line=($(cat $ALLPATHDIR/IDSide.txt | head -n $idx | tail -n 1))
  side=$(echo $Line | rev | cut -d _ -f 1 | rev)
  idraw=$(echo $Line | sed -e "s/_${side}//g")
  fnraw=$(ls $ASHSRUNDIR | grep $idraw | grep orig)
  iter=${4?}
  LANDMARKS=${5?}

  # directory
  WORK=$SHOOTDIR/template_${group}/$ID
  WITER=$WORK/iter_${iter}
  mkdir -p $WITER

  # Reference space (root node in cm-rep space)
  REFSPACE=$SHOOTDIR/template_${group}/refspace_${group}.nii.gz

  # Result meshes
  TARGET=$WITER/shooting_target_native.vtk
  LM_PROCRUSTES=$WITER/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$WITER/target_to_root_procrustes.mat
  SHOOTING_WARP=$WITER/shooting_warp.nii.gz

  # Mesh containing the momenta
  MOMENTA=$WITER/shooting_momenta.vtk

  if [[ ! -f $MOMENTA ]]; then

    time lmshoot -d 3 \
      -m $LANDMARKS $LM_PROCRUSTES \
      -s 2.0 -l 5000 -n 40 -i 300 0 \
      -o $MOMENTA \

  fi

  if [[ ! -f $SHOOTING_WARP ]]; then

    # Convert the shooting result into a warp
    lmtowarp -d 3 -n 40 -r $REFSPACE \
      -m $MOMENTA -o $SHOOTING_WARP \
      -s 2.0

  fi

  # Warp the native space image into the template
  if [ ! -f $WITER/reslice_${ID}_shooting_to_template_seg.nii.gz ]; then

  RM=""
  for sub in ${LABEL_IDS[*]}; do
    RM="$RM -rm $DATADIR/${ID}_${sub}_orig.nii.gz $WITER/reslice_${ID}_shooting_to_template_${sub}.nii.gz"
  done

  greedy -d 3 \
    -rf $REFSPACE \
    $RM \
    -r $SHOOTING_WARP \
       $LM_PROCRUSTES_MAT,-1

  # Create seg
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $WITER/reslice_${ID}_shooting_to_template_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $WITER/reslice_${ID}_shooting_to_template_seg.nii.gz

  fi
}

##########################################################
function ShapeStats()
{
  PREFIX=SS${expid}
  # submit jobs
  for group in ${groups[*]}; do
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=6.1G,s_vmem=6G \
         -N "${PREFIX}_${group}" $0 \
         ShapeStats_sub $group
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
     -q all.q,basic.q \
     -hold_jid "${PREFIX}_*" -sync y -b y \
     sleep 1
}

function ShapeStats_sub()
{
  group=$1
  VOLITERDIR=$SHOOTDIR/template_${group}/template/iter_$((GSHOOT_NITER-1))
  PCADIR=$VOLITERDIR/pca
  rm -rf $PCADIR/*
  mkdir -p $PCADIR
  IDs=$(cat $ALLPATHDIR/IDSide.txt)

  # Generate the momentum data in a format that is readable by R script
  idx=1
  echo "ID" > $PCADIR/IDs.csv
  for id in $IDs; do
    local MOMENTS=$SHOOTDIR/template_${group}/$id/iter_final/shooting_momenta.vtk
    echo $(echo $idx) $(dumpmeshattr $MOMENTS InitialMomentum) >> $PCADIR/initial_momenta.txt
    echo $idx >> $PCADIR/IDs.csv
    idx=$((idx+1))
  done

  # Generate an array of mesh coordinates for R
  dumpmeshpoints $VOLITERDIR/template_${group}_gshoot_MRGcombined_sampled.vtk | tail -n +3 > $PCADIR/template_points.txt

  # sample on probability maps
  cp $VOLITERDIR/template_${group}_gshoot_MRGcombined_sampled.vtk \
     $VOLITERDIR/template_${group}_gshoot_MRGcombined_sampled_labelprobs.vtk
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    mesh_image_sample \
      $VOLITERDIR/template_${group}_gshoot_MRGcombined_sampled_labelprobs.vtk \
      $VOLITERDIR/template_${group}_gshoot_${LABEL_IDS[i]}.nii.gz \
      $VOLITERDIR/template_${group}_gshoot_MRGcombined_sampled_labelprobs.vtk \
      Prob_${LABEL_IDS[i]}
  done

  # Generate the label probability data in the template space in a format that is readable by R script
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    echo Prob_${LABEL_IDS[i]} $(dumpmeshattr $VOLITERDIR/template_${group}_gshoot_MRGcombined_sampled_labelprobs.vtk Prob_${LABEL_IDS[i]}) >> $PCADIR/label_probs.txt
  done

  # Run the R statistics notebook
  #Rscript run_pca_oneroot.R
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/home/longxie/ASHS_PHC/thickness_newlabel/matlabcode');
    runPCA('$PCADIR/initial_momenta.txt', ...
           '$PCADIR/label_probs.txt', ...
           '$PCADIR/template_points.txt', ...
           '$PCADIR', '$PCADIR/template_allinfo_${group}.vtk');
MATCODE
}

##########################################################
function MakeMovies()
{
  PREFIX=MM${expid}
  for group in ${groups[*]}; do
    for what in mode1 mode2 mode3; do

      # submit jobs
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=6.1G,s_vmem=6G \
           -N "${PREFIX}_${group}_${what}" $0 \
           MakeMovies_sub $group $what

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
     -q all.q,basic.q \
     -hold_jid "${PREFIX}_*" -sync y -b y \
     sleep 1
}

function MakeMovies_sub()
{
  group=$1
  what=$2
  group=$1
  VOLITERDIR=$SHOOTDIR/template_${group}/template/iter_$((GSHOOT_NITER-1))
  PCADIR=$VOLITERDIR/pca
  MOVIEDIR=$PCADIR/movie_${what}
  mkdir -p $MOVIEDIR
  TEMPDIR=$PCADIR/tmp
  mkdir -p $TEMPDIR

  # Template
  local TEMPLATE=$VOLITERDIR/template_${group}_gshoot_seg.nii.gz
  local TEMPMESH=$VOLITERDIR/template_${group}_gshoot_MRG.vtk

  # Flow the mode vector forward and backward
  $LMTOWARPDIR/lmtowarp \
    -m $PCADIR/${what}_vector_neg.vtk -s 2.0 -n 40 -d 3 -a 4 \
    -M $TEMPMESH $TEMPDIR/neg_flow_%d.vtk

  $LMTOWARPDIR/lmtowarp \
    -m $PCADIR/${what}_vector_pos.vtk -s 2.0 -n 40 -d 3 -a 4 \
    -M $TEMPMESH $TEMPDIR/pos_flow_%d.vtk

  # Combine the flows into a single movie
  merge_neg_pos_movie $MOVIEDIR/movie_${group}_${what}_MRG_%02d.vtk $TEMPMESH $TEMPDIR
}

function merge_neg_pos_movie()
{
  local OUT_PATTERN=${1?}
  local CENTER=${2?}
  TEMPDIR=$3

  for ((i=4;i<=40;i+=4)); do
    local kneg=$((40-i))
    local kpos=$((40+i))
    cp $TEMPDIR/neg_flow_${i}.vtk $(printf $OUT_PATTERN $kneg)
    cp $TEMPDIR/pos_flow_${i}.vtk $(printf $OUT_PATTERN $kpos)
  done

  cp $CENTER $(printf $OUT_PATTERN 40)
}

##########################################################
function WarpTempToSubj()
{
  #mkdir -p $LABELWARPDIR

  # submit job for each subject
  PREFIX=WTTS${expid}
  for group in ${groups[*]}; do
    IDS=($(cat $ALLPATHDIR/IDSide.txt))
    for ((idx=1;idx<=${#IDS[*]};idx++)); do

      id=${IDS[$((idx-1))]}
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${group}_${idx}_${id}" $0 \
           WarpTempToSubj_sub $group $idx $id
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function WarpTempToSubj_sub()
{
  group=$1
  idx=$2
  id=$3
  SHOOTSUBJDIR=$SHOOTDIR/template_${group}/$id/iter_final
  VOLITERDIR=$SHOOTDIR/template_${group}/template/iter_$((GSHOOT_NITER-1))

  # convert ground truth label
  Line=($(cat $ALLPATHDIR/IDSide.txt | head -n $idx | tail -n 1))
  side=$(echo $Line | rev | cut -d _ -f 1 | rev)
  idraw=$(echo $Line | sed -e "s/_${side}//g")
  fnraw=$(ls $ASHSRUNDIR | grep $idraw | grep orig)
  GTSEG=$DATADIR/${id}_seg_orig.nii.gz

  # Volume template:  warp generate subject labels using mesh2img
  if [ ! -f $SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_seg.nii.gz ]; then
  TRANS=""
  TRANS1=""
  if [[ 1 == 1 ]]; then
  for ((i=0;i<${#MESHES[*]};i++)); do
    TRANS="$TRANS -M $VOLITERDIR/template_${group}_gshoot_${MESHES[i]}.vtk $SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_${MESHES[i]}.vtk"
    TRANS1="$TRANS1 -rs $SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_${MESHES[i]}.vtk $SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_${MESHES[i]}.vtk"
  done
  $LMTOWARPDIR/lmtowarp -d 3 -n 40 -s 2.0 \
      -m $SHOOTSUBJDIR/shooting_momenta.vtk \
      $TRANS

  greedy -d 3 \
    -rf $ASHSRUNDIR/$fnraw/mprage.nii.gz \
    $TRANS1 \
    -r $SHOOTSUBJDIR/target_to_root_procrustes.mat,-1

  for ((i=0;i<${#MESHES[*]};i++)); do

    mesh2img -f -vtk $SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_${MESHES[i]}.vtk \
      -a 0.3 0.3 0.3 4 \
      $TMPDIR/template_to_${id}_${group}_GSShoot_${MESHES[i]}.nii.gz

    c3d $GTSEG \
      $TMPDIR/template_to_${id}_${group}_GSShoot_${MESHES[i]}.nii.gz \
      -int 0 -reslice-identity \
      -o $SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_${MESHES[i]}.nii.gz

  done
  fi

  c3d $(for ((i=1;i<${#LABEL_IDS[*]};i++)); do \
    if [[ ${LABEL_IDS[i]} == "aCS" ]]; then \
      echo $SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_CS.nii.gz; \
    elif [[ ${LABEL_IDS[i]} == "pCS" ]]; then \
      continue; \
    else \
      echo $SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_${LABEL_IDS[i]}.nii.gz; \
    fi; done) \
    -foreach -sdt -scale -1 -endfor \
    -vote -shift 1 \
    $SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_NOBKG.nii.gz \
    -multiply -type ushort \
    -o $SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_seg.nii.gz
  fi
}

##########################################################
function EvalGSTempFit()
{
  EVALTMPDIR=$EVALDIR/tmp
  #rm -rf $EVALTMPDIR
  mkdir -p $EVALTMPDIR

  # submit jobs for evaluation
  if [[ 1 == 1 ]]; then
  PREFIX=EVGSTF${expid}
  for group in ${groups[*]}; do
    IDS=($(cat $ALLPATHDIR/IDSide.txt))
    for ((idx=1;idx<=${#IDS[*]};idx++)); do

      id=${IDS[$((idx-1))]}
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${group}_${idx}_${id}" $0 \
           EvalGSTempFit_sub $group $idx $id $EVALTMPDIR
      sleep 0.1

    done
  done
  fi

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1


  # Combine scores
  header=""
  for sub in ${EVALLABELS[*]}; do
    header="${header},${sub}"
  done

  echo "ID,SIDE,GRP${header}" \
    > $EVALDIR/overlap_GSTempFit.csv
  echo "ID,SIDE,GRP${header}" \
    > $EVALDIR/overlap_GSTempFitTemp.csv
  echo "ID,SIDE,GRP${header}" \
    > $EVALDIR/overlap_MSTTempFit.csv
  echo "ID,SIDE,GRP${header}" \
    > $EVALDIR/overlap_InitTemp.csv
  #for ((iter=0;iter<$GSHOOT_NITER;iter++)); do
  #  echo "ID,SIDE,GRP${header}" \
  #    > $EVALDIR/overlap_GSTempFit_iter${iter}.csv
  #done

  IDs=($(cat $SUBJ_TXT))
  for side in left right; do
    for id in ${IDs[*]}; do
      
      if [[ -f $EVALTMPDIR/${id}_${side}_GSTempFit_overlap.csv ]]; then
      cat $EVALTMPDIR/${id}_${side}_GSTempFit_overlap.csv \
        >> $EVALDIR/overlap_GSTempFit.csv
      else
      echo "" >> $EVALDIR/overlap_GSTempFit.csv
      fi

      if [[ -f $EVALTMPDIR/${id}_${side}_GSTempFitTemp_overlap.csv ]]; then
      cat $EVALTMPDIR/${id}_${side}_GSTempFitTemp_overlap.csv \
        >> $EVALDIR/overlap_GSTempFitTemp.csv
      else
      echo "" >> $EVALDIR/overlap_GSTempFitTemp.csv
      fi

      if [[ -f $EVALTMPDIR/${id}_${side}_MSTTempFit_overlap.csv ]]; then
      cat $EVALTMPDIR/${id}_${side}_MSTTempFit_overlap.csv \
        >> $EVALDIR/overlap_MSTTempFit.csv
      else
      echo "" >> $EVALDIR/overlap_MSTTempFit.csv
      fi

      if [[ -f $EVALTMPDIR/${id}_${side}_InitTemp_overlap.csv ]]; then
      cat $EVALTMPDIR/${id}_${side}_InitTemp_overlap.csv \
        >> $EVALDIR/overlap_InitTemp.csv
      else
      echo "" >> $EVALDIR/overlap_InitTemp.csv 
      fi 

      #for ((iter=0;iter<$GSHOOT_NITER;iter++)); do
      #  cat $EVALTMPDIR/${id}_${side}_GSTempFit_iter${iter}_overlap.csv \
      #    >> $EVALDIR/overlap_GSTempFit_iter${iter}.csv
      #done

    done
  done
}

function EvalGSTempFit_sub()
{
  group=$1
  idx=$2
  id=$3
  EVALTMPDIR=$4
  SHOOTSUBJDIR=$SHOOTDIR/template_${group}/$id/iter_final
  VOLITERDIR=$SHOOTDIR/template_${group}/template/iter_$((GSHOOT_NITER-1))

  # evaluate
  Line=($(cat $ALLPATHDIR/IDSide.txt | head -n $idx | tail -n 1))
  side=$(echo $Line | rev | cut -d _ -f 1 | rev)
  idraw=$(echo $Line | sed -e "s/_${side}//g")
  fnraw=$(ls $ASHSRUNDIR | grep $idraw | grep orig)
  GTSEG=$DATADIR/${id}_seg_orig.nii.gz
  GSSHOOTSEG=$SHOOTSUBJDIR/template_to_${id}_${group}_GSShoot_seg.nii.gz
  do_pair $GTSEG $GSSHOOTSEG
  echo "${idraw},${side},${group}${FULLOVL}" > \
    $EVALTMPDIR/${idraw}_${side}_GSTempFit_overlap.csv

  # GSHOOT in template space
  do_pair $SHOOTDIR/template_${group}/$id/iter_final/reslice_${id}_shooting_to_template_seg.nii.gz \
    $VOLITERDIR/template_${group}_gshoot_seg.nii.gz
  echo "${idraw},${side},${group}${FULLOVL}" > \
    $EVALTMPDIR/${idraw}_${side}_GSTempFitTemp_overlap.csv

  # MST template
  MSTWORKDIR=$MSTREGDIR/template_${group}
  do_pair $MSTWORKDIR/$id/final/${id}_to_MSTRoot_reslice_seg.nii.gz \
    $MSTTEMPDIR/template_${group}/template_${group}_seg.nii.gz
  echo "${idraw},${side},${group}${FULLOVL}" > \
    $EVALTMPDIR/${idraw}_${side}_MSTTempFit_overlap.csv

  # evaluate the overlap step by step
  INITWORKDIR=$INITTEMPDIR/template_${group}/work/
  do_pair $INITWORKDIR/${id}_totemp_reslice_seg.nii.gz \
    $INITWORKDIR/template_${group}_seg.nii.gz
  echo "${idraw},${side},${group}${FULLOVL}" > \
    $EVALTMPDIR/${idraw}_${side}_InitTemp_overlap.csv

  for ((iter=0;iter<$GSHOOT_NITER;iter++)); do
    VOLITERDIR=$SHOOTDIR/template_${group}/template/iter_$((iter))
    do_pair $SHOOTDIR/template_${group}/$id/iter_${iter}/reslice_${id}_shooting_to_template_seg.nii.gz \
      $VOLITERDIR/template_${group}_gshoot_seg.nii.gz
    echo "${idraw},${side},${group}${FULLOVL}" > \
      $EVALTMPDIR/${idraw}_${side}_GSTempFit_iter${iter}_overlap.csv
  done
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
    c3d \
      $seg_a -dup $seg_b -int 0 -reslice-identity \
      -foreach -replace $REPRULE -thresh 99 99 1 0 -endfor \
      -overlap 1 | tee $TMPDIR/ovl.txt

    # Get the full-extent overlap
    OVL=$(cat $TMPDIR/ovl.txt | grep OVL | awk -F '[ ,]+' '{print $6}')

    #echo $id ${LABELS[i]} full $OVL $DIST >> $out_file
    FULLOVL="${FULLOVL},${OVL}"

  done
}

##########################################################
function CopyShootMeshes()
{
  for group in ${groups[*]}; do
    IDS=($(cat $ALLPATHDIR/IDSide_${group}.txt))
    GRPDIR=$GSTEMPDIR/dataset/template_${group}
    mkdir -p $GRPDIR
    ITERDIR=$SHOOTDIR/template_${group}/shape_avg/iter_$((GSHOOT_NITER-2))
    VOLITERDIR=$SHOOTDIR/template_${group}/template/iter_$((GSHOOT_NITER-1))

    cp $ITERDIR/shavg_landmarks.vtk $GRPDIR/
    cp $ITERDIR/template_${group}_MRG.vtk $GRPDIR/

    for ((idx=1;idx<=${#IDS[*]};idx++)); do

      id=${IDS[$((idx-1))]}
      SUBJDIR=$GRPDIR/$id
      mkdir -p $SUBJDIR
      SHOOTSUBJDIR=$SHOOTDIR/template_${group}/$id/iter_$((GSHOOT_NITER-1))
      SHOOTSUBJGPUDIR=$SHOOTDIR/template_${group}/$id/iter_$((GSHOOT_NITER-1))_GPU
      cp $SHOOTSUBJDIR/shooting_target_procrustes.vtk $SUBJDIR/
      cp $SHOOTSUBJDIR/target_to_root_procrustes.mat $SUBJDIR/
      cp $SHOOTSUBJDIR/shooting_momenta.vtk $SUBJDIR
      cp $SHOOTDIR/template_${group}/$id/shooting_target_native.vtk $SUBJDIR/
      cp $SHOOTSUBJGPUDIR/shooting_momenta.vtk $SUBJDIR/shooting_momenta_gpu_cluster.vtk
      cp $SHOOTSUBJGPUDIR/template_${group}_to_${id}_MRG_procrustes.vtk $SUBJDIR/
      cp $SHOOTDIR/template_${group}/$id/template_${group}_to_${id}_MRG_native.vtk $SUBJDIR/


    done
  done
}

##########################################################
function reset_dir()
{
  rm -rf $DUMPDIR/*
}

##########################################################
if [[ $# -lt 2 ]]; then

  main

else

  cmd=$1
  shift
  $cmd $@

fi





