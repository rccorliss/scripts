#!/bin/bash


data=gpfs/mnt/gpfs04/sphenix/user/rcorliss/data #ross's path
#data=sphenix/u/mitay/Documents/latest/macros/macros/g4simulations/data_gen # must be absolute path, no slashes at beginning or end
parent=mom_scan
n_events=5000


#inst=/gpfs/mnt/gpfs04/sphenix/user/rcorliss/data/mom_scan/e-_pt0.9GeV_phi-180-180d_z0cm_eta-1.2-1.2/G4_sPHENIX_00zp.root_g4kalman_eval.root
nsubmitted=0

while read inst
do
  echo "$inst"

filename=${inst##*/}
fileprefix=${filename%%.root*}
config=${fileprefix##*_}

dirblock=${inst%/*}
dirname=${dirblock##*/}

particle=$(echo $dirname | cut -f1 -d_)

ptblock=$(echo $dirname | cut -f2 -d_)
pt=$(echo $ptblock | cut -f2 -dt | cut -f1 -dG)

phiblock=$(echo $dirname | cut -f3 -d_)
phi_min=-180
phi_max=180

zblock=$(echo $dirname | cut -f4 -d_)
z=0

etablock=$(echo $dirname | cut -f5 -d_)
eta_min=-1.2
eta_max=1.2

if [ nsubmitted == 0 ]
   then	 
       echo input $inst
       echo parses as:
       echo particle=$particle
       echo pt=$pt
       echo config=$config
       nsubmitted=1
fi

    directory="$particle"_pt"$pt"GeV_phi"$phi_min"-"$phi_max"d_z"$z"cm_eta"$eta_min"-"$eta_max"


    # check to see if there is a directory for this data set; create one if there isn't
    if [ ! -d "$data"/"$parent" ]; then
      mkdir "$data"/"$parent"
    fi

    if [ ! -d "$data"/"$parent"/"$directory" ]; then
      mkdir "$data"/"$parent"/"$directory"
    fi

    path="$data"/"$parent"/"$directory"

    echo condor_submit simple_job.job n_events=$n_events i=$i z_width=$z eta_min=$eta_min eta_max=$eta_max pt_min=$pt pt_max=$pt phi_min=$phi_min phi_max=$phi_max particle=$particle path=$path config=$config

done < "${1:-/dev/stdin}"
