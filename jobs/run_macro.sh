#!/bin/bash

#root.exe Fun4All_G4_sPHENIX.C"("$n_events","$i","$z_width","$eta_min","$eta_max","$pt_min","$pt_max","$phi_min","$phi_max","\"/"$data"/"$parent"/"$directory"/G4_sPHENIX_"$config".root\"")"

#root.exe Fun4All_G4_sPHENIX.C"("$1","$2","$3","$4","$5","$6","$7","$8","$9","\"${10}/G4_sPHENIX_${11}.root\"")"
root.exe /direct/star+u/rcorliss/sphenix/macros/macros/g4simulations/Fun4All_G4_sPHENIX.C"("$1","$2","$3","$4","$5","$6","$7","$8","$9","\"G4_sPHENIX_${11}.root\"")"

echo /${10}/G4_sPHENIX_${11}.root
mv G4_sPHENIX_${11}.root_g4svtx_eval.root /${10}/G4_sPHENIX_${11}.root_g4svtxeval.root
mv PHG4TrackKalmanFitter_eval.root /${10}/G4_sPHENIX_${11}.root_g4kalman_eval.root
exit
