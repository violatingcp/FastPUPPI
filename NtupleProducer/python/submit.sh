#!/bin/bash

#for x in `eos ls eos/cms/store/mc/TTI2023Upg14D/SinglePionPlusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/`; do
#    bsub -q 8nh run.sh /store/mc/TTI2023Upg14D/SinglePionPlusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/ $x $PWD
#done

#for x in `eos ls eos/cms/store/mc/TTI2023Upg14D/SinglePionMinusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/`; do
#    bsub -q 8nh run.sh /store/mc/TTI2023Upg14D/SinglePionMinusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/ $x $PWD
#done

#for  x in `eos ls  eos/cms/store/mc/TTI2023Upg14D/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/`; do
#     bsub -q 8nh run.sh /store/mc/TTI2023Upg14D/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/ $x $PWD
#done

#for  x in `eos ls  eos/cms/store/mc/TTI2023Upg14D/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/`; do
#     bsub -q 8nh run.sh /store/mc/TTI2023Upg14D/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/ $x $PWD
#done

dir=/store/mc/TTI2023Upg14D/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/
for  x in `eos ls  $dir`; do
     bsub -q 8nh run.sh $dir $x $PWD Top
done


dir=/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/
for  x in `eos ls  $dir`; do
     bsub -q 8nh run.sh $dir $x $PWD Zmm
done