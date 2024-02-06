#cd ../
source ../source_at_start.sh
#cd mPMTmapping
#make clean
#make
bash make_mPMTmap.sh 2024
#bash make_solidAngleMap_LB.sh 1330
#bash make_solidAngleMap_LB.sh 1331
#
#bash make_solidAngleMap_LB.sh 1324
#bash make_solidAngleMap_LB.sh 1325
#for i in {1100..1149}
#do
#	a=$(ls ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID""$i""*)
#	echo $a
#        ./bin/make_referenceMap_PMT-basedBins -f $a -o Maps/AttenuationTest_PMT-basedBins/AttenuationTest_hadded_wcsim_mPMTmapping_401nm_FileID""$i""_all.txt
#
#done
#
#for i in {1160..1189}
#do
#	a=$(ls ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID""$i""*)
#	echo $a
#        ./bin/make_referenceMap_PMT-basedBins -f $a -o Maps/AttenuationTest_PMT-basedBins/AttenuationTest_hadded_wcsim_mPMTmapping_401nm_FileID""$i""_all.txt
#done
#
#bash makeAttenuationFitScans.sh
#
#./bin/make_referenceMap_PMT-basedBins -f ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID800_Absff1.0e+11_Rayff1.0e+11_R10.00_all.root -o Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID800_Absff1.0e+11_Rayff1.0e+11_R10.00_all.txt

#./bin/splitReference_intoMaps -f Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID800_Absff1.0e+11_Rayff1.0e+11_R10.00_all.txt
#
#./bin/make_referenceMap_PMT-basedBins -f ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID801_Absff1.0e+11_Rayff1.0e+11_R20.00_all.root -o Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID801_Absff1.0e+11_Rayff1.0e+11_R20.00_all.txt
#
#./bin/splitReference_intoMaps -f Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID801_Absff1.0e+11_Rayff1.0e+11_R20.00_all.txt
#
##./bin/make_referenceMap_PMT-basedBins -f ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID802_Absff1.0e+11_Rayff1.0e+11_R40.00_all.root -o Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID802_Absff1.0e+11_Rayff1.0e+11_R40.00_all.txt
#
#./bin/splitReference_intoMaps -f Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID802_Absff1.0e+11_Rayff1.0e+11_R40.00_all.txt
#
##./bin/make_referenceMap_PMT-basedBins -f ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID803_Absff1.0e+11_Rayff1.0e+11_R80.00_all.root -o Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID803_Absff1.0e+11_Rayff1.0e+11_R80.00_all.txt
#
#./bin/splitReference_intoMaps -f Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID803_Absff1.0e+11_Rayff1.0e+11_R80.00_all.txt
#
##./bin/make_referenceMap_PMT-basedBins -f ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID804_Absff1.0e+11_Rayff1.0e+11_R120.00_all.root -o Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID804_Absff1.0e+11_Rayff1.0e+11_R120.00_all.txt
#
#./bin/splitReference_intoMaps -f Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID804_Absff1.0e+11_Rayff1.0e+11_R120.00_all.txt
#
##./bin/make_referenceMap_PMT-basedBins -f ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID805_Absff1.0e+11_Rayff1.0e+11_R140.00_all.root -o Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID805_Absff1.0e+11_Rayff1.0e+11_R140.00_all.txt
#
#./bin/splitReference_intoMaps -f Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID805_Absff1.0e+11_Rayff1.0e+11_R140.00_all.txt
#
##./bin/make_referenceMap_PMT-basedBins -f ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID806_Absff1.0e+11_Rayff1.0e+11_R160.00_all.root -o Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID806_Absff1.0e+11_Rayff1.0e+11_R160.00_all.txt
#
#./bin/splitReference_intoMaps -f Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID806_Absff1.0e+11_Rayff1.0e+11_R160.00_all.txt
#
##./bin/make_referenceMap_PMT-basedBins -f ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID807_Absff1.0e+11_Rayff1.0e+11_R180.00_all.root -o Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID807_Absff1.0e+11_Rayff1.0e+11_R180.00_all.txt
#
#./bin/splitReference_intoMaps -f Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID807_Absff1.0e+11_Rayff1.0e+11_R180.00_all.txt
#
##./bin/make_referenceMap_PMT-basedBins -f ../../WCSim/mPMTmapping/merged-data/hadded_wcsim_mPMTmapping_401nm_FileID808_Absff1.0e+11_Rayff1.0e+11_R210.00_all.root -o Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID808_Absff1.0e+11_Rayff1.0e+11_R210.00_all.txt
#
#./bin/splitReference_intoMaps -f Maps/AttenuationRef_PMT-basedBins/AttenuationRef_hadded_wcsim_mPMTmapping_401nm_FileID808_Absff1.0e+11_Rayff1.0e+11_R210.00_all.txt
#
