This directory contains PU Reweighting files created on 4 November 2015.
The command used for their creation is the following:

./CalcPU.sh CMSSW_5_3_7_patch6 https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-208686_All_2012_pixelcorr_v2.txt crabJSON_SingleMu_22Jan2013.json true 60 60 3427

The directory where this command should be executed is:
/user/aolbrech/PUReweighting/

They correspond to the files with timestamp: 1446642523

The given TopTreeId (3427) is deceptive since it corresponds to the entire ABCD rereco sample.
The PileUp JSON file used is the latest one and also contains the pixelcorrections.

The main difference with respect to the distributions used by Gerrit and Stijn is the tail visible at high values.
However since this only seems to be relevant for log-scale, it will probably  not matter that much ...
In case strange results are obtained for the DATA-MC agreement, this is maybe the place to have a look!
