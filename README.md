AnomalousCouplings
==================
Three directories defined within this repository:

*** MadWeightOutput: 
Contains all the output files obtained when running the AnomalousCouplingsTreeCreator file (which should still be updated to write this file to this directory)
Also the python script developed to perform the matching between the generator and the reconstructed events can be found here.

(
Maybe it is a possibility to first create the output file for the reconstructed events and only then rerun the generator events. 
This can avoid having to much generator event information in the .lhco output which will only slow down MadWeight!
--> Maybe add a boolean allowing the entire generator events collection to be used when only generator info is needed!

Also try to get a .lhco file for every different file within the config file.
This will automatically reduce the number of events which go into one .lhco file.
)

*** ProducerTestFiles:
Contains the files developed when doing the initial tests for the generator event information.
This was first obtained without any TopTree information, but using TFiles in the TopTreeProducer/test directory!
The used macro is based on Macro.c.

*** config:
Contains the config file which will be used for the analysis.
(The AnomalousCouplingsTreeCreator should still be updated to find the config file in the correct directory!)

°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Command to run the analyzer:

g++ -m64 -g -L /user/aolbrech/lib -I ../ -l TopTreeAnaContent53 -l TopTreeAna53 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` AnomalousCouplingsTreeCreator.cc -o AnomalousCouplingsTreeCreator

source /jefmount_mnt/jefmount/cmss/slc5_amd64_gcc462/external/dcap/2.47.5.0-cms/etc/profile.d/init.sh

./AnomalousCouplingsTreeCreator

----------------
push command which need to be used: git push origin master
--> Important to make sure that origin is defined correctly in the .git/config file
Should be as followed:

[core]
	repositoryformatversion = 0
	filemode = true
	bare = false
	logallrefupdates = true
[remote "origin"]
	fetch = +refs/heads/*:refs/remotes/origin/*
	url = git@github.com:/TopBrussels/AnomalousCouplings.git
[branch "master"]
	remote = origin
	merge = refs/heads/master

Git command to merge development branch 'test' with master branch:
 * git checkout master
 * git pull origin master
 * git merge test
 * git push origin master
