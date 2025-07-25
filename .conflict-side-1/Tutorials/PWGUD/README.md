# UD tutorial, April 27th, 2023
- from 09:00 to 12:30 (CERN time)
- zoom link: [https://cern.zoom.us/j/63176005388?pwd=MjNVTUlrTVJwMGcyZ25xeUxVanVndz09](https://cern.zoom.us/j/63176005388?pwd=MjNVTUlrTVJwMGcyZ25xeUxVanVndz09)

## 1. Get started

For this tutorial it is assumed that you have a working installation of a recent version of O2Physics and that you know how to modify, compile and run c++ code in the O2 framework. An introduction into this topic can be found [here](https://indico.cern.ch/event/1267433/contributions/5359315/attachments/2629584/4547945/Git,%20GitHub%20and%20aliBuild%20intro.pdf).

To get started with this UD tutorial move into a working directory. I recommend to create a new directory for this, which can be deleted again, when all work is done.

Now download [UDTutorials_0423.tar.gz](https://cernbox.cern.ch/s/1DbjlSoVyQ7RlCE) (by clicking on the link) into your working directory. The tar-file contains a README and a few other files, which are needed to run the tutorial examples.

Unpack the tar file with
> \> tar -xvzf UDTutorials_0423.tar.gz

This will create the directory ./UDTutorials which includes the files

- .getDerivedData
- plotHistograms.C
- README
- UDTutorialsConfig.json

To continue, change into directory ./UDTutorials.
> \> cd UDTutorials

The tutorial examples will need skimmed data as input. What *skimmed data* means and how it is produced is explained during the tutorial session.
To download the data - roughly 900 MB - to your computer you first need to load the O2Physics environment. Then use the scripts .getDerivedData to retrieve the data

> \> ./.getDerivedData

The skimmed data is copied into directory ./derivedData and the downloaded files will be listed in ./rootfiles.txt.

When this is done you are ready to run the tutorial examples!


## 2. Running the tutorial examples

There are three example tasks which all can be found in the O2Physics repository at [O2Physics/Tutorials/PWGUD](https://github.com/AliceO2Group/O2Physics/tree/master/Tutorials/PWGUD).

- [UDTutorial_01.cxx](https://github.com/AliceO2Group/O2Physics/tree/master/Tutorials/PWGUD/UDTutorial_01.cxx)
- [UDTutorial_02a.cxx](https://github.com/AliceO2Group/O2Physics/tree/master/Tutorials/PWGUD/UDTutorial_02a.cxx)
- [UDTutorial_02b.cxx](https://github.com/AliceO2Group/O2Physics/tree/master/Tutorials/PWGUD/UDTutorial_02b.cxx)

During the tutorial session we will talk about the tasks and implement some modifications. However, to test whether your
installation is working fine you can run the first example as is by executing the following commands

> \> copts="--configuration json://UDTutorialsConfig.json -b"

> \> o2-analysis-ud-udtutorial-01 $copts > UDTutorial01.log

> \> mv AnalysisResults.root AR_UDTutorial01.root

> \> root -b -l -q 'plotHistograms.C("AR_UDTutorial01.root","LHC22q, UD Tutorials")';  mv UDTutorials.ps UDTutorial01.ps

This should run the task on the downloaded skimmed data, produce the result file AR_UDTutorial01.root with some histograms, and finally create UDTutorial01.ps with plots of the histograms.

Running the task over all data files can take a while. On my not so new laptop, it takes around 15 minutes. For testing you can reduce this time by removing some file names from rootfiles.txt. But don't forget to add them later again!
