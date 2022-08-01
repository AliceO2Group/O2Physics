# O2Physics Tutorial for CF
Table of contents: <br>
1. [Prerequisites](#prerequisites) <br>
2. [Create your task](#create-your-task-in-o2physics) <br>
3. [Tutorial checkpoints](#tutorial-checkpoints) <br>
_________________
_______________
# Prerequisites
## How to install O2Physics Framework?
1. Build prerequisities for your operating system, configure aliBuild (but don't build packages!) - [here](https://alice-doc.github.io/alice-analysis-tutorial/building/custom.html).<br>
2. Get a GRID certificate - [here](https://alice-doc.github.io/alice-analysis-tutorial/start/cert.html) <br>
3. Follow instructions on [ALICE O2 documentation](https://aliceo2group.github.io/analysis-framework/docs/installing/) - Prepare your source code, check prerequisites, build and rebuild and check if everything works! <br>

Keep in mind that it might take quite a long time.
______________
## How to update ALICE O2Physics?

To update the software, you need to do the following: <br>
1. in `alice/alidist`, type `git pull` <br>
2. in `alice/O2Physics`, type `git pull --rebase` <br>
3. in `alice/`, type `aliBuild build O2Physics --defaults o2`. Here you can add `--debug`  (or `-d`) for more information.<br>

You need to update it frequently, because O2 is constantly evolving. Keep in mind that it might take quite a long time.
_________
## How to obtain data files for the analysis?

<!-- markdown-link-check-disable-next-line -->
This tutorial is made and tested for pilot-beam data (reconstruction pass4). You can download the reference `AO2D.root` file from [here](https://alimonitor.cern.ch/catalogue/index.jsp?path=%2Falice%2Fdata%2F2021%2FOCT%2F505669%2Fapass4%2FAOD#/alice/data/2021/OCT/505669/apass4/AOD/006). We use file `006`. Be aware that due to the file size (around 400 MB), it could take a while.

<!-- markdown-link-check-disable-next-line -->
A lot of Run 2 data converted into the Run 3 format can be found [here](https://alimonitor.cern.ch/trains/train.jsp?train_id=132). A detailed description of how to download converted Run 2 data is reported in the [offifial documentation](https://aliceo2group.github.io/analysis-framework/docs/download/)

In summary, to download a bunch of data, do the following: <br>
1. scroll all the way down, choose the train number you are interested in and click on run number. <br>
2. A pop-up window appears: click on the **Test Results**. <br>
3. Scroll down, find the **Full Train** option and then click on **output**. <br>
4. Search for the `AO2D.root` file: click on it to download it. <br>

Be aware that it might take a while, because the size is around 1 GB file.<br>
<!-- markdown-link-check-disable-next-line -->
Information about trains (job details) can be found [here](https://alimonitor.cern.ch/job_details.jsp).

____________________________
____________________________
# Create your task in O2Physics

There are many examples of simple introductory tasks in O2Physics that can be found in `O2Physics/Tutorials`.
All the information about how to write a task can be found in the [official documentation](https://aliceo2group.github.io/analysis-framework/docs/tutorials/).

We hope that this hands-on tutorial will help you learning to write your analysis tasks starting from scratch.

_________________
## A O2Physics task in summary

A task is a C++ `struct`, which must contain a `process` method. The `process` is the equivalent of what in `AliPhysics` was `UserExec`. It is used to access analysis objects, such as collisions, tracks, etc.

In addition, the task can also have an `init` funcion, which is used to initialise analysis-objects such as histogram. It is the equivalent of `UserCreateOutputObjects` in AliPhysics.

In the same source file where the `struct` is defined, one has to add a `defineDataProcessing` method, which is necessary to include the task to the workflow, assigning it a specific name.

The task must be added to the `CMakeLists.txt` file in the same directory. For our tutorial, the first task has name `cf-tutorial-0` and it is defined in the source file `CFTutorialTask0.cxx`. The corresponfing lines in the `CMakeLists.txt` file are:

```
o2physics_add_dpl_workflow(cf-tutorial-0
                    SOURCES CFTutorialTask0.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGCFCore
                    COMPONENT_NAME Analysis)
```
Everything will be (hopefully) clear in the hands-on session.

_____________
## How to build a task?

The safest way to build your task is to use `aliBuild`, following the same instructions reported in the section [_How to update O2Physics?_](#how-to-update-alice-o2physics)

If you want to save time and opt for a more advance method, read the following section

 **Advanced method**

Go to `~/alice/sw/BUILD/O2Physics-latest-master/O2Physics`.
<details><summary>If you have properly installed direnv, as soon as you enter a lot of text should appear and it should look like this:</summary>
<p>

![image](https://user-images.githubusercontent.com/87480906/162129203-4a4b833b-fefc-48c6-9229-908354cf0620.png)

</p>
</details>

In `~/alice/sw/BUILD/O2Physics-latest-master/O2Physics` enter ninja and O2Physics environment (`alienv load ninja/latest O2Physics/latest`). Then build your task using `ninja stage/bin/<your-analysis-file>`. If you don't know what you should type in `<your-analysis-file>` place, open `CMakeList.txt` and see how your analysis is called there. <br>
Keep in mind that if you add a new file or modify CMakeList you need to use `cmake .`
Then, after building part, copy this builded file to directory with AOD file (`cp stage/bin/<your-analysis-file>`). You can skip this step (copying) when you use `ninja install` instead of `ninja` or you can use alibuild.
____________________
## How to run your code?

In your directory with the AOD file, do the following: <br>
1. enter O2Physics envinronment, e.g with `alienv enter O2Physics/latest`
2. You need to connect to the GRID to acces remote configuration files (CCDB): type `alien.py` or `alien-token-init`. You will be asked for you password.<br>
3. Run your code as:
`./o2-analysistutorial-simple-analysis --aod-file <aod_file_name> [-b]`. Option `-b` stops the GUI from showing (batch mode). <br>

It is very likely that your task has dependencies: you have to run also the corresponding tasks (_helper tasks_).
For example, for our first task `cf-tutorial-task0` you would need to run:

```
o2-analysis-timestamp | \
o2-analysis-track-propagation | \
o2-analysis-cf-cf-tutorial-0 --aod-file AO2D.root
```

**Tip**: create a bash script with all the commands that you need!

**Tip**: tasks have parameters that can be provided via command line. For example:
```
o2-analysis-mytask --aod-file <aod_file_name> [--<par_name> <par_value>]
```
When you have many parameters, it is better to provide them in a configuration file in json format:
```
o2-analysis-mytask --configuration json://<config_file_name>
```
In the following, we will use configuration files.

### Possible errors
<details><summary>Couldn't get TTree.</summary>
<p>

Sometimes you may get an error, for example:
 ```c
 [ERROR] Exception caught: Couldn't get TTree "DF_2853960297589372650/O2v0dataext from <your AOD file>"
 ```
It means that `v0dataext` couldn't be found in your AOD file or it has not been produced by any attached helper task.
So now you know which table is missing. There are two paths you can choose; <br>
1. **EASIEST**: follow the instructions reported in  [ALICE O2 documentation- Tree not found](https://aliceo2group.github.io/analysis-framework/docs/troubleshooting/treenotfound.html). After entering [ALICE O2 documentation - helper task tables](https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html), you search the table that you are missing and you will find the task you need to attach as dependency.<br>
2. More difficult: enter directory `alice/O2Physics` and look for the missing table. </br>
For example, if you are missing `pidtofka` table, you should type:
```c
grep -rnw . -e "PIDTOFKA"
grep -rnw . -e "PIDTOFKa"
grep -rnw . -e "pidTOFKa"
```
You should see a list of files where this expression occurs: you need to find out the name of the file which produces this table. After you find it, you need to go to the corresponding `CMakeLists.txt` file. There, you will find the name of the task you are looking for. You can see that it has many sections but they all look alike. For example:
```c
o2physics_add_dpl_workflow(<task-name>
                    SOURCES <C++ file>
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)
```
So you need to find the **task-name** which will correspond to the **C++ file** you've found using grep ;) <br>


</p>
</details>
<details><summary>Fatal errors with non-existing CCDB entries.</summary>
<p>
While trying to run various O2 tasks/tutorials, you can run into fatal errors with non-existing CCDB entries for specific data. <br>
The reason is often the incorrect configuration of the paramenters. For example, for some executables such as `o2-analysis-timestamp` duging the option `--isRun2MC` could solve the problem. Or, simialary, using the oprion `--isMC` for o2-analysis-event-selection. <br>
For example:

```c
o2-analysis-mm-dndeta --aod-file AODmc2.root --aod-memory-rate-limit 100000000000000 -b | \
o2-analysis-pid-tpc -b | \
o2-analysis-pid-tof -b | \
o2-analysis-trackselection -b | \
o2-analysis-trackextension -b | \
o2-analysis-event-selection --isMC -b | \
o2-analysis-timestamp --isRun2MC -b
```

</p>
</details>

____________________
## How to obtain results of the analysis?
Running command `./o2-analysistutorial-simple-analysis --aod-file <aod_file_name> -b` (and others of that type) creates .root file and .json file. <br>
To enter AnalysisResults.root file enter O2Physics environment (`enter O2Physics/latest`) and type: `root -l` and then `new TBrowser`. (`-l` is optional but it runs root without additional information about it).<br>
You can change histograms manually or you can write macros and run them on .root file. If you've never written any root macros consider [_this repository_](https://github.com/zchochul/KADDre). You can find more on root and how to use it here -> [_root.cern/manual/_](https://root.cern/manual/first_steps_with_root/). <br>
_____________________________
_____________________________
# Tutorial checkpoints
## First task
 We will be using [CFTutorialTask0.cxx](https://github.com/CF-tutorials/O2Physics/blob/tutorial/PWGCF/Tutorial/CFTutorialTask0.cxx). It's an example task illustrating how to create histograms and fill them with basic information. We will also apply basic event selection.<br>

 Here you can see a structure of a typical analysis task in the O2:<br>
 ```c
 struct myTask{
    Partition<something> partition_name;
    Filter filter_name;
    Configurable<something> name{};
    void init(){
    }
    void process(){
    }
 }
 WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
    WorkflowSpec workflow{};
    return workflow;
 }
 ```
Each part will be described in greater detail further on during this tutorial.<br>

The first element of the task is `init` (equivalent to _UserCreateOutputObjects_ in _AliPhysics_). It looks like this:

  ```c
   void init(o2::framework::InitContext&)
  {
    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    // Add histograms to histogram manager (as in the output object in AliPhysics)
    histos.add("hZvtx_before_sel", ";Z (cm)", kTH1F, {vtxZAxis});
    histos.add("hZvtx_after_sel", ";Z (cm)", kTH1F, {vtxZAxis});
    histos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.5, 4.}});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
  }
 ```
`init` is used to initialise objects, e.g. histograms, which will be used later in the process.

The second element is `process`, equivalent to the _UserExec_ in _AliPhysics_. It looks like this:

```c
    void process(aod::Collision const& coll, aod::Tracks const& inputTracks)
  {
    // Performing the event selection
    histos.fill(HIST("hZvtx_before_sel"), coll.posZ());
    if (fabs(coll.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());

    for (auto track : inputTracks) { // Loop over tracks
      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
    }
  }
 ```

The arguments of `process` are the data _tables_ that we analyse, i.e. the tables to which we _subscribe_. In the snippet above, the subscripion is to the tables of collisions and tracks.

Finally, `WorkflowSpec` is necessay to add the task to the DPL workflow . And it looks like:

 ```c
   WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<CFTutorialTask0>(cfgc)};
  return workflow;
}
```

## Second task

We will be using [CFTutorialTask1.cxx](https://github.com/CF-tutorials/O2Physics/blob/tutorial/PWGCF/Tutorial/CFTutorialTask1.cxx). In this part of the tutorial, we will focus on how to access information from different tables. Sets of information stored in different tables can be put together using `Join`:

```c
namespace o2::aod
{
using MyCollision = soa::Join<aod::Collisions,
                              aod::EvSels,
                              aod::Mults>::iterator;
using MyTracks = soa::Join<aod::FullTracks,
                           aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                           aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe>;
} // namespace o2::aod
 ```
The data model provides some predifined joins, suche as `FullTracks = soa::Join<Tracks, TracksExtra>`. The complete list of predefined joins can be found in [ALICE O2 documentation - The Data Model, Joins and iterators](https://aliceo2group.github.io/analysis-framework/docs/datamodel/joinsAndIterators.html).

Joined tables can be normally used as agruments of process.

```c
 // CFTutorialTask1.cxx version
 void process(aod::MyCollision const& coll, aod::MyTracks const& inputTracks){...}
 // CFTutorialTask0.cxx version
 void process(aod::Collision const& coll, aod::Tracks const& inputTracks){...}
```

## Third task

We will be using [CFTutorialTask2.cxx](https://github.com/CF-tutorials/O2Physics/blob/tutorial/PWGCF/Tutorial/CFTutorialTask2.cxx). In this part of the tutorial, we will focus on how to use **Configurables** and **Filters**. <br>

A `Configurable` is a parameter of the task, which can be set from command line or in the dedicated space in Hyperloop.
In the following, an example for the configurable used for the selection on the vertex z coordinate:
```c
Configurable<float> ConfZvtxCut{"ConfZvtxCut", 10, "Z vtx cut"};
```
In the task, tha value of the configurable can be accesses simply by calling the Configurable itself.

A `Filter` is used to select the rows of a table which satisfy particular requirements. For example:
```c
Filter trackFilter = (nabs(aod::track::eta) < ConfEtaCut) && (aod::track::pt > ConfMinPtCut) && (aod::track::pt < ConfMaxPtCut);
```
Note that configurables can be used in the definition of a filters.

## Fourth task

We will be using [CFTutorialTask3.cxx](https://github.com/CF-tutorials/O2Physics/blob/tutorial/PWGCF/Tutorial/CFTutorialTask3.cxx). In this part of the tutorial, we will focus on how to create and use **partitions**.<br>

A `Partition` is a subset of a given table that satisfies particular requirements. In the following, a partition obtained from the table _Tracks_:

```c
Partition<o2::aod::Tracks> positive = (nabs(aod::track::eta) < ConfEtaCut) && (aod::track::pt > ConfMinPtCut) && (aod::track::pt < ConfMaxPtCut) && (aod::track::signed1Pt > ConfChargeCut);
```
Note that configurables can be used in the definition of a partitions.

Partition can be used in the following way:
 ```c
void process(MyFilteredCollision const& coll, o2::aod::Tracks const& tracks)
{

  ...

  auto groupPositive = positive->sliceByCached(aod::track::collisionId, coll.globalIndex());

  for (auto track : groupPositive) {
    histos.fill(HIST("hChargePos"), track.sign());
  }

  ...
}
```
**WARNING**: using a `Partition` of the table `Tracks`, you are considering **ALL** the tracks (which satisfy the requirements) of **ALL** the collisions. In order to select the elements corresponding to the relevant collision, one must use `sliceByCached`, as shown in the example above.

</p>
</details>

## Fifth task

We will be using [CFTutorialTask4.cxx](https://github.com/CF-tutorials/O2Physics/blob/tutorial/PWGCF/Tutorial/CFTutorialTask4.cxx). In this part of the tutorial, we will focus on how to use and combine elements of different **partitions**. In the following, you can see the invariant mass of two pions, a negative one and a positive one, belonging to two different partitions.

`combinations` returns a pair of elements taken from different partitions. `CombinationsFullIndexPolicy` is needed to combine all the element of the first set with all the element of the second set.

```c
for (auto& [pos, neg] : combinations(soa::CombinationsFullIndexPolicy(groupPositive, groupNegative))) {

      if (fabs(pos.tpcNSigmaPi()) > 3 or fabs(neg.tpcNSigmaPi()) > 3) {
        continue;
      }

      TLorentzVector posVec;
      posVec.SetPtEtaPhiM(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassPionCharged);
      TLorentzVector negVec;
      negVec.SetPtEtaPhiM(neg.pt(), neg.eta(), neg.phi(), o2::constants::physics::MassPionCharged);

      TLorentzVector sumVec(posVec);
      sumVec += negVec;
      histos.fill(HIST("hInvariantMass"), sumVec.M());
    }
```


## Sixth task

We will be using [CFTutorialTask5.cxx](https://github.com/CF-tutorials/O2Physics/blob/tutorial/PWGCF/Tutorial/CFTutorialTask5.cxx). In this part of the tutorial, we will show how to use multiple processes in the same task with **process switches** and how to access and combine elements of different **partitions** and **different events**.

With diss task, we compute the invariant mass of two pions, one positive and one negative.
We have two separate process functions: <br>

1. `processSameEvent` computes the invariant mass for two pions produced in the same collision.

2. `processMixed` computes the invariant mass for two pions produced in the different collisions.

```c
void processSameEvent(){}
PROCESS_SWITCH(<name-of-your-structure>, processSameEvent, "Enable processing same event", true);

void processMixed(){}
PROCESS_SWITCH(<name-of-your-structure>, processMixed, "Enable processing mixed event", true);
```

The process functions can be activated using `PROCESS_SWITCH`, setting the last parameter as `true`. This parameter can be configured via command line, json file or in the hyperloop interface.
