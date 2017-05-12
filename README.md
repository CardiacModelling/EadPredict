# EadPredict - Early Afterdepolarisation Prediction

This project is an extension of Chaste that is intended to be used 
for simulation of drug action in cardiac electrophysiology models in 
combination with EADs.

## Prerequisites

Before using this code you will need to download and install Chaste's
dependencies and the Chaste source code itself.

Please see [Getting Started] for details of how to do this 
(follow instructions for "Development Code User" to keep up to date with the latest code, or a release version if you want longer-term stability).

Chaste version 3.4 is recommended for this project.

The bolt-on ApPredict project is also needed for the APD simulations. Please see the [ApPredict] github page for installation instructions.

## Installation

This repo must be cloned into
```sh
<chaste source directory>/projects/EadPredict
```
so that all the file paths can be picked up correctly (replacing ```<chaste source directory>``` with the place you have put the Chaste source code). Alternatively, you can put a sim link from the above folder to wherever you clone this repo.

The following instructions should do the cloning, note that this project pulls in a submodule from the Chaste/cellml repository, so cloning it requires the ```--recursive``` option:
```sh
$ cd <chaste source directory>/projects
$ git clone --recursive https://github.com/teraspawn/EadPredict.git
```

### Older git

N.B. on really old git versions (<1.6.5), `--recursive` doesn't work and you need to do:
```sh
$ cd <chaste source directory>/projects
$ git clone https://github.com/Chaste/EadPredict.git
$ cd EadPredict
$ git submodule init
$ git submodule update
```

[Getting Started]: <https://chaste.cs.ox.ac.uk/trac/wiki/GettingStarted>
[ApPredict]: <https://github.com/Chaste/ApPredict/releases>

### Simulations

To find EAD thresholds for the drug test pack, run

```sh
cd <chaste source directory>
scons projects/EadPredict/TestEADs.hpp
scons projects/EadPredict/TestControlEADs.hpp
```

Then run the simulation using the following command line options:

```
--model
	1 = shannon_wang_puglisi_weber_bers_2004
	2 = fink_noble_giles_model_2008
	3 = ten_tusscher_model_2006_epi
	4 = ten_tusscher_model_2006_endo
	5 = ten_tusscher_model_2006_M
	6 = ohara_rudy_2011
	7 = grandi_pasqualini_bers_2010_ss
--intervention
	1 = cytosolic_sodium_concentration
	2 = membrane_rapid_delayed_rectifier_potassium_current_conductance
	3 = membrane_fast_sodium_current_shift_inactivation
	4 = membrane_L_type_calcium_current_conductance
	5 = membrane_L_type_calcium_current_conductance_double_extracellular_potassium_concentration
	6 = membrane_persistent_sodium_current_conductance
	7 = membrane_rapid_delayed_rectifier_potassium_current_conductance_double_extracellular_potassium_concentration
	8 =  membrane_fast_sodium_current_shift_inactivation_double_extracellular_potassium_concentration
	9 = membrane_L_type_calcium_current_conductance_three_quarters_extracellular_potassium_concentration
	10 = membrane_rapid_delayed_rectifier_potassium_current_conductance_three_quarters_extracellular_potassium_concentration
	11 = membrane_fast_sodium_current_shift_inactivation_three_quarters_extracellular_potassium_concentration
	12 = membrane_fast_sodium_current_reduced_inactivation
--ll [lower limit]
--hl [higher limit]
--IP
	0 = do not create intervention profile 
	1 = do create intervention profile 
```

where "model" is the cell model to use, "intervention" is the strategy for inducing EADs, "ll" is the lower limit for the interval bisection protocol and "hl" is the upper limit (if in doubt, use 0 and 50), and "IP" creates an "intervention profile", which shows at which value of the intervention parameter an EAD is created over the whole range. For speed use "--IP 0".

To get the threshold values used for Table 2 in the paper, run: 

```sh
cd <chaste source directory>
./projects/EadPredict/build/debug/TestEADsRunner --model 6 --intervention 2 --ll 0 --hl 1 --IP 0
./projects/EadPredict/build/debug/TestEADsRunner --model 6 --intervention 3 --ll 0 --hl 30 --IP 0
./projects/EadPredict/build/debug/TestEADsRunner --model 6 --intervention 4 --ll 1 --hl 30 --IP 0
./projects/EadPredict/build/debug/TestControlEADsRunner --model 6 --intervention 2 --ll 0 --hl 1 --IP 0
./projects/EadPredict/build/debug/TestControlEADsRunner --model 6 --intervention 3 --ll 0 --hl 30 --IP 0
./projects/EadPredict/build/debug/TestControlEADsRunner --model 6 --intervention 4 --ll 1 --hl 30 --IP 0
```

To get the other metrics:

```sh
cd <chaste source directory>
scons projects/EadPredict/test/TestAPD90.hpp
scons projects/EadPredict/test/TestGrandiLancasterSobie.hpp
scons projects/EadPredict/test/TestOharaLancasterSobie.hpp
scons projects/EadPredict/test/TestCqinwardMetric.hpp
scons projects/EadPredict/test/TestControlAPDs.hpp
```

To collate the data into one file for analysis:

```sh
cp projects/EadPredict/collate_data.py $CHASTE_TEST_OUTPUT
cd $CHASTE_TEST_OUTPUT
python collate_data.py
```

To scale metrics for combination use the Matlab script `OutputCombinedMetrics.m`.

To classify data into categories use the Matlab scripts `ClassifyByCqinward.m` and `Compare_Classifiers.m`.

To do five-fold validation use the Matlab scripts `Cqinward_Five_Fold_Validation.m` and `Five_Fold_Validation.m`.

To create the dendrograms, run `OrderDendrogram.R` in R followed by `SideBySideDendrograms.m` in Matlab.

To create Figures 2 and 3, run the Matlab scripts `PlotEADClassificationExamples.m` and `EADComparisonFigure.m`.

