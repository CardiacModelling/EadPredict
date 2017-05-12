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
