# Summer Dual Site Observations
Code for scheduling dual site observations from Birr and Onsala LOFAR stations.

It produces schedules in the iLiSA format and the RÉALTA format.

Exoplanets are queried from the NASA Exoplanet Archive and filtered to only target exoplanets discoverd by TESS. The list of pulsars and their properties are pulled from the ATNF catalogue.

Running ```shcedule.sh``` in a bash terminal will create a schedule using the default parameters - observe from 20:00 UTC on the next Tuesday to 04:00 UTC the next Wednesday.

In order to schedule your own custom parameters run ```schedule_options.sh``` in a bash terminal.

## Outputs
```schedule.sh``` creates a folder ```./outputs/``` which contains the schedules in both formats as .txt files. Also contains graphs of the optimal observation times for the targets as well as the actual observation times (in RA/LST).

```schedule_options.sh``` creates a folder ```./outputs_options/``` which contains the schedules in both formats as .txt files. Also contains graphs of the optimal observation times for the targets as well as the actual observation times (in RA/LST).

## Requirements
The scheduling script requires a number of python modules in order to successfully query the respective databases. They can all be installed using ```pip```.

Astroquery can be installed using
```
python -m pip install -U --pre astroquery
```

psrqpy is installed by
```
python -m pip install psrqpy
```

It also requires more common modules: ```numpy```, ```astropy```, ```matplotlib``` and ```datetime```.

## To do
- Add ability to automatically choose a flare star to observe.
