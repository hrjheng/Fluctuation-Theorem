#!/bin/bash

root -l -q -b Simulation_fluctuation.C\(0.6,0.1,0.1,120\)
root -l -q -b Simulation_fluctuation.C\(0.6,0.3,0.3,120\)
root -l -q -b Simulation_fluctuation.C\(0.6,0.5,0.5,120\)
root -l -q -b Simulation_fluctuation.C\(0.6,0.7,0.7,120\)
root -l -q -b Simulation_fluctuation.C\(0.6,0.9,0.9,120\)
root -l -q -b ana.C\(0.6,0.1,0.1,120\) 
root -l -q -b ana.C\(0.6,0.3,0.3,120\)
root -l -q -b ana.C\(0.6,0.5,0.5,120\)
root -l -q -b ana.C\(0.6,0.7,0.7,120\)
root -l -q -b ana.C\(0.6,0.9,0.9,120\)
    