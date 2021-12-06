# RepeatingEQs
Scrach Scripts for idetifying Repeating Earthquakes

# simple description of tools 
A tentatively selection would be using long-term observations
a. recurrence pattern to evaluate the aseismic(creep) rate 
b. velocity/material property change (given the "identical" source location)
c. change in the source area (assuming the "identical" ray path) -- fault healing 

# waveform cross-correlation 
A few scripts for measuring the cross-correlation (or using obspy built-in tool)
check the following two commands: 
make
then type: sac_wfcc and wfcc

# coherence 
A modified script from K. Materna (2018, GRL). The difference between coherence and cross-correlation: 
coherence generally describes the similarity (or cc) at all different frequencies (hence, a vector)
cross-correlation is the similarity at a given frequency band 

# aided by relocation (hypoDD -- input: waveform CC differential travel times)
I tested the relocation using different starting locations (as expected, different starting locations 
generally converges)
