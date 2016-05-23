# Analysis_waveforms



The code reads out waveforms taken with an oscilloscope
and characterizes them in order to quantify the noise
of a SiPM.


To run here is an example:
./output -d /Users/Analysis_waveforms/ov_example/ -S /Users/Analysis_waveforms/config_file.txt -o /Users/Analysis_waveforms/Results/ -a

There is an already compiled version of the program in this repository,
but if desired the file "compile.sh" together with the cpp files
can be used to compile it. The program ROOT is necessary: -> root.cern.ch

There is a folder with waveforms examples to be analyzed.
The folder Results has the output of the analyzed waveforms
using the configure file also in found in this repository.
