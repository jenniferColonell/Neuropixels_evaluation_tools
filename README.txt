This is a set of MATLAB tools for some basic evaluation of in vivo data from Neuropixels 1.0 or Neuropixels 2.0 probes.

There are two MATLAB scripts to run:

detect_merge_peaks.m
	-takes as input a .bin file from SpikeGLX and a tab delimited file of the site coordinates
	-calculates median absolute deviation for each channel in the data
	-thresholds the data (use a constant threshold for running probe comparisons)
	-merges peaks that occur close in space and time into "events" using the JRClust peak_merge algorithm
	-histograms the maximum negative amplitudes for each channel
	-stores all these results in a MATLAB matrix file, <input filename>_peaks.mat

allSpike_dist.m
	-takes as input the _peaks.mat file and the XY coordinate file
	-rebins the amplitude histograms to a specified size (2.34375 uV, the size of 1-bit in NP 1.0 data)
	-outputs by-channel stats in <input filename>_chan.txt
	-outputs summary stats about the amplitude distribution, estimated rms, etc, in <input filename>_sum.txt
	-outputs the histogram of amplitudes in <input filename>_hist.txt

To run these tools:
	
(1) If you have long recordings, use the SpikeGLX offline viewer export function to make a 10 minute subset.

(2) Filter the data (and concatentate trials) using the CatGT tool. Typical parameters are: -aphipass=300 -aplopass=9000 -gbldmx -gfix=5,8

(3) Set parameters in detect_merge_peaks.m, specifically dataType. If your recording only has a subset of the channels, set nchan and dataChan accordingly. Also, note that the XY coordinate file needs to match your recording. Run detect_merge_peaks.m.

(4) Set parameters in allSpike_dist.m. Set dataType to 0 for NP 1.0 or 1 for NP2.0 data. When analyzing NP1.0 data to compare to NP2.0 data, set zMax to 3000 microns. In exChan, list reference and "dead" channels which will have abnormally low rms. Run allSpike_dist.m. Check the "est RMS" plot for channels with anomalously low signal, if there are any, add those to the exChan list and rerun.



