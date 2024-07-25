# BPplus-Reservoir
Reservoir analysis for BPplus files *bpp_Res2*.
Matlab scripts to read BPplus files (\*.xml) and perform reservoir analysis, pressure-only wave intensity analysis and estimate some other hemodynamic parameters from blood pressure waveforms.

**This is beta 4 version and hasn't been subjected to extensive testing**

The method for wave intensity analysis is based on: Alun Hughes, Chloe Park, Anenta Ramakrishnan, Jamil Mayet, Nish Chaturvedi, Kim Parker.
Feasibility of estimation of aortic wave intensity using non-invasive pressure recordings in the absence of flow velocity in man.
Front. Physiol., 2020; 11: 550. https://doi.org/10.3389/fphys.2020.00550

I am grateful to Richard Scott for information about variables in the BP plus XML file and some code suggestions.

Known bugs
* identification of end-systole unreliable. SEVR estimates dubious when this fails. 
