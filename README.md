# BPplus-Reservoir
Reservoir analysis for BPplus files (bpp_Res2).

#### Project Status: 
This is beta 4 version in development and hasn't been subjected to extensive testing

### Known bugs
* none at present!

## Project Description
The blood pressure waveform contains valuable information about the heart and circulation over and above simple measurement of systolic and diastolic blood pressure. The aim of this project is to use brachial blood pressure waveforms recorded by a BPplus device (https://www.uscom.com.au/products/bp/overview/) and perform reservoir analysis, pressure-only wave intensity analysis and estimate some other hemodynamic parameters from blood pressure waveforms. The general approach is similar to a method developed previously to perform the same analysis on radial blood pressure waveforms recorded with a Sphygmocor tonometer device (https://atcormedical.com/). 

The method for reservoir analysis is based on the approach used in doi: 10.1152/ajpheart.00875.2009. and the method for wave intensity analysis is based on https://doi.org/10.3389/fphys.2020.00550. 

The code is written in Matlab R2021a and aims to read BPplus files (\*.xml) obtained from version 1.0 (Cardioscope software version SW.R7.VME.037) to version 5.0 (software version 3.0.0.0) and exports results to Microsoft Excel. 

### Requirements
- Matlab R2021a or later versions (code may work on earlier versions but this hasn't been tested). 
- Microsoft Excel

### Getting started
Use of the code is described in the manual.

## Useful references
 - Davies JE, Baksi J, Francis DP, Hadjiloizou N, Whinnett ZI, Manisty CH, Aguado-Sierra J, Foale RA, Malik IS, Tyberg JV, Parker KH, Mayet J, Hughes AD. The arterial reservoir pressure increases with aging and is the major determinant of the aortic augmentation index. Am J Physiol Heart Circ Physiol. 2010 Feb;298(2):H580-6. doi: 10.1152/ajpheart.00875.2009. Epub 2009 Dec 11. PMID: 20008272
 - Alun Hughes, Chloe Park, Anenta Ramakrishnan, Jamil Mayet, Nish Chaturvedi, Kim Parker. Feasibility of estimation of aortic wave intensity using non-invasive pressure recordings in the absence of flow velocity in man. Front. Physiol., 2020; 11: 550. https://doi.org/10.3389/fphys.2020.00550

## Authors

* Alun Hughes (https://github.com/adh30)

## License

This project is licensed under the GNU General Public License v3.0 (https://github.com/adh30/BPplus-Reservoir/blob/Version-1-beta-4/LICENSE)

## Acknowledgments

* I am grateful to Kim Parker who wrote much of the early code for reservoir analysis
* I am grateful to Richard Scott for information about variables in the BP plus XML file. 
