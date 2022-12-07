# retinaCaImaging
Code repo for analysing calcium imaging data collected from the FLAME system at Newcastle University

Complete Guide to Code Installation and Use

This toolbox is designed to allow the user to extract calcium imaging data from the FN1 FLAME system. This toolbox also contains a number of analysis scripts to extract stimulus driven Calcium transient activity from imaging recordings. 

Software requirements:
1.	Download toolbox and add to Matlab path
2.	If you want to be able to use no rigid motion correction please clone the non-rigid motion correction toolkit (https://github.com/flatironinstitute/NoRMCorre) 
3.	Install an up to date version of FIJI (https://fiji.sc/) 
4.	Connect FIJI with matlab as explained here (http://bigwww.epfl.ch/sage/soft/mij/) NB, instead of using ij.jar, place the up to date version from your FIJI package (FIJI.app/jars), it will be named something like ij-1.52g.jar into the MATLAB folder.
5.	Install Cell Magic Wand into FIJI (https://www.maxplanckflorida.org/fitzpatricklab/software/cellMagicWand/) OR (https://github.com/GrimmSnark/Cell_Magic_Wand)
6.	You may need to increase your java heap size for FIJI and matlab to work with large images see (https://www.mathworks.com/matlabcentral/answers/92813-how-do-i-increase-the-heap-space-for-the-java-vm-in-matlab-6-0-r12-and-later-versions) NB use the java.opts method. 
7.	You will need to modify the "intializeMIJ.m" to your local FIJI path.
