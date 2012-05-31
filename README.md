SN-lines
========

Code for line ID

snlines.py
----------
This is the main code. Run by executing ./snlines.py _filename_, where _filename_ is a 2 column (wavelength, intensity) text file.

pickledata.py
-------------
This transforms a text file of line information like kurucz_cd23_cut.dat into a python dict object and pickles it. The main code loads the pickle file, not the .dat file.