# PyKMOD
A Kurucz model atmosphere interpolator for Python!

This is essentially a re-writing of the original tool, [KMOD](http://hebe.as.utexas.edu/stools/), which was written in IDL, but modified such that the output files are immediately readable by the popular spectroscopic analysis program [MOOG](http://www.as.utexas.edu/~chris/moog.html).
Python being a more flexible language, PyKMOD can be run either by importing it into a Python script and using the function pykmod(), or it can be run directly from the command line via your Python interpreter.

Included for convenience are ATLAS9 model atmospheres which were downloaded from [this page](http://research.iac.es/proyecto/ATLAS-APOGEE/) on the ATLAS-APOGEE website. PyKMOD can be adapted for reading any model atmosphere file by modifying its rd_kmod function.

(Should be run with the latest versions of SciPy and Numpy)

### Installation
Almost none required. Simply download the files provided and extract them to a directory of your choice. If you want to import it as a python module, you can also copy the files to your Python lib folder to make this easier.

### Using PyKMOD
#### Command Line
The command line syntax is as follows:
```
python pykmod.py (Teff (K)) (logg ([cm/s^2])) (vmicro (km/s)) ([Fe/H] (dex)) [output file path]
```
The output file path is optional. If none is provided, it will default to a file named "modelatmosphere.txt" in the PyKMOD directory.

#### As Python Module
```python
from pyKMOD.pykmod import pykmod

pykmod(Teff,logg,vmicro,FeH,outfile)
```
In this case an output file path is required

