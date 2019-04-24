# A finite element model for calcium reaction-diffusion in a cardiomyocyte 

This finite element model simulates reaction-diffusion in an eight-sarcomere representation of
a cardiomyocyte during the rising phase (first 30ms) of the calcium transient. This model builds on a previous half-sarcomere model (see [Rajagopal
et al.
2015](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004417)
and [accompanying code](https://github.com/CellSMB/cardiac_ecc)). Calcium (Ca2+)
is released from ryanodine receptor (RyR) ion channels into the intracellular
volume with an exponential rise and decay model specified in the CellML model
found in the input directory. RyR cluster locations, release times (RyRTimelag),
and densities (RyRDensity) were simulated for each z-disk location using the
[RyR-Simulator](https://github.com/CellSMB/RyR-simulator). 

The finite element model solves the
reaction-diffusion mechanics of diffusing Ca2+ in several fields: calcium (Ca),
fluorophore (F), fluorophore-bound calcium (FCa), calmodulin (CaM),
calmodulin-bound calcium (CaMCa), adenosine triphosphate (ATP), ATP-bound
calcium (ATPCa), troponin-c (CaTnC).

Two sets of parameters may be explored with the available supporting files:
* The effect of RyR cluster spacing in low (51 clusters per z-disk) or high (123
clusters per z-disk) cluster density configurations. This is controlled with the
'lowRyRDensity' boolean parameter.
* The effect of mitochondria acting as barriers to diffusion. This is controlled
  with the 'mitochondria' boolean parameter.

## Usage

This model is written in Python 3. It relies on the Iron library from the OpenCMISS mathematical modelling
environment to create, constrain, and solve the finite element model. See the
[OpenCMISS website](https://www.opencmiss.org) and [GitHub repository](https://github.com/OpenCMISS/iron)
for details on installation, usage, and source code. 

The main script checks for the expected mesh files in the input folder and
downloads them from online repositories using the Python
[Requests](http://docs.python-requests.org/en/master/) module if they are not
yet present. This works around tracking large binary datasets with git. To
run this program without internet access, download the mesh files from the URLs
indicated in the resourceLinks dictionary in the main script.

Python modules used in this project:
* [OpenCMISS-Iron](https://www.opencmiss.org/)
* [NumPy](https://www.numpy.org/)
* [SciPy](https://www.scipy.org/)
* [pandas](https://pandas.pydata.org/)
* [Requests](http://docs.python-requests.org/en/master/)
* [gc](https://docs.python.org/3/library/gc.html)
* [time](https://docs.python.org/3/library/time.html)
* [os](https://docs.python.org/3/library/os.html)

## Authors

* **David Ladd** - [dladd](https://github.com/dladd)

## Contributor(s)

* **Vijay Rajagopal** - [vraj004](https://github.com/vraj004)

## License

This project is licensed under the Apache 2.0 License - see the LICENSE file for details.

## Acknowledgments

* [Systems Biology Laboratory](https://systemsbiologylaboratory.org.au/) at the
  University of Melbourne
* [Centre of Excellence in Convergent Bio-Nano Science and Technology](https://www.cbns.org.au/)
* [Cell Structure and Mechanobiology Group](https://cellularsmb.org/) at the
  University of Melbourne
