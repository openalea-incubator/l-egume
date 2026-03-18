=========================
README for L-egume
=========================

This is L-egume model, a generic model of forage legume morphogenesis.

See :
* Louarn, G., Faverjon, L. (2018). A generic individual-based model to simulate morphogenesis, C–N acquisition and population dynamics in contrasting forage legumes. Annals of botany, 121(5), 875-896.
* Faverjon, L. (2018). Calibration et evaluation d’un modele individu-centre generique de morphogenese des legumineuses fourrageres – Application a la prediction des equilibres inter-specifiques dans des communautes prairiales  experimentales. PhD Thesis. Univ. Poitiers.

<img src="legume/logo_l-egume.png" alt="alt text" width="60%"/>

## 1. Getting Started

These instructions will get you a copy of *L-egume* up and running on your local 
machine.

### 1.1 Prerequisites

To install and use *L-egume*, you need first to install the dependencies.

*L-egume* has been tested on Windows 10 64bit.
 
#### 1.1.1 Install the dependencies on Windows 10 64 bit

1. Install Python 3.7 or 3.9 using Anaconda 

    * go to https://www.anaconda.com/download/ 
    * click on "64-Bit Graphical Installer", 
    * download "Anaconda3-2020.02-Windows-x86_64.exe" and install it selecting the following options:
        * install for all users,
        * default destination directory,
        * install all subfeatures, including subfeature "Add python.exe to Path".

		
2. Create and Activate a conda environment using  'Anaconda Prompt':
	* Open an 'Anaconda Prompt' console
	* Create a new environment (e.g. *envtest*) using the following command lines:
		```bash
		conda create -n envtest python=3.9 xlrd=2.0.1 numpy=1.20.3 scipy=1.7.3 pandas=1.3.4 openalea.lpy openalea.mtg alinea.caribu -c conda-forge -c fredboudon
		```
		
	* Activate the new environment using the following command line:
		```bash
		activate *envtest*
		```


	
### 1.2 Installing VGL submodels

__Note__: We suppose you already installed the dependencies for your operating system. Otherwise follow these [instructions](prerequisites "Prerequisites").


#### 1.2.1 Install *riri5* and *soil3ds* environmental models


To install *riri5* :

* open and activate the *envtest* conda environment with installed dependencies ,
* go to your local copy of project *riri5* (from https://github.com/glouarn/riri5),
* or get a copy of the latest model version (from a Git console: git clone https://github.com/glouarn/riri5.git)
* run command: 
	```bash
	python setup.py develop
	```


To install *soil3ds* :

* open and activate the *envtest* conda environment with installed dependencies ,
* go to your local copy of project *soil3ds* (from https://github.com/glouarn/soil3ds),
* or get a copy of the latest model version (from a Git console: git clone https://github.com/glouarn/soil3ds.git)
* run command: 
	```bash
	python setup.py develop
	```


#### 1.2.2 Install *L-egume* plant model in "develop" mode (recommended: will handle shortcuts)

Install *L-egume* in "develop" mode if you want to get *L-egume* installed and then 
be able to frequently edit the code and not have to re-install *L-egume* to have the 
changes to take effect immediately.

To install *L-egume* in "develop" mode:

* open and activate a conda environment with installed dependencies,
* go to your local copy of project *L-egume* (you can get the latest version from https://github.com/glouarn/l-egume/),
* or get a copy of the latest model version (from a Git console: git clone -b Develop https://github.com/glouarn/l-egume.git)
* run command: 
	```bash
	python setup.py develop
	```


### 1.3 Running


* open and activate the *envtest* conda environment with installed models

To run a simulation example, three options:

* 1. Run l-egume from the L-py GUI,
	 launch 'lpy' from the *envtest* conda environment 
	 open/load 'l-egume.lpy' file from l-egume folder,
	 Use Run or Animate button to launch a simulation from within L-py GUI
	 
  2. Run l-egume from the command line: 
		- default example:
		```bash
		python run_legume_usm.py
		```
		- run of a specific Unit of Simulation (USM):
		```bash
		python run_legume_usm.py -f 'usm_xlsfile' -i 'inputs_folder' -b 'usm_spreasheet_name' -u 'usmID' -o 'outputs_folder'
		```
		
  
  3. Run multiple simulations: see l-egume_batch.py in multisim folder for an example (require mutiprocessing)

See the user guide for a step by step explanation of how to set and run model *L-egume* (https://github.com/glouarn/TD_VGL).



[AFTER: TO BE COMPLETED!!!]



## 2. Reading the docs

To build the user and reference guides:

* install the model (see [Installation of the model](installing "Installing")), 
* open and activate the *envtest* conda environment
* to install sphinx, run command: 
	```bash
	conda install pytest sphinx sphinx-rtd-theme -c conda-forge
	```
* move to the *docs* folder within *l-egume* project
* run command:
	```bash
	make html
	```
* and direct your browser to file `docs/_build/html/index.html`.
* (To be done...`),



## 3. Testing

The test allows to verify that the model implementation accurately 
represents the developer’s conceptual description of the model and its solution.


To run the test :

* install the model (see [Installation of the model](installing "Installing")), 
* open a command line interpreter,
* go to the directory `test` of your local copy of the project,
* (To be done: and run this command: `python test_legume.py`).

## Built With

* [Python](http://www.python.org/), [NumPy](http://www.numpy.org/), [SciPy](http://www.scipy.org/), 
* [Sphinx](http://sphinx-doc.org/): building of the documentation, 



## Contact

For any question, send an email to <gaetan.louarn @ inrae.fr>.

## Versioning

We use a Git repository of OpenAlea on [GitHub](https://github.com/openalea-incubator/) for 
versioning: https://github.com/openalea-incubator/l-egume  
If you need an access to the current development version of the model, please send 
an email to <gaetan.louarn @ inrae.fr>.
For versionning, use a git client and get git clone git+git@github.com:openalea-incubator/l-egume.git SSH will is required

## Authors

**Gaetan LOUARN**, **Lucas FAVERJON** - see file [AUTHORS](AUTHORS) for details

## License

This project is licensed under the CeCILL-C License - see file [LICENSE](LICENSE) for details
