=========================
README for L-egume
=========================

This is L-egume model, a generic model of forage legume morphogenesis.

See 
Louarn, G., & Faverjon, L. (2018). A generic individual-based model to simulate morphogenesis, C–N acquisition and population dynamics in contrasting forage legumes. Annals of botany, 121(5), 875-896.
Faverjon, L. (2018). Calibration et évaluation d’un modèle individu-centré générique de morphogenèse des légumineuses fourragères – Application à la prédiction des équilibres inter-spécifiques dans des communautés prairiales  expérimentales. PhD Thesis. Univ. Poitiers.



## 1. Getting Started

These instructions will get you a copy of *L-egume* up and running on your local 
machine.

### 1.1 Prerequisites

To install and use *L-egume*, you need first to install the dependencies.

*L-egume* has been tested on Windows 7/10 64bit.
 
#### 1.1.1 Install the dependencies on Windows 10 64 bit

1. Install Python  2.7 using Anaconda 

    * go to https://www.anaconda.com/download/ 
    * click on "64-Bit Graphical Installer", 
    * download "Anaconda2-2018.12-Windows-x86_64.exe" and install it selecting the following options:
        * install for all users,
        * default destination directory,
        * install all subfeatures, including subfeature "Add python.exe to Path".

		
2. Create and Activate a conda environment using  'Anaconda Prompt':
	* Open an 'Anaconda Prompt' console
	* Create a new environment (e.g. py27_64) using the following command lines:
		set CONDA_FORCE_32BIT=
		conda create -n py27_64 python=2.7
	* Activate the new environment using the following command line:
		activate py27_64

3. Install OpenAlea PlantGL:  

    * go to https://github.com/openalea/plantgl to check latest install recommendations for PlantGL
    * in your acivated conda environment, install successively : 
		conda install -c openalea openalea.core
		conda install -c openalea openalea.deploy
		conda install -c openalea qhull
		conda install -c openalea ann
		conda install -c openalea/label/unstable openalea.plantgl

4. Install OpenAlea L-py:  

    * go to https://github.com/openalea/lpy to check latest install recommendations for PlantGL
    * in your acivated conda environment, install successively : 
		conda install -c openalea openalea.lpy
		
	* To correct lpy console display bug:
		- go to envs\nom_environnement\library\bin
		- create a new file 'mylpy.bat'
		- write inside the following command line:
			python -c "from openalea.lpy.gui.lpystudio import main ; main()"
		- you can now launch lpy studio from your acivated conda environment using  “mylpy”

		
5. Install xlrd package:  
    * in your acivated conda environment, install successively : 
		conda install -c anaconda xlrd



	
### 1.2 Installing

__Note__: We suppose you already installed the dependencies for your operating system. Otherwise follow these [instructions](prerequisites "Prerequisites").

You can install *L-egume* either in "install" or "develop" mode.

#### 1.2.1 Install *L-egume* in "install" mode

Install *L-egume* in "install" mode if you're not going to develop, edit or debug 
it, i.e. you just want to used it as third party package.

To install *L-egume* in "end-user" mode:

* open and activate a conda environment with installed dependencies,
* go to your local copy of project *L-egume* (you can get the latest version from https://sourcesup.renater.fr/projects/l-egume/),
* run command: `python setup.py install --user`.

#### 1.2.2 Install *L-egume* in "develop" mode

Install *L-egume* in "develop" mode if you want to get *L-egume* installed and then 
be able to frequently edit the code and not have to re-install *L-egume* to have the 
changes to take effect immediately.

To install *L-egume* in "develop" mode:

* open and activate a conda environment with installed dependencies,
* go to your local copy of project *L-egume* (you can get the latest version from https://sourcesup.renater.fr/projects/l-egume/),
* run command: `python setup.py develop --user`.

### 1.3 Running

To run a simulation example, two options:

* 1. open Lpy platform,
	 load l-egume.lpy file from legume folder,
	 Use Run or animate button to launch a simulation
  2. Run l-egume from a python script (see l-egume_batch_mixture8.py in multisim folder for an example)

See the user guide for a step by step explanation of how to set and run model *L-egume*.



[AFTER: TO BE COMPLETED!!!]



## 2. Reading the docs

To build the user and reference guides:

* install the model (see [Installation of the model](installing "Installing")), 
* open a command line interpreter,
* go to the top directory of your local copy of the project,
* run this command: `python setup.py build_sphinx`,
* and direct your browser to file `doc/_build/html/index.html`.

## 3. Testing

The test allows to verify that the model implementation accurately 
represents the developer’s conceptual description of the model and its solution.

The test:

* runs the model on 200 steps,
* concatenates the outputs of the model in pandas dataframes,
* writes the outputs dataframes to CSV files,
* compares actual to expected outputs,
* raises an error if actual and expected outputs are not equal up to a given tolerance.     

To run the test :

* install the model (see [Installation of the model](installing "Installing")), 
* open a command line interpreter,
* go to the directory `test` of your local copy of the project,
* and run this command: `python test_lgrass.py`.

## Built With

* [Python](http://www.python.org/), [NumPy](http://www.numpy.org/), [SciPy](http://www.scipy.org/), 
* [Sphinx](http://sphinx-doc.org/): building of the documentation, 

## Contributing

First, send an email to <lgrass-request@groupes.renater.fr> to be added to the project.  

Then,
 
* check for open issues or open a fresh issue to start a discussion around a
  feature idea or a bug: https://sourcesup.renater.fr/tracker/?group_id=3957.
* If you feel uncomfortable or uncertain about an issue or your changes, feel
  free to email <lgrass@groupes.renater.fr>.

## Contact

For any question, send an email to <lgrass-request@groupes.renater.fr>.

## Versioning

We use a Git repository on [SourceSup](https://sourcesup.renater.fr) for 
versioning: https://sourcesup.renater.fr/projects/l-egume/.  
If you need an access to the current development version of the model, please send 
an email to <lgrass-request@groupes.renater.fr>.
For versionning, use a git client and get git clone git+ssh://git@git.renater.fr:2222/l-egume.git SSH will is required

## Authors

**Gaëtan LOUARN**, **Lucas FAVERJON** - see file [AUTHORS](AUTHORS) for details

## License

This project is licensed under the CeCILL-C License - see file [LICENSE](LICENSE) for details
