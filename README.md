# CryoGrid
This is the modular version of CryoGrid which is based on an object-oriented programming paradigm.

## Get started as a first time user

### Getting the code and examples by download

Here we describe how to set up the code to run an initial test model as a first time user with no knowledge of GIT.
(You may also retrieve the code using git, see section [below](#git_download))

#### Get the test example

First download the test example from the separate github repository [CryoGrid/CryoGridExamples](https://github.com/CryoGrid/CryoGridExamples/tree/develop/main). Make sure to retrive the `develop/main` branch. Follow download instructions given in the CryoGridExamples repository.

Unzip the the contents to your preferred folder, f.ex. `c:\my_matlab_code\cryogrid\`.
You will now have the following folder structure: `c:\my_matlab_code\cryogrid\CryoGridExamples-develop-main`, which will contain example run files and model definitions.
Rename this folder:

```
rename c:\my_matlab_code\cryogrid\CryoGridExamples-develop-main c:\my_matlab_code\cryogrid\CryoGridExamples
```

#### Get the main CryoGrid code

Download the main CryoGrid code as zip-file. From the main GitHub repository page, click the green `code` button and choose `Download Zip`. See screenshot below.
(In the upper left hand corner, make sure the `develop/main` branch is selected).

![image-20201203092734144](./readme_im1.png)

Now unzip the contents of the zip file  to the folder `c:\my_matlab_code\cryogrid\CryoGridExamples`.
You will now have the following folder structure: `c:\my_matlab_code\cryogrid\CryoGridExamples\CryoGrid-develop-main`
Rename this folder:
```
rename c:\my_matlab_code\cryogrid\CryoGridExamples\CryoGrid-develop-main c:\my_matlab_code\cryogrid\CryoGridExamples\CryoGrid
```

#### <a name="run_model"></a>Run the test model 

Open MatLab and change the path to `c:\my_matlab_code\cryogrid\CryoGridExamples\`. This can be done either using the MatLab path selector dialog, or by typing in the command window:

`cd c:\my_matlab_code\cryogrid\CryoGridExamples\`

Change the paths according to your actual install directory you chose for the code.

To run the example model, in the MatLab terminal run the file run_CG.m file by typing:

`run_CG`

The code will start running. It will produce one ouput file per year, which is written to the disc at a specific date (defined in the parameter Excel file `c:\my_matlab_code\cryogrid\CryoGridExamples\results\test\test.xlsx`). You can stop the code any time after a full year has been calculated (to ensure you have an output file written to disc).

The first output file written to disc will have the name `c:\my_matlab_code\cryogrid\CryoGridExamples\results\test\test_19800901.mat`.

#### View the output of the model

To plot the results, change the path to `c:\my_matlab_code\cryogrid\CryoGridExamples\CryoGrid\analyze_display\`:

```
cd c:\my_matlab_code\cryogrid\CryoGridExamples\CryoGrid\analyze_display\
load(`c:\my_matlab_code\cryogrid\CryoGridExamples\results\test\test_19810901.mat`)
read_display_out()
```

Plots will be generated for several parameters. Not all of them are meaningful for all model configurations. Find and inspect the figure showing the temperature field.

#### To change the model parameters

The model is defined in the file `c:\my_matlab_code\cryogrid\CryoGridExamples\results\test\test.xlsx`.
You may play around with the model parameters and se how the output changes.

For example, you could change the thickness of layer 1:

- Open the excel file
- Find the section `STRAT_layers`
- First column in the defined matrix lists the depth to the bottom of each layer (so row 1 has the depth to the bottom of layer 1)
- The first layer is by default from 0 - 0.5 m (0.5 m thick).
- To change the thickness of the first layer to 1 m, simply change the value from 0.5 to 1.

Rerun the model to see the changes. (Be aware that the out files are overwritten, back them up if you want to store for comparison.)

### Full documentation

There is a pdf file available with more in depth explanation of the model: [CryoGrid_documentation.pdf](./CryoGrid_documentation.pdf)


### <a name="git_download"></a>Getting the code and examples using the git commandline tool

1. Clone the CryoGridExamples repository to a new directory (fx `c:\my_matlab_code\cryogrid`): 

```
git clone --single-branch --branch develop/main https://github.com/CryoGrid/CryoGridExamples.git
```

2. Navigate into the new directory `c:\my_matlab_code\cryogrid\CryoGridExamples`

```
cd c:\my_matlab_code\cryogrid\CryoGridExamples
```

3. Clone the main CryoGrid model code

```
git clone --single-branch --branch develop/main https://github.com/CryoGrid/CryoGrid.git
```

Continue with running the model as described [above](#run_model)



# THE TEXT BELOW HERE IS NOT YET UPDATED

### From the command line

First clone the code with 
```
cd [folder of choice]

# if you have github setup with SSH
git clone git@github.com:CryoGrid/CryoGrid.git

# or if you are not setup with SSH
git clone https://github.com/CryoGrid/CryoGrid.git
```
Github has a [Desktop GUI]() software that can be used to clone the repository.

## Get started as a developer


## Documentation

### Participate to the documentation
The documentation is written in markup language reStructuredText `.rst` or [Markdown](https://www.markdownguide.org/basic-syntax/)(`.md`). rST is the preferred markup language, but both are supported as argued in [this blog post](https://www.ericholscher.com/blog/2016/mar/15/dont-use-markdown-for-technical-docs/).

You can directly edit a page on our [GitHub]() documentation repository, or you can clone the repository, modify it, verify it builds locally, and finally push it back to GitHub.

Creating a new page is as easy as 

### Steps to create documentation from scratch
1. install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)`
2. Then in your terminal (assuming you're on Linux or MacOS)
```shell
# creates a virtual environment
conda create -n docu_env

# Activate the virtual environment
conda activate docu_env

# install pip. Pip is a Python package installer
conda install pip

# install Sphinx which is the engine building the html documentation from the .rst or .md files
pip install sphinx`

# Go to Path
cd /path/to/project
mkdir docs
cd docs

# Build original documentation
sphinx-quickstart
```
From now you can open the documentation build locally with the file `_build/html/index.html`

3. To setup your documentation for markdown, and also use the ReadTheDoc template install
```shell
pip install sphinx-rtd-theme
pip install recommonmark
```

Then in the file `config.py` insert the following
```python
extensions = ['recommonmark', 
				'sphinx_rtd_theme'
]

# and replace
#html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"
```
4. finally you can build the new version of the documentation with
```shell
make html
```
Open the file `_build/html/index.htm` with your browser.

5. Push the project to the github documentation
6. Create a new page
```shell
# create a directory called source that will contain all files
mkdir source
nano source/intro.md
```
Indicate **Sphinx** to seek for this new file by adding in the file `index.rst`:
```rst

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   source/intro
   source/ add here you next content. One file per page. 
```



