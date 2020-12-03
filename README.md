# CryoGrid
This is the modular version of CryoGrid which is based on an object-oriented programming paradigm.

## Get started as a user

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
1. install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
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



