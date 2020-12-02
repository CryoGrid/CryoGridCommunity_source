# CryoGrid
This is the modular version of CryoGrid which is based on an object-oriented programming paradigm.

## Get started

First clone the code with 
```
cd [folder of choice]

# if you have github setup with SSH
git clone git@github.com:CryoGrid/CryoGrid.git

# or if you are not setup with SSH
git clone https://github.com/CryoGrid/CryoGrid.git
```
Github has a [Desktop GUI]() software that can be used to clone the repository.


## Documentation

1. install Miniconda
2. `conda create -n docu_env`
3. `conda activate docu_env`
4. `conda install pip`
5. `pip install sphinx`

```
cd /path/to/project
mkdir docs
cd docs
sphinx-quickstart
