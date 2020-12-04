# Inheritance diagrams

This folder contains inheritance diagrams for selected CryoGrid classes. The diagrams are created in Unified Markup Language (UML), and exported to the Scalable Vector Graphics (SVG) format, which can be opened/viewed in all modern browsers.

The folder contains both the *.uml and the *.svg versions.
The diagrams contains links to the classes, methods and properties in the actual code.
If the SVG diagrams are opened in the MatLab web browser, clicking one of these links will direct the MatLab editor to the corresponding definition in the code. This functionality only works in the MatLab browser, however, the diagrams are viewable in any browser.

## Getting started

For the links in the diagrams to function, the *.svg files must contain absolute paths to the CryoGrid class files.
Use the function `localize_svgs` in the 'UMLs' folder to convert all paths in all *.svg files to the correct path for your particular installation folder.
The code currently assumes that the 'UMLs' folder is at the same level as the 'modules' folder:

```
- My_CryoGrid_folder
   - ...
   - modules
   - UMLs
   - ...
```

Navigate to the 'UMLs' folder in your MatLab console (or use the path selector):

```
cd c:\path\to\cryogrid_install\UMLs
```

The run the `localize_svgs` script:

```
localize_svgs
```

To open the Matlab browser and view the index of diagrams available, issue the following command in the Matlab console:

```
open_index
```

From the index, you can open inheritance diagrams of the classes you are interested in, and use the links in the diagrams to jump to the corresponding code in the MatLab editor.

## Commiting changes to the GIT repository

When committing changes to the GIT repository, the links in the svg files should be relative (`./modules/*`), so that you do not commit links with the absolute paths to your particular install directory.

For this purpose, a `relativize_svgs` script is supplied, which will convert absolute paths to relative paths in the svg files.
Run this script before committing, or simply avoid to commit changes to the UML diagrams (unless you specifically changed them and know what you are doing).

## Regenerating inheritance diagrams

For generation (regenerating) of the inheritance diagrams, use the script `create_all_uml`. This script will generate all the uml diagrams (which ones are defined in separate scripts), and regenerate the index file.

The code relies on the toolbox `m2uml` and a `plantuml.jar` file.
(Explanation of how to install and run this should be added...)



