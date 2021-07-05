
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![LinkedIn][linkedin-shield]][linkedin-url]
[![License][license-shield]][license-url]




<!-- PROJECT LOGO -->
<br />
<p align="center">
  
  <h3 align="center">An Exact Algorithm for Large-Scale Non-Convex Facility Location</h3>

  <p align="center">
    Source code for the paper authored by Ivan Contreras, Carlos A. Zetina, and Sachin Jayaswal 
    <br />
    <a href="https://sites.google.com/view/carloszetina/Research/Publications"><strong>Paper available here</strong></a>
  </p>
</p>




<!-- ABOUT THE PROJECT -->
## About The Project
This project contains, to the best of our knowledge, the fastest solver of a well-benchmarked set of non-convex quadratic facility location problems. The repository is structured for fast setup both in Windows (using Visual Studio 2019) and Linux (using bash scripts).   

To achieve this, our algorithm uses:
* A fast and efficient matheuristic that solves a "linearized" problem on a "restricted" set of variables.
* A formulation of reduced size obtained by a Benders reformulation of a strong linearization of the quadratic formulation.
* Variable fixing techniques such as Reduced Cost Variable Elimination and Partial Enumeration

To avoid excess overhead incurred when using higher-level languages, the source code of the solution algorithm is coded in C while python is used only in a postprocessing step to create Latex Tables from the output files of the solution algorithm.

## Dependencies
* C/C++ compiler such as GCC for Linux (this usually already comes pre-installed) or use Visual Studio (preferably 2019) for Windows.
* [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) [version 12.9]- a general mixed integer programming solver. This is a choice made by the authors based on their level of comfort with this software. CPLEX is called throughout the program via its Callable Library API.

<!-- Quick Start -->
## Quick Start

To get a local copy of the algorithm you can open your preferred terminal and, assuming you have previously installed [GIT](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git), you can run the following command:
```
git clone https://github.com/czet88/Benders_4_Quadratic_Facility_Location.git
```
This will copy all the contents of the repository into a folder named *Benders_4_Quadratic_Facility_Location*.

### Windows
Open the visual studio project named CHLPA.sln located in the root folder. 

The project is currently setup to read the **input file** (located in the "Inputs" directory) named "in_u_f.txt". This file contains the execution parameters to solve the Uncapacitated Single Allocation Hub Location Problem with Fixed Costs of the well-known AP instances (located in the "Data" directory). Upon execution, the code will write the **output file** "R_in_u_f.txt" into the relative path ```\Results\Raw_Results``` which can then be processed by the python script located in the Results folder.

Once the project has opened you can run it in release mode by pressing ```Ctrl```+```F5```.

The project is setup assuming:
* The CPLEX include directory is located at C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\include\ilcplex
* The CPLEX library file used in the Linker is located C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\lib\x64_windows_vs2017\stat_mda\cplex1290.lib

### Linux
First, modify the Library path in the Makefile found in the src\ folder by changing the
```
LIB = /encs/pkg/cplex-12.9.0/root/cplex
```
to
```
LIB = ${Your_Path_To_CPLEX_root_Folder}
```

Next, in the root folder, execute the following command

```
./run_tests.sh
```

This will execute the code reading the **input file**  (located in the "Inputs" directory) named "in_u_f.txt". This file contains the execution parameters to solve the Uncapacitated Single Allocation Hub Location Problem with Fixed Costs of the well-known AP instances (located in the "Data" directory). Upon execution, the code will write the **output file** "R_in_u_f.txt" into the relative path ```\Results\Raw_Results``` which can then be processed by the python script located in the Results folder.

<!-- Code details -->
## Code details
The following section contains more details about the compilation process and usage of the code. It is in particular for those who wish to run different instances and/or use this codebase in their own research.

<!-- FOLDER STRUCURE -->
### Folder Structure

To keep the root folder as clean as possible, the authors have used the following structure to organize files.

* Data - contains the files with the data of the instances to be executed (AP instances and I instances).
* Inputs - contains the files that are read as inputs.
* Results - Folder containing both raw and processed output of the code.  
  * Latex - contains the result of the python script that analyzes the raw output of the code.
  * Raw_Results -contains the raw files produced by the code.
  * Reported_Results - contains the files with the results reported in the paper.
* Scripts - contains shell scripts used in the LINUX setup
* src - contains the C source code of the algorithm used to solve Quadratic Facility Location Problems

At the root folder the following visual studio files for quick _Windows setup_. (Note: For the smoothest process we suggest Visual Studio 2019 community edition) <br />
+-- CHLPA.sln <br />
+-- CHLPA.vcxproj <br />
+-- CHLPA.vcxproj.user <br />

The licence file (GNU GLP3.0) and README are also located in the root folder <br />
+-- LICENCE <br />
+-- README.MD <br />

Finally, a shell script is also located at the root folder for quick _LINUX setup_. (Note: The path to CPLEX must be adjusted in the Makefile) <br />
+-- run_tests.sh <br />

The code is setup to exploit the before mentioned folder structure keeping the instance data, execution parameters, and results in their corresponding folders. This is done for the command line arguments to only be file names instead of file paths.

<!-- Compilation -->
### Compilation
The code is setup for execution from the root folder. Hence, if trying to compile in your own terminal directly, we recommend moving the executable file to the root folder.

There are two ways in which you can compile the code in Linux.
* From the root folder run the command ```./Scripts/update.sh```.
* Enter the src folder and run the command ```make```. Then cut and paste it into the root folder.

In windows, as mentioned before, one can simply use the visual studio project that comes in this repo.

<!-- Usage -->
### Usage

``` 
./B4QFL input_file_name output_file_name [Heuristic_Level] 
```
By default, the code will use the maximum Heuristic_Level (2) to find a feasible solution at the beginning of the solving process to be able to perform variable elimination

Supported options are:
* 0- Do not use any heuristic method
* 1- Use the matheuristic based on ignoring the quadratic costs
* 2- Use the matheuristic along with a local search algorithm




<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[license-shield]: https://img.shields.io/github/license/czet88/Benders_4_Quadratic_Facility_Location.svg?style=for-the-badge
[license-url]: https://github.com/czet88/Benders_4_Quadratic_Facility_Location/blob/master/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/cazetina/




