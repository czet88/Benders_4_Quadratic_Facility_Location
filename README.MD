
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
This project contains, to the best of our knowledge, the fastest solver of a well-benchmarked set of non-convex quadratic facility location problems.  

To achieve this, our algorithm uses:
* A fast and efficient matheuristic that solves a "linearized" problem on a "restricted" set of variables.
* A formulation of reduced size obtained by a Benders reformulation of a strong linearization of the quadratic formulation.
* Variable fixing techniques such as Reduced Cost Variable Elimination and Partial Enumeration

To avoid excess overhead incurred when using higher-level languages, the source code of the solution algorithm is coded in C while python is used only in a postprocessing step to create Latex Tables from the output files of the solution algorithm.

## Dependencies
* C/C++ compiler such as GCC for Linux (this usually already comes pre-installed) or use Visual Studio (preferably 2019) for Windows.
* [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) a general mixed integer programming solver. This is a choice made by the authors based on their level of comfort with this software. CPLEX is called throughout the program via its Callable Library API.

<!-- GETTING STARTED -->
## Quick Start

To get a local copy of the algorithm you can open your preferred terminal and, assuming you have previously installed [GIT](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git), you can run the following command:
```
git clone https://github.com/czet88/Benders_4_Quadratic_Facility_Location.git
```
This will copy all the contents of the repository into a folder named *Benders_4_Quadratic_Facility_Location*.



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
 




<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[license-shield]: https://img.shields.io/github/license/czet88/Benders_4_Quadratic_Facility_Location.svg?style=for-the-badge
[license-url]: https://github.com/czet88/Benders_4_Quadratic_Facility_Location/blob/master/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/cazetina/




