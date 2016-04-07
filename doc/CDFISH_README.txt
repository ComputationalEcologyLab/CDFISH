======
README
======

-------------------  
CDFISH 0.55 release
-------------------
       
Welcome to the CDFISH v0.55 release!  CDFISH is a program to simulate gene flow in complex stream landscapes for a wide range of biological scenarios of aquatic species.  This release includes installation instructions. For more information, please see the CDFISH user manual.
       
Program Contributors: Erin Landguth
Link: http://cel.dbs.umt.edu/software/CDFISH/
Version: 0.55 
Python: 2.7.2
Release Date: 2013.10.09
README Update: 2013.10.09 (ell)
       
--------
Contents
--------
       
Included in this release are the following:       
CDFISH_README.txt - this file
CDFISH_user_manual.pdf - more complete guide

Example Test Files:
cdfish_test.cd - sample cost distance matrix for 4 subpopulations
cdfish_test.xy - sample n-(x,y) for 64 individuals (16 per subpopulation)
cdfish_test.csv - input run file with parameters corresponding to the sample files

Source Code:				
CDFISH.py - Python driver code and run file
CDFISH_Disperse.py - Python library for dispersal functions
CDFISH_GetMetrics.py - Python library for metric functions
CDFISH_Mate.py - Python mating library
CDFISH_Modules.py - Python main library
CDFISH_Offspring.py - Python library for offspring functions
CDFISH_PostProcess.py - Python post-processing library
CDFISH_PreProcess.py - Python pre-processing library
       
---------------------------------------
Requirements and Pre-requisite Software
---------------------------------------

1. Baseline Requirements.  CDFISH requires the Python2.7.x interpreter and NumPy package.  Remember that Python modules usually require particular Python interpreters, so be sure the version ID for any external Python module or package (e.g. NumPy or others) matches the version of your Python interpreter (normally v2.7.x). 

-----------------------
CDFISH Installation
----------------------- 

1. Download CDFISH.zip (for LINUX or WINDOWS OS) or cdfish_setup.exe (for windows OS only).

2. Install CDFISH.  Next, install the CDFISH software itself by unpacking the zip archive supplied (or clicking on the cdfish_setup.exe file and following the installation instructions). At this point you should be able to execute the supplied test inputs.

------------------
Example CDFISH Run
------------------

The example run is for 64-points representing 16 individuals in 4 subpopulations (cdfish_test.xy) with an example cost distance matrix (cdfish_test.cd). To run the following example, follow these steps:

1. Double check that the files provided in the archive are in the same directory.  

2. The included file cdfish_test.csv specifies the parameters that can be changed and used in a sample CDFISH run.  Open cdfish_test.csv in your editor of choice. 
 
3. There will be 3 lines of information in cdfish_test.csv: a header line and 2 lines of information corresponding to 2 separate CDFISH runs (batch process).  See the CDFISH_user_manual.pdf that contains a breakdown for each column header and the parameters that can be changed.    

4. Save cdfish_test.csv in the same format, a comma delimited file.  

5. Start the program with a command line: For example, if you use python from the command line, then open a terminal window and change your shell directory to the CDFISH home directory.  

6. Run the program: There are a number of ways to run this program. If you are using a command shell you can run the program by typing “python CDFISH.py cdfish_test.csv”.

7. Check for successful model run completion: The program will provide step-by-step output in the Shell window.  Once completed, a simulation time will be printed out and folders batchrun0mcrun0, batchrun1mcrun0, and batchrun1mcrun1 will be created in your CDFISH home directory to store output from the separate batch and/or Monte-Carlo runs. Each of these folders will have a unique date/time stamp preceding ‘batchrun0mcrun0’ in case you want to run multiple CDFISH runs in this same directory.  The program will also provide a log file with program steps in your CDFISH home directory (cdfish.log).

Happy Simulations!

Erin.

Contact Information
Erin Landguth
Computational Ecology Laboratory
The University of Montana
32 Campus Drive
Missoula MT, 59812-1002
(406) 243-2393
Erin.landguth@mso.umt.edu
