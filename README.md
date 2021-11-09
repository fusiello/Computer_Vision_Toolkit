# Computer Vision Toolkit (CVT)

### Andrea Fusiello, 2018

![alt text](http://www.diegm.uniud.it/fusiello/demo/toolkit/banner.jpg)

This toolkit is a collection of Matlab functions implementing
several Computer Vision algorithms.

**SETUP**: cd to the main directory of CVT, where the file
`cvt_setup.m` is and run: `cvt_setup`. This will add the relevant
folders to the search path.

The main directory containing the relevant function is `m-files`.
Its content will be displayed by `cvt_setup`.

This code runs indifferently on Octave and Matlab. It has been
tested on Octave 4.2.1 and Matlab R2017a. It does not have any
external dependency on toolboxes or packages.

The `test` directory contains scripts that test (mostly) all the
functions in the toolkit.  After setup cd to `test` and
run `testAll`. This will call all the test scripts and report
coverage.

The `m-files/aux_fun` directory contains auxiliary, helper functions. 

The directory `cherubino12` contains 12 JPEG images of a statue,
together with camera matrices (.pm) and masks (.png) (curtesy of
3Dflow srl).  They are used in the testIMG ad testCARVE scripts.

The `License` file applies to all material but the `thirdparty` 
directory, which contains auxiliary functions from other authors, 
where the respective copyright notices apply.

The functions in this toolkit are described in the book:

* Fusiello, A. Visione computazionale. Tecniche di ricostruzione tridimensionale.  Franco Angeli, Milano, 2018, 2a ed.

The book is written in italian, although the code documentation
is in english.


---
Andrea Fusiello                
Dipartimento Politecnico di Ingegneria e Architettura (DPIA)  
Università degli Studi di Udine, Via Delle Scienze, 208 - 33100 Udine  
email: <andrea.fusiello@uniud.it>

---

