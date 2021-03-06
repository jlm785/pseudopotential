# pseudopotential
Modernized version of a pseudopotential generation code

This is a fortran90 version of the venerable atomic configuration and pseudopotential generation code, with extensive documentation and some new capabilities.

The initial code was written by Sverre Froyen at the University of California at Bekeley in the early 1980s.  Norm Troullier modified it to generate the Troullier-Martins pseudopotental at the University of Minnesota at Minneapolis.  It has been maintained by José Luís Martins and its collaborators.  It has been on the web since early 1998 (https://web.archive.org/web/19980111014249/http://bohr.inesc.pt/~jlm/pseudo.html)

The 6.x.y version is a conversion from the old fortran77 to modern fortran90 with minimal changes to the core algorithms.  

Early references of the code are:
Nonlinear ionic pseudopotentials in spin-density-functional calculations, Steven G. Louie, Sverre Froyen, and Marvin L. Cohen, Phys. Rev. B 26, 1738 (1982)
Efficient pseudopotentials for plane-wave calculations, N. Troullier and José Luís Martins, Phys. Rev. B 43, 1993 (1991).
