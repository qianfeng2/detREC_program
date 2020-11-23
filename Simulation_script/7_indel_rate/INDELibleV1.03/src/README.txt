UNIX
----

To compile the source code on UNIX use a suitable compiler such as g++.

  e.g.  g++ -o indelible -O4 indelible.cpp

but insert the relevant compiler flags for your system.

I use:
  
   g++ -o indelible -m64 -O4 -march=opteron indelible.cpp -lm

but there are many other flags that can be added.


MAC - (N.B. pre-compiled executables in "bin" folder)
-----------------------------------------------------

The process is similar to that of UNIX.  

Just open a terminal window, use 'cd' to navigate to the folder where the 
source code is, then type the compile command.

  e.g. g++ -o indelible -O4 -march=pentium-m indelible.cpp -lm

for new intel-based MACs.

The following two seem to work on older G4 and G5 systems:

  e.g. g++ -o indelible -O4 -mcpu=G4 indelible.cpp -lm
  e.g. g++ -o indelible -O4 -mcpu=G5 indelible.cpp -lm


WINDOWS - (N.B. pre-compiled executables in "bin" folder)
---------------------------------------------------------

Using Microsoft Visual Studio 6.0 simply load the file indelible.cpp making 
sure it is in the same folder as the other source files.
* Click on the compile button.  
* Once compiled go to the build menu and choose "set active configuration".
* Then choose "RELEASE" otherwise you will get a very slow version with debug code in.
* Then re-compile and build.
* The INDELible executable will be in the "release" folder 

Please note if you do not want to use the executable and wish
to compile INDELible yourself on Windows systems then please
include the line

 #define WINDOWS

at the top of the file indelible.cpp