The physunits toolbox is an attempt to confer dimensional
"awareness" to the MATLAB environment. The motivation for this, as
well as a suggested way to go about it, are explained in
% Petty, G.W., 2001: Automated computation and consistency checking of
% physical dimensions and units in scientific programs.  Software -
% Practice and Experience (Volume 31, Issue 11, 19 June 2001)

The author of the same has also made available for download a
FORTRAN90 module that implements this ides in the FORTRAN language.
The module and paper can be downloaded from
http://rain.aos.wisc.edu/~gpetty/physunits.html

The physunits toolbox is based on the above module, and expands
it, while trying to adhere consistently to MATLAB standards and
practices.

The physunits toolbox works by defining the PREAL class, with a PREAL
data type and overloaded functions to support it. An object of type
PREAL represents a physical quantity. It has two fields:
'value' - The numerical value, which must be a numerical type with
zero imaginary part.
'units' - A vector of 7 numbers, representing the physical dimensions
associated with 'value'. The format for the 'units' vector is [length,
weight, time, temperature, electric current, amount of matter,
illumination, luminous intensity],
but this format is usually transparent to the user, who defines her
variables via an interface structure.

The interface structure contains predefined variables of type preal,
representing the most useful SI units as well as many other units,
derived units, constants of nature, parameters, etc. Get this
structure by calling the setUnits function.

The physunits distribution contains the following files:

@preal (directory) -
In MATLAB, classes are defined by a directory of the class name
preceded by a @ sign. The @preal directory contains the preal class
definition and methods. To use the preal class, the directory
IMMEDIATELY ABOVE @preal should be added to the MATLAB path.

setUnits.m -
MATLAB .m file with predefined units and constants. See documentation
for usage.

info.xml - Pointer to the html documentation file. Used by the MATLAB
help browser.

html (directory) -
Containing the html help and documentation files.

readme.txt - this file.

physunits.m - a function for switching the physunits on or off.

iceball.m - a short sample program that demonstrates the preal class.