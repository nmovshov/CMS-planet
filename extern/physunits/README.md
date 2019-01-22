# PHYSUNITS - dimensioned MATLAB variables with enforced dimensional consistency

The `physunits` toolbox is an attempt to confer "dimensional awareness" to the
MATLAB environment. The motivation and a suggested algorithm are explained in:
Petty, G.W., 2001. Automated computation and consistency checking of physical
dimensions and units in scientific programs. _Software: Practice and Experience_,
31(11), pp.1067-1076. The author also made available for download a FORTRAN90
module that implements this idea in the FORTRAN language. The module and paper can
be downloaded from http://rain.aos.wisc.edu/~gpetty/physunits.html. The
`physunits` toolbox is based on that module and expands it while trying to adhere
consistently to MATLAB standards and practices (_from 2006_).

The `physunits` toolbox works by defining the `preal` data type and overloaded
functions to support it. A variable of type `preal` represents a physical
quantity. It has two _private_ fields: 'value' - the numerical value, which must
be a numeric type with zero imaginary part, and 'units' - a vector of 7 numbers,
representing the physical dimensions associated with 'value'. The format for the
'units' vector is [length, mass, time, temperature, electric current, amount of
matter, illumination, luminous intensity]. But these details are never needed by
the user, who defines his or her variables via an interface structure. The
interface structure contains predefined variables of type `preal`, representing
the most useful SI units as well as many other units, derived units, constants of
nature, parameters, etc. Get this structure by calling the `setUnits` function.
Create your own dimensioned variables by multiplying a number with an existing
`preal` variable.

See the HTML documentation for more details and examples.

### Setup

The physunits toolbox is implemented as an old-style (pre-2008) MATLAB class. In
this style MATLAB classes are defined by a directory of the class name preceded by
a `@` sign. The `@preal` directory contains the `preal` class definition and
methods. **To use the `preal` class, the directory immediately above `@preal` must
be on the MATLAB search path.**
