Grid Specification 
==============================================

Rayleigh solves the fluid equations in spherical-shell geometry.  As the poles are included, the grid is fully specified by providing four pieces of information:

- The coordinates of the computational domain's radial boundaries of the domain, :math:`r_\mathrm{min}` and :math:`r_\mathrm{max}`
- The number of radial grid points, :math:`N_r`
- The number of latitudinal grid points, :math:`N_\theta`

The number of longitudinal grid points, :math:`N_\phi` , is always twice :math:`N_\theta`.
The total number of gridpoints for a Rayleigh simulation is then given by :math:`2N_rN_\theta^2`.
Note that both :math:`N_r` and :math:`N_\theta` must be even.   Rayleigh's computational grid is specified using the problemsize namelist in the main_input file.
A quick reference for all problemsize-namelist variables is provided in the :ref:`namelist documentation <namelists>`.  In this section, we discuss in detail how to define
Rayleigh's grid using these variables.



Standard grid specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We begin by discussing how to define a grid employing a single Chebyshev domain in radius, meaning that a single Chebyshev expansion is carried out over the domain
:math:`r_\mathrm{min} \le r \le r_\mathrm{max}`.  This is probably the most common grid setup employed in Rayleigh.

The problemsize variables *n_r* and *n_theta* provide values
for :math:`N_r` and :math:`N_\theta` respectively.  Similarly, *rmin* and *rmax* define the value for :math:`r_\mathrm{min}` and `r_\mathrm{max}`.
If we wanted to define a spherical shell extending from r=1.0 to r=2.0, with :math:`N_r=48` and :math:`N_\theta=96`,
out problemsize namelist should look like:

::

   &problemsize_namelist
    n_r = 48
    n_theta = 96
    rmin = 1.0
    rmax = 2.0    
   /

Note that :math:`N_r` and :math:`N_\theta` may also be specified at the command
line using the flags -nr and -ntheta, e.g.:

::

   mpiexec -np 8 ./rayleigh.opt -nr 48 -ntheta 96

Doing so will override any values supplied via main_input.  This can be particularly useful when scripting performance analyses on a new machine.

If desired, a user may instead specify the radial domain bounds in terms of the shell aspect ratio :math:`\chi=r_\mathrm{min}/r_\mathrm{max}`, and
the shell depth :math:`r_\mathrm{max}-r_\mathrm{min}`.  This is accomplished using the using the aspect_ratio and shell_depth problemsize variables.  The example below describes a
grid equivalent to the one described above.

::

   &problemsize_namelist
    n_r = 48
    n_theta = 96
    aspect_ratio = 0.5
    shell_depth = 1.0
   /



Rayleigh's horizontal resolution (:math:`N_\theta\times N_\phi`) may alternatively be described in terms of spherical harmonics.  The maximum Legendre degree employed in Rayleigh's truncated spherical harmonic expansion is denoted by :math:`\ell_\mathrm{max}`, and the total number of degrees by :math:`N_\ell`.  
These two variables are related to :math:`N_\theta` via

.. math:: N_\ell = \ell_\mathrm{max}+1 = \frac{2}{3}N_\theta,

and they are described by the problemsize variables n_l and l_max.  Thus, the examples

::

   &problemsize_namelist
    n_r = 48
    l_max = 63
    rmin = 1.0
    rmax = 2.0       
   /

and

::

   &problemsize_namelist
    n_r = 48
    n_l = 64
    aspect_ratio = 0.5
    shell_depth = 1.0
   /

both describe a grid extending from r=1.0 to r=2.0, with :math:`N_r=48` and :math:`N_\theta=96`.

Defining multiple Chebyshev domains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some instances, it may be advantageous to describe the radial grid using
multiple Chebyshev domains.  The most common use case probably occurs
when the system under consideration is characterized by layers
subject to different physical conditions.  For instance, models that
include regions that are both superadiabatically and subadiabatically stratified
might employ a different Chebyshev expansion within each domain.  Similarly so
for geodynamo models that include the solid inner core.

When describing a grid with *N* Chebyshev domains, the main_input file must first
supply *N+1* points :math:`r_i` that define the bounds of these domains.
The *ith* Chebyshev domain will span the interval :math:`r_i \le r \le r_{i+1}`, and the global domain bounds
are defined such that

.. math:: r_0 \equiv r_\mathrm{min}\,\,\,\,\mathrm{and}\,\,\,\,r_{N+1}\equiv r_\mathrm{max}.


I was here.  The next example is good.  Note that we use ncheby.  Note that a radial point will be repeated.


It is possible to run Rayleigh with multiple, stacked domains in the
radial direction. Each of these is discretized using their own set of
Chebyshev polynomials. The boundaries and number of polynomials can be
set for each domain indiviadually, which makes it possible to control
the radial resolution at different radii.

To use this feature the problem size has to be specified using
``domain_bounds`` and ``ncheby`` instead of ``rmin``, ``rmax``, and
``n_r``. ``ncheby`` takes a comma-separated list of the number of radial
points to use in each domain. ``domain_bounds`` takes a comma-separated
list of the radii of the domain boundaries, starting with the smallest
radius. It has one element more than the number of domains. This is an
example of two radial domains, one covering the radii 1 to 2 with 16
radial points, the other the radii 2 to 4 with 64 radial points.

::

   &problemsize_namelist
    domain_bounds = 1.0, 2.0, 4.0
    ncheby = 16, 64
   /

Radial values in the diagnostic output will be repeated at the inner
domain boundaries. Most quantities are forced to be continuous at these
points.

Controlling radial dealiasing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
