# i16sim
Simulation package for the i16 6-axis kappa diffractometer at Diamond Light Source Ltd

Diffcalc - A Diffraction Condition Calculator for Diffractometer Control
========================================================================
bv
Diffcalc is a python/jython based diffraction condition calculator used for
controlling diffractometers within reciprocal lattice space. It performs the
same task as the fourc, sixc, twoc, kappa, psic and surf macros from SPEC.

There is a `user guide <https://diffcalc.readthedocs.io/en/latest/youmanual.html>`_ and `developer guide <https://diffcalc.readthedocs.io/en/latest/developer/contents.html>`_, both at `diffcalc.readthedocs.io <https://diffcalc.readthedocs.io>`_

|Travis| |Read the docs|

.. |Travis| image:: https://travis-ci.org/DiamondLightSource/diffcalc.svg?branch=master
    :target: https://travis-ci.org/DiamondLightSource/diffcalc
    :alt: Build Status

.. |Read the docs| image:: https://readthedocs.org/projects/diffcalc/badge/?version=latest
    :target: http://diffcalc.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. contents::

.. section-numbering::

Software compatibility
----------------------

- Written in Python using numpy
- Works in Jython using Jama
- Runs directly in `OpenGDA<http://www.opengda.org>`
- Runs in in Python or IPython using minimal OpenGda emulation (included)
- Contact us for help running in your environment

Diffractometer compatibility
----------------------------

Diffcalc’s standard calculation engine is an implementation of [You1999]_ and
[Busing1967]_. Diffcalc works with any diffractometer which is a subset of:

 .. image:: https://raw.githubusercontent.com/DiamondLightSource/diffcalc/master/doc/source/youmanual_images/4s_2d_diffractometer.png
     :alt: 4s + 2d six-circle diffractometer, from H.You (1999)
     :width: 50%
     :align: center

Diffcalc can be configured to work with any diffractometer geometry which is a
subset of this. For example, a five-circle diffractometer might be missing the
nu circle above.

Note that the first versions of Diffcalc were based on [Vlieg1993]_ and
[Vlieg1998]_ and a ‘Vlieg’ engine is still available.  There is also an engine
based on [Willmott2011]_. The ‘You’ engine is more generic and the plan is to
remove the old ‘Vlieg’ engine once beamlines have been migrated.

Installation
------------

Check it out::

   $ git clone https://github.com/DiamondLightSource/diffcalc.git
   Cloning into 'diffcalc'...

At Diamond Diffcalc may be installed within an OpenGDA deployment and is
available via the 'module' system from bash.
