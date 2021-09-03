========================================================================
i16sim
========================================================================
Simulation package for the I16 6-circle kappa diffractometer at Diamond Light Source Ltd

Implements `diffcalc <https://github.com/DiamondLightSource/diffcalc>`_ functionality in  `Blender <https://www.blender.org/>`_ and uses it to animate a model of the  `diffractometer <https://www.diamond.ac.uk/Instruments/Magnetic-Materials/I16/layout.html>`_.

If used for your research, please credit Aurys Silinga.

.. contents::
.. section-numbering::

Features
=======================

- **3D model of the diffractometer with mesh error < 1 mm**
- **Collision testing**
- **User interface for visualisation of (pseudo-)motor rotations**
- **Crystal calculations and movements in hkl space via diffcalc commands**
- **Console and script editor for testing experiment scripts**
- **Reading nexus data or GDA state files to show diffractometer state during experiment**
- **Visualisation of reciprocal lattice vectors and azimuthal reference in the laboratory coordinate system**
- **Perspective view from beamline cameras or any point in space**

Installation
=======================
#. Install `Blender <https://www.blender.org/download/>`_ 2.93.1 or newer.
#. Download ``i16sim main.zip`` from GitHub and extract.
#. Open ``diffractometerXX.blend`` in Blender.
#. Bring up the terminal if not already open so sript output can be viewed.
#. Open ``./i16sim main/install/install_i16sim_environment.py`` in Blender's internal script editor.
#. Run the script and wait untill it prints 'Installation finished' in the terminal.
#. Restart Blender.
#. Navigate via the menu to ``Edit > Preferences > Add-ons``.
#. Install the ``./i16sim main/i16sim.zip`` addon and enable it.
#. To access diffcalc commands in the console, run ``from i16sim.commands import *``.

Troubleshooting
----------------------
If steps 5-7 do not work, try installing the ``./i16sim main/install/install_environment.zip`` addon. Enable it via the menu. When it finishes executing and is fully enabled, disable it. Installing this addon forces Blender to create a folder it can write to and and installs all necessary modules to that folder.

If the neither of the above methods worked, you need to manually install the python packages ``numpy``, ``scipy``, and ``h5py`` in a directory that Blender's internal python environment has a path to. Also, copy ``./i16sim main/install/userpref.blend`` to Blender's config folder.


Using the Features
====================
Starting the simulation
--------------------
Enable the ``i16sim`` addon and open up the i16sim tab in the top right corner of the user interface.

Run ``from i16sim.commands import *`` in the console and at the start of every script

Each feature
-------------------
- Collision testing
    Press ``Check for collisions`` button or run ``intersect()`` command.
    If collisions are detected, highlights coliding objects and returns a list of coliding meshes.

- User interface for visualisation of (pseudo-)motor rotations
    Move the sliders in the UI, run ``pos(...)`` commands, or unhide axes of rotation by clicking the 'eye' icon in the outliner.
    The simulation moves in real time.
    
- Crystal calculations and movements in hkl space via diffcalc commands
    See examples folder, `diffcalc <https://github.com/DiamondLightSource/diffcalc>`_ documentation, `diffcalc-core <https://github.com/DiamondLightSource/diffcalc-core>`_ documentation, commands section below, or import nexus data and GDA state files.
    All calculatoins present in diffcalc are possible, but some can require accessing internal diffcalc-core objects.
    
- Console and script editor for testing experiment scripts
    Start each script with command
    ``enable_lm(True)``.
    The scripting environment at the beamline detects syntax errors but cannot predict runtime errors, such as moving the diffractometer to a forbidden position.
    If limits are enabled, all diffcalc commands should throw the same errors as in GDA.
    
- Reading nexus data or GDA state files to show diffractometer state during experiment
    Select the file with the ``file`` widget in the UI and press the ``import form file`` button.
    If the data is present in a file, importing it sets position, diffcalc constraints, crystal lattice, UB matrix, azimuthal reference, and energy. 
    Then it test for collisions.
    
- Visualisation of reciprocal lattice vectors and azimuthal reference in the laboratory coordinate system
    Unhide reciprocal vectors by clicking the 'eye' icon in the outliner. Reciprocal vectors are show if a UB matix is set.

- Perspective view from beamline cameras or any point in space
    Select the camera in the UI 'cameras' tab and click the ``Toggle the camera view`` button.


Commands
====================
Orientation Commands
--------------------

+-----------------------------+---------------------------------------------------+
| **STATE**                                                                       |
+-----------------------------+---------------------------------------------------+
| **-- newub** ({'name'})     | start a new ub calculation, name                  |
+-----------------------------+---------------------------------------------------+
| **-- loadub** ('name'|num)  | load an existing ub calculation                   |
+-----------------------------+---------------------------------------------------+
| **-- lastub** ()            | load the last used ub calculation                 |
+-----------------------------+---------------------------------------------------+
| **-- listub** ()            | list the ub calculations available to load        |
+-----------------------------+---------------------------------------------------+
| **LATTICE**                                                                     |
+-----------------------------+---------------------------------------------------+
| **-- setlat** ()            | interactively enter lattice parameters (Angstroms |
|                             | and Deg)                                          |
+-----------------------------+---------------------------------------------------+
| **-- setlat** ( name, a)    | assumes cubic                                     |
+-----------------------------+---------------------------------------------------+
| **-- setlat** ( name, a, b) | assumes tetragonal                                |
+-----------------------------+---------------------------------------------------+
| **-- setlat** (name, a, b,  | assumes ortho                                     |
| c)                          |                                                   |
+-----------------------------+---------------------------------------------------+
| **-- setlat** (name, a, b,  | assumes mon/hex with gam not equal to 90          |
| c, gamma)                   |                                                   |
+-----------------------------+---------------------------------------------------+
| **-- setlat** (name, a, b,  | arbitrary                                         |
| c, alpha, beta, gamma)      |                                                   |
+-----------------------------+---------------------------------------------------+
| **-- c2th** ([h, k, l])     | calculate two-theta angle for reflection          |
+-----------------------------+---------------------------------------------------+
| **REFERENCE (SURFACE)**                                                         |
+-----------------------------+---------------------------------------------------+
| **-- setnphi** ({[x, y, z]})| sets or displays n_phi reference                  |
+-----------------------------+---------------------------------------------------+
| **-- setnhkl** ({[h, k, l]})| sets or displays n_hkl reference                  |
+-----------------------------+---------------------------------------------------+
| **REFLECTIONS**                                                                 |
+-----------------------------+---------------------------------------------------+
| **-- showref** ()           | shows full reflection list                        |
+-----------------------------+---------------------------------------------------+
| **-- addref**  ()           | add reflection interactively                      |
+-----------------------------+---------------------------------------------------+
| **-- addref** ([h, k, l],   | add reflection with current position and energy   |
| {'tag'})                    |                                                   |
+-----------------------------+---------------------------------------------------+
| **CRYSTAL ORIENTATIONS**                                                        |
+-----------------------------+---------------------------------------------------+
| **-- showorient** ()        | shows full list of crystal orientations           |
+-----------------------------+---------------------------------------------------+
| **-- addorient** ()         | add crystal orientation interactively             |
+-----------------------------+---------------------------------------------------+
| **-- addorient** ([h, k, l],| add crystal orientation in laboratory frame       |
| [x y z], {'tag'})           |                                                   |
+-----------------------------+---------------------------------------------------+
| **UB MATRIX**                                                                   |
+-----------------------------+---------------------------------------------------+
| **-- checkub** ()           | show calculated and entered hkl values for        |
|                             | reflections                                       |
+-----------------------------+---------------------------------------------------+
| **-- calcub**               | (re)calculate u matrix from ref1 and ref2         |
| ( num1|'tag1', num2|'tag2') |                                                   |
+-----------------------------+---------------------------------------------------+
| **-- trialub** ()           | (re)calculate u matrix from ref1 only (check      |
|                             | carefully)                                        |
+-----------------------------+---------------------------------------------------+

Motion Commands
---------------

+-----------------------------+---------------------------------------------------+
| **CONSTRAINTS**                                                                 |
+-----------------------------+---------------------------------------------------+
| **-- con** ()               | list available constraints and values             |
+-----------------------------+---------------------------------------------------+
| **-- con** (<name>, {val})  | constrains and optionally sets one constraint     |
+-----------------------------+---------------------------------------------------+
| **-- con** (<name>,{val},   | clears and then fully constrains                  |
| <name>,{val}, <name>,{val}) |                                                   |
+-----------------------------+---------------------------------------------------+
| **HKL**                                                                         |
+-----------------------------+---------------------------------------------------+
| **-- allhkl** ([h, k, l])   | print all hkl solutions ignoring limits           |
+-----------------------------+---------------------------------------------------+
| **HARDWARE**                                                                    |
+-----------------------------+---------------------------------------------------+
| **-- showlm** ()            | show diffcalc limits and cuts                     |
+-----------------------------+---------------------------------------------------+
| **-- enable_lm** (bool)     | enable or disable all limits                      |
+-----------------------------+---------------------------------------------------+
| **MOTION**                                                                      |
+-----------------------------+---------------------------------------------------+
| **-- sim** (scn, val)       | simulates moving scannable (hkl or sixc)          |
+-----------------------------+---------------------------------------------------+
| **-- sixc** ()              | get Eulerian position                             |
+-----------------------------+---------------------------------------------------+
| **-- pos** (sixc [phi, chi, | move to Eularian position(None holds an axis      |
| eta, mu, delta, gam]        | still)                                            |
+-----------------------------+---------------------------------------------------+
| **-- sim** (sixc, [phi, chi,| simulate move to Eulerian positionsixc            |
| eta, mu, delta, gam])       |                                                   |
+-----------------------------+---------------------------------------------------+
| **-- hkl** ()               | get hkl position                                  |
+-----------------------------+---------------------------------------------------+
| **-- pos** (hkl, [h, k, l]) | move to hkl position                              |
+-----------------------------+---------------------------------------------------+
| **-- pos** ({h  |k | l},    | move h, k or l to val                             |
| val)                        |                                                   |
+-----------------------------+---------------------------------------------------+
| **-- sim** (hkl, [h, k, l]) | simulate move to hkl position                     |
+-----------------------------+---------------------------------------------------+



