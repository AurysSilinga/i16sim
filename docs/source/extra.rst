i16sim developer extras
===================================

How to extend the code
-------------------------

* scannables
	Any new scannables should be defined in the user scannable factory file at ``i16sim.user.additional_scannables.py``.
	See documentation and examples in code on how to use the factory. 
	Any new scannables added to the dictionary will be added to the namespace
 
* functions
	New functions should be defined in ``i16sim.user.additional_functions.py``.
	If defined inside the ``get_additional_functions():`` code block, any new functions will be added to the namespace.
 
* limits
	New limits should be defined in the ``i16sim.util.scannables.create_composite_limits(dc)`` factory.
	See documentation and examples in code on how to use the factory. 
	Limits defined here will be tested every time the diffractometer moves. 
	It is good practice to store your limit parameters in `i16sim.util.params``
 
Known bugs and todos
---------------------

* Directory problems. Some scannables and functions do not load if the simulation is started from a script and its ``os.curdir`` points to a different directory than the ``diffractometerXX.blend`` file.    

* If the coputer falls asleep, the simulation sometimes crashes on Windows.

* The simulation has a sphere error of 1 micron.

* If the UB matrix and the real lattice do not match, diffcalc-core will throw an error that hkl calculations give different results.

* hkl scan should varify each h,k,l value separately and stop incrementing some while others still increase. Currently it stops when any value excedes the maximum.

* Implement con command with dynamic number of parameters. Need to call con(bisect, True), while GDA needs con(bisect).

* Add lattice fixing commands. e.g. get real vectors from UB.

* Add delta offsets. e.g. do do.pil.

* Add pipe physics. Blender has physics simulations for ropes, but the pipes are implemented as simple Bezier curves.



Feedback and support
-----------------------

Please submit an issue to https://github.com/AurysSilinga/i16sim.
