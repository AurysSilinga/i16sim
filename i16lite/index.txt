i16sim command reference
************************

Created on 19/07/2021

@author: Aurys Silinga

Diffcalc emulator for the Python console in Blender.

animate commands should have options to visualise an eta,chi,phi axes

add lattice fixing commands. get real vectors from UB delta offsets

#documentation mention what happens if lattice is wrong

bug testing scripts

con setting with 5 constraint parameters

intesect and animate parameters for scans

turn off python tooltips

sphere error of 1 micron

check user preferences tab for colour of grids

falls asleep crash

edit preferences so trackball is enabled and drag to move is disabled

scan hkl check if any are ourside limit. Should varify each separately
and stop at maximum

class i16sim.diffcalc_emulator.DiffcalcEmulator(position=None, ubcalc=None, cons=None, hklcalc=None, wl=1, k_angles=[0, 0, 0], ubname='test', limits=False, collision_test=False, verbose=False)

   Moves the simulation in Blender and provides commands that mimic
   diffcalc in GDA

   addorient(hkl, xyz, name='', position=None)

      Add a reference orientation.

      Adds a reference orientation in the diffractometer coordinate
      system.

      Example:

         #if at current position c* is in the same direction as +z in the sample coordinate system
         addorient([0,0,1], [0,0,1], '+z||c*') ```

      Parameters:
         * **hkl** (*Tuple**[**float**, **float**, **float**]*) – hkl
           index of the reference orientation

         * **xyz** (*Tuple**[**float**, **float**, **float**]*) – xyz
           coordinate of the reference orientation in the diffcalc
           sample coordinate frame alternatively, a [x,y,z] vector
           parallel to the (h * a* + k * b* + l * c*) vector

         * **name** (*str*) – identifying tag for the reflection

         * **position** (*Position**, **optional**, **default is
           current position*) – position of reflection as diffcalc
           position object

   addref(hkl, name='', position=None)

      Add a reference reflection.

      Adds a reflection position in degrees and in the systems
      internal representation.

      -[ Example ]-

      addref([1,1,1], “111_first_try”) #adds a reflection

      Parameters:
         * **hkl** (*Tuple**[**float**, **float**, **float**]*) – hkl
           index of the reflection

         * **name** (*str**, **optional**, **default = None*) –
           identifying tag for the reflection

         * **position** (*Position**, **optional**, **default is
           current position*) – position of reflection as diffcalc
           position object

   allhkl(hkl, *args)

      prints all potential positions for a certain hkl regardless of
      limits.

      -[ Example ]-

      allhkl([1,1,1]) #prints all theoretically possible eulerian
      positions.

      Parameters:
         **hkl** (*[**float**,**float**,**float**]*) – Miller indices.

   c2th(hkl)

      Calculate two-theta scattering angle for a reflection

      Parameters:
         **hkl** (*[**float**,**float**,**float**]*) – Miller indices.

      Returns:
         **tth** – 2*theta.

      Return type:
         float

   calcub(*args)

      Calculate UB matrix. Enables hkl calculations.

      Calculate UB matrix using two reference reflections and/or
      reference orientations.

      By default use the first two reference reflections when
      provided. If one or both reflections are not available use one
      or two reference orientations to complement mission reflection
      data.

      -[ Example ]-

      calcub(1,2) #uses reflections 1 and 2 to calculate

      Parameters:
         * **idx1** (*int** or **str**, **optional*) – The index or
           the tag of the first reflection or orientation.

         * **idx2** (*int** or **str**, **optional*) – The index or
           the tag of the second reflection or orientation.

   checkub(*args)

      Print the hkl of added reflections next to the calculated hkl
      for corresponding positions.

      -[ Example ]-

         checkub():
         # id, hkl,    hkl_calculated
         # 1, (1,1,1), (1,1,1)
         # 2, (0,0,1), (0,0,0.98) #means the lattice parameters are not exactly right

   clear(keep_scannables=True)

      Clear previous calculations

      -[ Example ]-

      clear()

      Parameters:
         **keep_scannables** (*bool**, **optional*) – Copies
         scannables over. if keep==false, need to reload namespace
         afterwards. The default is True.

      Returns:
      Return type:
         None.

   con(*args)

      diffcalc con function. If no arguments given, print all
      constraints in a table elif the the only argument is a
      scannable, print its constraint. elif there are two arguments, a
      scannable and a value, set the scannable constraint to that
      value. elif there are multiple scannable and value pairs, set
      each scannable’s constraint to the corresponding value.

      -[ Example ]-

      con() # show constraint table con(psi) # show value of psi
      constraint con(bisect, True) # enable bisect constraint

      Parameters:
         ***args** (*[**scannable1**, **value1**, **scannable2**,
         **value2**, **scannable3**, **value3**]*) –

      Returns:
         value of constraint if con is called with one parameter.

      Return type:
         float or bool

   enable_collision_test(if_use=True)

      Enable or disable checking for collisions after every movement.
      Slow on most computers.

      -[ Example ]-

      enable_collision_test(True) #will test for collisions all the
      time

      Parameters:
         **if_use** (*boolean**, **optional*) – if enable. The default
         is True.

   enable_lm(if_use=True)

      Enable or disable safety limits

      -[ Example ]-

      enable_lm(False) #safety limits are disabled for the simulation

      Parameters:
         **if_use** (*bollean**, **optional*) – if enable. The default
         is True.

   enable_verbose(if_use=True)

      Enable or disable printing the full position string every time
      the diffractometer moves

      -[ Example ]-

      enable_verbose(True) #will print every movement

      Parameters:
         **if_use** (*boolean**, **optional*) – if enable. The default
         is True.

   fdict_print(d)

      Print a dictionary of floats nicely

      Parameters:
         **d** (*dict {str:float}*) –

   fk(**args)

      disable inverse kinematics

   get_base_namespace(other_diffcalc_emulator=None)

      returns a dictionary with the base namespace.

      -[ Example ]-

      dc=DiffcalcEmulator() globals().update(dc.get_namespace()) #to
      enable shorthand sytax and scannables

      Returns:
         **dict**

      Return type:
         {str:Scannable, str:method}

   get_constraints()

      get current diffcalc-core constraints object

      Returns:
      Return type:
         Position

   get_dist(pos1, pos2)

      Get estimate of how far away two position are in eulerian space

      Parameters:
         * **pos1** (*Position*) – A position object describing
           eulerian angles.

         * **pos2** (*Position*) – A position object describing
           eulerian angles.

      Returns:
         **dist** – Geometric sum of all eulerian rotations if going
         from pos1 to pos2.

      Return type:
         float

   get_hklcalc()

      get current diffcalc-core hkl calculation object

      Returns:
      Return type:
         HklCalculation

   get_namespace()

      returns a dictionary with the full namespace (base and user
      defined).

      -[ Example ]-

      dc=DiffcalcEmulator() globals().update(dc.get_namespace()) #to
      enable shorthand sytax and scannables

      Returns:
         **dict**

      Return type:
         {str:Scannable, str:method}

   get_position(asdict=False)

      get current position eulerian angles as a tuple or a dictionary
      in i16 sixc ordering

      -[ Example ]-

      get_postion(asdict=True) #returns a dictionary of eulerian
      angles

      Parameters:
         **asdict** (*bool**, **optional*) – if return as a
         dictionary. The default is False.

      Returns:
         **[phi** – i16 sixc order of eulerian angles

      Return type:
         float, chi:float, eta:float, mu:float, delta:float,
         gam:float] or dict {str:float}

   get_position_ob()

      Get current position as diffcalc object

      Returns:
         current diffcalc position object.

      Return type:
         Position

   get_scannables()

      Get dicionary of all scannables

      Returns:
         **dict**

      Return type:
         {str:Scannable}

   get_sixc()

      get eulerian angles in the same order as in GDA at i16

      -[ Example ]-

      get_sixc()

      Returns:
         **[phi** – eulerian angles in intuitive i16 order.

      Return type:
         float, chi:float, eta:float, mu:float, delta:float,
         gamma:float]

   get_ubcalc()

      get current diffcalc-core ub calculation object

      Returns:
      Return type:
         UBCalculation

   ik(**args)

      enable inverse kinematics

   inc(key, val)

      increments a scannable’s value

      Parameters:
         * **key** (*Scannable*) – the scannable you want to
           increment.

         * **val** (*float** or **list of floats*) – The value by
           which to increment

      Returns:
      Return type:
         None.

   init_position_ob(e_angles)

      Create a position object

      Parameters:
         **e_angles** (*[**mu:float**, **delta:float**,
         **gamma:float**, **eta:float**, **chi:float**,
         **phi:float**]*) – position tuple.

      Returns:
         **pos** – position object.

      Return type:
         Position

   inlimits(position=None, raise_error=False)

      if given position is within safety limits. If no parameters
      given, checks current position.

      -[ Example ]-

      inlimits(init_position_ob([0,0,0,0,90,0])) #check if chi 90 is
      allowed by the limits

      Parameters:
         * **position** (*Position**, **optional*) – A position object
           describing eulerian angles. The default is current
           position.

         * **raise_error** (*bool**, **optional*) – if function should
           raise error instead of returning false if position is
           outside the safety limits. The default is False.

      Returns:
         **ret** – if position is within safety limits.

      Return type:
         bool

   intersect(popups=False, **kwargs)

      Check if meshes in the simulation are intersecting

      -[ Example ]-

      coliding_meshes=intersect()

      Parameters:
         * **popups** (*bool**, **optional*) – Draw popup window if
           collision is detected. The default is False.

         * **verbose** (*bool**, **optional*) – print every check. The
           default is False.

      Returns:
         **intersections** – list of intersecting mesh name pairs.

      Return type:
         list [[mesh1,mesh2],…,[mesh69,mesh420]]

   lastub()

      load first diffractomter state file it can find.

      -[ Example ]-

      lastub()

   latt()

      get current crystal lattice

      -[ Example ]-

      latt() #prints current latice

      Returns:
         **lattice** – Current lattice.

      Return type:
         [str,float,float,float,float,float,float]

   listub(ret=False)

      Print or return a list of diffractomter state files in the
      current directory.

      -[ Example ]-

      listub() #print all state files in the current directory

      Parameters:
         **ret** (*bool**, **optional*) – if return the list of files.
         The default is False.

      Returns:
         list of state files.

      Return type:
         list

   loadub(file)

      load a diffractomter state from a file. Can read .i16sim.txt or
      .nxs files. Loads UB, U, position, constraints, energy, crystal
      lattice, and azimuthal references if they are given.

      -[ Example ]-

      loadub(0) #load first state file it can find

      Parameters:
         **file** (*str** or **index:int*) – path or index to file.

   moveto(e_angles, use_limits=True, UI_call=True)

      Move the simulation and update the global state

      -[ Example ]-

      moveto([0,0,0,0,90,0]) #chi 90

      Parameters:
         * **e_angles** (*list** [**mu:float**, **delta:float**,
           **gamma:float**, **eta:float**, **chi:float**,
           **phi:float**]*) – eulerian rotations.

         * **use_limits** (*bool**, **optional*) – if safety limits
           should be used if they are enabled. The default is True.

         * **UI_call** (*bool**, **optional*) – If UI update call
           should be tried. The default is True.

      Returns:
      Return type:
         None.

   new_hklcalc()

      Creates a new hkl calculation object

   newub(name)

      Starts a new UB calculation

      -[ Example ]-

      newub(‘example1’)

      Parameters:
         **name** (*str*) –

   pos(*args)

      The diffcalc pos command. If no arguments given, print the
      values of all scannables. elif the the only argument is a
      scannable, print its value. elif there are two arguments, a
      scannable and a value, set the scannable to that value. elif
      there are multiple scannable and value pairs, set each scannable
      to the corresponding value.

      -[ Examples ]-

      pos() #prints everything pos(sixc) #prints
      [phi,chi,eta,mu,delta,gamma] pos(hkl,[1,1,1]) #moves simulation
      to position hkl = [1,1,1]

      Parameters:
         ***args** (*[**scannable1**, **value1**, **scannable2**,
         **value2**, **...****]*) –

      Returns:
      Return type:
         None.

   pos_from_hkl(hkl)

      Get the closest position and virtual angles corresponding to a
      certain hkl. If limits are enabled, makes sure the position is
      within limits. Raises an exception if no possition.

      -[ Example ]-

      position_ob, virtual_angles = pos_from_hkl([1,1,1])

      Parameters:
         **hkl** (*[**float**,**float**,**float**]*) – Miller idex
         list.

      Returns:
         * **pos_best** (*Position*) – the optimal position for this
           hkl.

         * **va_best** (*TYPE*) – the corresponding virtual angles.

   pos_to_sixc(e_angles)

      Parameters:
         **e_angles** (*[**mu:float**, **delta:float**,
         **gamma:float**, **eta:float**, **chi:float**,
         **phi:float**]*) – diffcalc order

      Returns:
         **[phi** – i16 sixc order.

      Return type:
         float, chi:float, eta:float, mu:float, delta:float,
         gamma:float]

   read_visual_pos()

      Calculate motor rotations from visual roatation of meshes

      Returns:
         **[mu** – eulerian angles as seen on screen. Might be
         different from emulator position if desynced.

      Return type:
         float, delta:float, gamma:float, eta:float, chi:float,
         phi:float]

   scan(*args, w=0.01)

      scan ‘scannable’ from ‘start’ to ‘stop’ in increments of ‘step’.
      .. rubric:: Example

      scan(eta,0,5,1) #moves to eta = 0, 1, 2, 3, 4, and 5 in turn

      Parameters:
         * **[****scannable** (**args =*) –

         * **start** –

         * **stop** –

         * **step** –

         * **...****]** –

      Returns:
      Return type:
         None.

   scancn(*args, w=0.01)

      centered scan on current position. scan ‘scannable’ around the
      starting position by going to ‘numer_of_steps’ number of
      positions that are separated by value ‘step’ Then return to
      starting positoin.

      -[ Example ]-

      pos(eta,0) scancn(eta,1,5) # moves to positions eta = -2, -1, 0,
      1, and 2

      Parameters:
         * **[****scannable** (**args =*) –

         * **step** –

         * **number_of_steps** –

         * **...****]** –

      Returns:
      Return type:
         None.

   set_reciprocal_vectors(*args)

      Calculates and displays reciprocal lattice vectors in the lab
      frame

      -[ Example ]-

      set_reciprocal_vectors()

   setlat(*args)

      Set crystal lattice parameters using shortform notation.

      Following combinations of system and lattice parameters are
      supported:

      (‘Cubic’, a) – sets Cubic system (‘Tetragonal’, a, c) – sets
      Tetragonal system (‘Hexagonal’, a, c) – sets Hexagonal system
      (‘Orthorhombic’, a, b, c) – sets Orthorombic system
      (‘Rhombohedral’, a, alpha) – sets Rhombohedral system
      (‘Monoclinic’, a, b, c, beta) – sets Monoclinic system
      (‘Triclinic’, a, b, c, alpha, beta, gamma) – sets Triclinic
      system

      Crystal system can be inferred from the lattice parameters for
      the following cases:

      (a,) – assumes Cubic system (a, c) – assumes Tetragonal system
      (a, b, c) – assumes Orthorombic system (a, b, c, angle) –
      assumes Monoclinic system with beta not equal to 90 or

         Hexagonal system if a = b and gamma = 120

      (a, b, c, alpha, beta, gamma) – sets Triclinic system

      -[ Example ]-

      #5 Angstrom base and all lattice vectors are at right angles.
      setlat(‘orthorombic 5x5x10 A’, 5, 5, 10, 90, 90, 90)

      Parameters:
         * **name** (*str*) – Crystal name

         * **system** (*Optional**[**float**]**, **default = None*) –
           Crystal lattice type.

         * **a** (*Optional**[**float**]**, **default = None*) –
           Crystal lattice parameter.

         * **b** (*Optional**[**float**]**, **default = None*) –
           Crystal lattice parameter.

         * **c** (*Optional**[**float**]**, **default = None*) –
           Crystal lattice parameter.

         * **alpha** (*Optional**[**float**]**, **default = None*) –
           Crystal lattice angle.

         * **beta** (*Optional**[**float**]**, **default = None*) –
           Crystal lattice angle.

         * **gamma** (*Optional**[**float**]**, **default = None*) –
           Crystal lattice angle.

   setlm(name, vals)

      set the limits of a scannable of composite limit

      -[ Example ]-

      setlm(eta,[None,90]) #the safety limits for eta are now ‘eta <
      90’

      Parameters:
         * **name** (*Scannable**, **str*) – an id for a scannable or
           composite limit

         * **vals** (*float** or **list of floats*) – new limit
           values.

      Returns:
      Return type:
         None.

   setnhkl(vector=None)

      Return reference vector property represented using miller
      indices. If a parameter is given, set the vector to the given
      value.

      Parameters:
         **vector** (*[**float**,**float**,**float**]**, **optional*)
         – The value to be set. The default is None.

      Returns:
         Reference vector represented as (3,1) NumPy array.

      Return type:
         np.ndarray

   setnphi(vector=None)

      Return reference vector property represented using sample frame
      coordinates. If a parameter is given, set the vector to the
      given value.

      Parameters:
         **vector** (*[**float**,**float**,**float**]**, **optional*)
         – The value to be set. The default is None.

      Returns:
         Reference vector represented as (3,1) NumPy array.

      Return type:
         np.ndarray

   showlm(*args)

      Prints limits. If no parameters given, prints all safety limits
      If for every parameter, prints its limits.

      -[ Example ]-

      showlm(eta) # shows the limits and cut of eta.

      Parameters:
         * **[****id1** (**args =*) – scannables, composite limits, or
           their ids.

         * **id2** – scannables, composite limits, or their ids.

         * **id3** – scannables, composite limits, or their ids.

         * **.****]** – scannables, composite limits, or their ids.

      Returns:
      Return type:
         None.

   showorient()

      shows added orientations

      -[ Example ]-

      showorient() #prints a table of orientations

   showref()

      shows added reflections

      -[ Example ]-

      showref() #prints a table of reflections

   sim(scannable, value)

      Simulates a move to a certain sixc or hkl position. prints all
      the position values of that position. Checks for collisions and
      errors if trying to move to that positoin.

      -[ Example ]-

      sim(hkl,[1,1,1]) #predicts what would happen if moving to hkl =
      [1,1,1]

      Parameters:
         * **scannable** (*Scannable** or **str*) – hkl or sixc (or
           their keys)

         * **value** (*list of floats*) – corresponding position to
           simulate.

      Returns:
      Return type:
         None.

   sixc_to_pos(sixc_e_angles)

      Parameters:
         **sixc_e_angles** (*[**phi:float**, **chi:float**,
         **eta:float**, **mu:float**, **delta:float**,
         **gamma:float**]*) – i16 sixc order..

      Returns:
         **[mu** – diffcalc order

      Return type:
         float, delta:float, gamma:float, eta:float, chi:float,
         phi:float]

   surfnhkl(vector=None)

      Return surface vector property represented using miller indices.
      If a parameter is given, set the vector to the given value.

      Parameters:
         **vector** (*[**float**,**float**,**float**]**, **optional*)
         – The value to be set. The default is None.

      Returns:
         Reference vector represented as (3,1) NumPy array.

      Return type:
         np.ndarray

   surfnphi(vector=None)

      Return surface vector property represented using sample frame
      coordinates. If a parameter is given, set the vector to the
      given value.

      Parameters:
         **vector** (*[**float**,**float**,**float**]**, **optional*)
         – The value to be set. The default is None.

      Returns:
         Reference vector represented as (3,1) NumPy array.

      Return type:
         np.ndarray

   trialub(id=1)

      Estimate UB from one reflection

      -[ Example ]-

      trialub() # will approximate the UB matrix if one reflection was
      added

   ub()

      prints the state of the current UB calculation

      -[ Example ]-

      ub() #prints the ub table

      Returns:
         The UB matrix.

      Return type:
         numpy matrix (3,3)

   update_pos(command=None)

      updates position and state to synchronise UI and emulator.

      -[ Example ]-

      update_pos() #synchronises

      Parameters:
         **command** (*int**, **optional*) – 0 - go to sixc() ==
         [0,0,0,0,0,0] positoin. 90 - go to sixc() == [0,90,0,0,0,0]
         position. The default is None.

   varifykey(key)

      get the id of a scannable after comparing to alternative names
      and checking if this is not the object itself.

      -[ Example ]-

      varifykey(‘sixcircle’) #returns ‘sixc’ - the real scannable id

      Parameters:
         **key** (*str** or **Scannable*) – The scannable id or the
         scannable itself.

      Returns:
         **key** – The scannable id.

      Return type:
         str

* i16sim command reference diffcalc

* copy 2

  * i16sim c2

  * i16sim c3
