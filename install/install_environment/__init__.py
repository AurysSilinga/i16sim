bl_info = {
    "name": "i16sim installation Add-on",
    "blender": (2, 80, 0),
    "category": "Object",
    "description": "Only use if the install_i16sim_environment.py script does not work"
}
"""Installs modules and settings for i16sim into Blender's addon folder. 
Only use if the install_i16sim_enviromnet.py script does not work. 
Install this addon via the 'preferences > add-ons' menu, open a terminal, 
and then enble this addon. Once it has finished executing, disbale this addon and restart Blender"""

import subprocess
import sys
import os
import shutil
import bpy

#print(bpy.utils.user_resource('SCRIPTS', "addons"))

packages=['numpy','scipy','h5py']
python_path=sys.executable
script_file = os.path.realpath(__file__)
directory = os.path.dirname(script_file) #the addons folder
addon_folder = os.path.dirname(directory)

def install(package,location):
    subprocess.check_call([sys.executable, "-m", "pip", 'install', '--target', location, package])

def register():
    print('Installing packages for Blender in',directory)
    subprocess.check_call([sys.executable, "-m", "ensurepip", ])
    for package in packages:
        install(package,addon_folder)
        pass
        
    
    print("Done installing modules")
    print("Installing user preferences at",directory+'/userpref.blend')
    config_folder=bpy.utils.user_resource('CONFIG')
    shutil.copy(directory+'/userpref.blend',config_folder)
    print()
    print("Finished installation")
    print("Please restart Blender")
    print()
    
def unregister():
    pass