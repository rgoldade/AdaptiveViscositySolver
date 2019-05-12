# Adaptive Viscosity Solver Plugin for Houdini

## Build
To build this project in Houdini (Linux):

1. Install Houdini 17.0 or higher.

2. Go to install folder (/opt/hfs.xx).

3. Type "source houdini_setup" to get the necessary environment variables.

4. Make a build folder the top level of the repository.

5. Run cmake .. in the build folder (or cmake -DUSEEIGEN .. if Eigen 3 is preferred over the native Houdini solver).

6. Run make in the build folder.

7. Verify that it was added to Houdini by:
	7.1 Launch Houdini.
	7.2 Press "tab" in the Network Editor and select a "DOP Network" and place it in the editor.
	7.3 Jump into the DOP Network, press "tab" again and verify that "HDK Adaptive Viscosity" is searchable.

## Usage

1. Build your liquid simulation scene.

2. Unlock the FLIP Solver by right clicking on the DOP and selecting "Allow Editing of Contents".

3. Dive into the FLIP Solver, find the viscosity solver DOP called "gasviscosity" and disable it by selecting the "Bypass" radial menu option.

4. Place the "HDK Adaptive Viscosity" solver, wire it into the disabled "gasviscosity" DOP.


## Legal

Copyright 2019 Ryan Goldade

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

