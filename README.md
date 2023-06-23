# Bezier, BSpline Curves & Surfaces
Creating BSpline Curves &amp; Surfaces in OpenGL, computer graphics related

![BSpline 2](https://github.com/juandiegovil/BSplineCurves-Surfaces2/assets/66028457/030833b5-7c88-450f-81ec-e62c62c19704)

This program produces four different scenes, each with control points (not drawn in scenes 3 and 4) and Bezier and/or B-Spline.
The first scene is in 2D only, but you can move the camera around in the other scenes for 3D viewing. Scene 3 shows a 3d object made
by revolving the curve in scene 1 among the y-axis. Scene 4 shows 2 different B-Spline surfaces constructed from a given set of control points and an arbitrary (but fixed) set of control points.

## Running Code:
Open the 453-skeleton.exe file in \SpaceShipGame\out\build\x64-Debug folderpath

## Controls:
Switching Scenes: Use keys 1, 2, 3, and 4 to switch from scenes 1-4 respectively

Scene 1:
	
  - Mouse Left Click (and release): Adds a control point for the Bezier or B-Spline curve
	- Mouse Right Click (and release): Deletes the control point that the cursor is on top of
	- Mouse Left Click (hold) and drag: Selects the control point that the cursor is on top of and move it, following the cursor until released
	- Backspace: Deletes the entire set of control points in the scene
	- Spacebar: Toggles from Bezier curve to B-Spline and vice versa

Scene 2-4:
  - W Key: Zooms in to the focal point of the camera, will not zoom in if it has reached such point
	- S Key: Zooms out of the focal point of the camera
	- A Key: Rotates the camera left
	- D Key: Rotates the camera right
	- Mouse Left Click (hold) and drag: Moves the drawing along the window in the same direction as the cursor, done so by moving the location and focal point of the camera

Scene 2:

  - Spacebar: Toggles from Bezier curve to B-Spline and vice versa

Scene 3:
	
  - Spacebar: Toggles from wireframe to fill drawing mode

Scene 4:
  
  - Spacebar: Toggles from wireframe to fill drawing mode
  - Q: Toggles from the surface with the given set of control points and the surface with arbitrary values (fixed)

## Compiler and Platform
Compiler: Clang++

Build Tools: Cmake

![BSpline 1](https://github.com/juandiegovil/BSplineCurves-Surfaces2/assets/66028457/021a379d-8941-400c-a969-cccaa428da3a)
