#define _USE_MATH_DEFINES
#include <glad/glad.h>
//#include <GLFW/glfw3.h>
#include <cmath>
#include <glad/glad.h>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"


// We gave this code in one of the tutorials, so leaving it here too
void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry const &cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);
}

// EXAMPLE CALLBACKS
class Assignment3 : public CallbackInterface {

public:
	Assignment3(int screenWidth, int screenHeight) :
		screenDim(screenWidth, screenHeight)
	{
	}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (action == GLFW_PRESS || action == GLFW_REPEAT) {
			if (key == GLFW_KEY_1) { //Curves in 2D with Clearing Option
				scene = 1;
			}
			else if (key == GLFW_KEY_2) { //In 2D with no editing, Camera movement. With both bezier and spline curves
				scene = 2;
			}
			else if (key == GLFW_KEY_3) { //3D rotation with camera movement, only B-spline
				scene = 3;
			}
			else if (key == GLFW_KEY_4) { //Surface thing with camera movement and 2 iterations, one given, one arbitrary. B-spline
				scene = 4;
			}

			if (key == GLFW_KEY_BACKSPACE) {
				pointState = 4; //might need to input some geometries for this, or maybe just clear shader and do some shit (main startup, change a reset value in main for reset starting code
			}
			if (scene > 1) { //Do nothing when movement is 0
				if (key == GLFW_KEY_W) { //UP
					movement = 1;
				}
				else if (key == GLFW_KEY_S) { //DOWN
					movement = 2;
				}
				else if (key == GLFW_KEY_A) { //LEFT
					movement = 3;
				}
				else if (key == GLFW_KEY_D) { //RIGHT
					movement = 4;
				}
			}
			if (key == GLFW_KEY_SPACE) {
				if (scene > 2) { //TOGGLING WIREFRAME & FULL 
					if (wireframe) {
						wireframe = false;
					}
					else {
						wireframe = true;
					}
				}
				else { //TOGGLING BEZIER & BSPLINE 
					if (spline) {
						spline = false;
					}
					else {
						spline = true;
					}
				}
			}
			if (key == GLFW_KEY_Q && scene == 4) {//TOGGLING TWO SURFACES
				if (surf) {
					surf = false;
				}
				else {
					surf = true;
				}
			}
		}
	}

	virtual void mouseButtonCallback(int button, int action, int mods) {
		if (scene == 1) {
			if (button == GLFW_MOUSE_BUTTON_LEFT && action != GLFW_RELEASE) { //Move point
				pointState = 1;
			}
			else if (button == GLFW_MOUSE_BUTTON_LEFT) { //New point
				pointState = 2;
			}
			if (button == GLFW_MOUSE_BUTTON_RIGHT && action != GLFW_RELEASE) { //Delete point
				pointState = 3;
			}
		}
		if (scene > 1) {
			if (button == GLFW_MOUSE_BUTTON_LEFT && action != GLFW_RELEASE) { //Mouse Released
				movement = 5;
			}
			else if (button == GLFW_MOUSE_BUTTON_LEFT) { //Drag camera around
				movement = 0;
			}
		}
	}

	glm::vec2 mouseGL() {
		glm::vec2 startingVec(xScreenPos, yScreenPos);
		glm::vec2 shiftedVec = startingVec + glm::vec2(0.5f, 0.5f);
		glm::vec2 scaledToZeroOne = shiftedVec / glm::vec2(screenDim);

		glm::vec2 flippedY = glm::vec2(scaledToZeroOne.x, 1.0f - scaledToZeroOne.y);

		glm::vec2 final = flippedY * 2.0f - glm::vec2(1.0f, 1.0f);

		return final;
	}

	virtual void cursorPosCallback(double xpos, double ypos) {
		xScreenPos = xpos;
		yScreenPos = ypos;
		cursorPosition = mouseGL();
	}
	virtual void scrollCallback(double xoffset, double yoffset) {
	}
	virtual void windowSizeCallback(int width, int height) {
		// The CallbackInterface::windowSizeCallback will call glViewport for us
		CallbackInterface::windowSizeCallback(width,  height);
	}

	int getState() {
		return pointState;
	}
	void updateState() {
		pointState = 0;
	}

	int getScene() {
		return scene;
	}

	int getMovement() {
		return movement;
	}
	void updateMovement() {
		movement = 0;
	}

	glm::vec2 getCursorPosition() {
		return cursorPosition;
	}

	bool getWireframe() {
		return wireframe;
	}

	bool getSpline() {
		return spline;
	}

	bool getSurf() {
		return surf;
	}


private:
	glm::ivec2 screenDim;
	glm::vec2 cursorPosition;
	double xScreenPos;
	double yScreenPos;
	int pointState = 0;
	int scene = 1;
	int movement = 0;
	bool wireframe = true;
	bool spline = false;
	bool surf = false;

};

void originalSquare(std::vector<glm::vec3> const& originQuad, CPU_Geometry& geom) {
	for (auto i = originQuad.begin(); i < originQuad.end(); ++i) {
		geom.verts.push_back(glm::vec3(*i));
	}
	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
}

//FACES OF SQUARE

//void positiveZFace(std::vector<glm::vec3> const& originQuad, CPU_Geometry& geom) {
//	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, 1.0f)); ///CHANGED IT TEST WITHOUT A Z
//	for (auto i = originQuad.begin(); i < originQuad.end(); ++i) {
//		geom.verts.push_back(glm::vec3(T * glm::vec4((*i), 1.0f)));
//	}
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//}
//
//void positiveXFace(std::vector<glm::vec3> const& originQuad, CPU_Geometry& geom) {
//	const glm::mat4 R = glm::rotate(glm::mat4(1.0f), glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
//	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(0.5, 0.0, 0.0f));
//	for (auto i = originQuad.begin(); i < originQuad.end(); ++i) {
//		geom.verts.push_back(glm::vec3(T * R * glm::vec4((*i), 1.0f)));
//	}
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//}
//
//void negativeZFace(std::vector<glm::vec3> const& originQuad, CPU_Geometry& geom) {
//	const glm::mat4 R = glm::rotate(glm::mat4(1.0f), glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
//	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, -0.5f));
//	for (auto i = originQuad.begin(); i < originQuad.end(); ++i) {
//		geom.verts.push_back(glm::vec3(T * R * glm::vec4((*i), 1.0f)));
//	}
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//}
//
//void negativeXFace(std::vector<glm::vec3> const& originQuad, CPU_Geometry& geom) {
//	const glm::mat4 R = glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
//	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(-0.5f, 0.0f, 0.0f));
//	for (auto i = originQuad.begin(); i < originQuad.end(); ++i) {
//		geom.verts.push_back(glm::vec3(T * R * glm::vec4((*i), 1.0f)));
//	}
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//}
//
//void positiveYFace(std::vector<glm::vec3> const& originQuad, CPU_Geometry& geom) {
//	const glm::mat4 R = glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
//	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.5f, 0.0f));
//	for (auto i = originQuad.begin(); i < originQuad.end(); ++i) {
//		geom.verts.push_back(glm::vec3(T * R * glm::vec4((*i), 1.0f)));
//	}
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//}
//
//void negativeYFace(std::vector<glm::vec3> const& originQuad, CPU_Geometry& geom) {
//	const glm::mat4 R = glm::rotate(glm::mat4(1.0f), glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
//	const glm::mat4 T = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, -0.5f, 0.0f));
//	for (auto i = originQuad.begin(); i < originQuad.end(); ++i) {
//		geom.verts.push_back(glm::vec3(T * R * glm::vec4((*i), 1.0f)));
//	}
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//	geom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
//	geom.cols.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
//}


float getDistance(glm::vec3 firstPoint, glm::vec3 secondPoint) {
	return sqrtf(pow((secondPoint.x - firstPoint.x), 2) + pow((secondPoint.y - firstPoint.y), 2));
}

auto checkClosestPoint(CPU_Geometry& geom, glm::vec3 pos) {
	auto i = geom.verts.begin();
	auto closest = i;
	float d = getDistance(*i, pos);
	i++;
	while (i != geom.verts.end()) {
		float newD = getDistance(*i, pos);
		if (newD < d) {
			closest = i;
			d = newD;
		}
		i++;
	}
	return closest;
}

void changePoints(int state, float time, CPU_Geometry& cpuGeom, glm::vec2 position, bool &stuck, std::vector<glm::vec3>::iterator iter, float d) {
	if (state == 1) {
		if (position.x > 0.99f) { //not going out of bounds
			position.x = 0.99f;
		}
		else if (position.x < -0.99f) {
			position.x = -0.99f;
		}
		if (position.y > 0.99f) {
			position.y = 0.99f;
		}
		else if (position.y < -0.99f) {
			position.y = -0.99f;
		}

		if (time > 0.25f) { //move point
			if (d < 0.03f || stuck) {
				stuck = true;
				glm::mat4 T1 = glm::translate(glm::mat4(1.0f), -(*iter));
				glm::mat4 T2 = glm::translate(glm::mat4(1.0f), glm::vec3(position, 0.0f));
				*iter = glm::vec3(T2 * T1 * glm::vec4((*iter), 1.0f));
			}
		}
	}
	else if (state == 2) {
		if (!(position.x > 0.99f || position.x < -0.99f || position.y > 0.99f || position.y < -0.99f)) { //add point
			cpuGeom.verts.push_back(glm::vec3(position, 0.0f));
			cpuGeom.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
		}
	}
	else if (state == 3) { //delete point
		if (d < 0.02f) {
			if (cpuGeom.verts.size() > 1) {
				cpuGeom.verts.erase(iter);
			}
			else {
				cpuGeom.verts.clear();
			}
		}
	}
}

void bezierCurve(CPU_Geometry & curveGeom, CPU_Geometry const & squareGeom) {
	//Do I need to do something with degree??
	
	std::vector<glm::vec3> square;

	for (auto iter = squareGeom.verts.begin(); iter != squareGeom.verts.end(); ++iter) {
		square.push_back(*iter);
	}
	std::vector<glm::vec3> origSquare = square;
	int d = square.size();

	//column i = the points in the square //an array of points
	for (float u = 0; u < 1; u += 0.01) {
		square = origSquare;
		for (int i = 1; i < d; i++) { //goes through every column dependent on degrees chosen
			for (int j = 0; j < d-i; j++) { //Does every simple 2 by one addition, puts in in d-i elements in column i
				square[j] = (1 - u) * (square[j]) + u * (square[j + 1]);
			}
		} //only p at 0 is the final iteration at the end
		//column i at 0 is what you want to plot, depending on degree is how accurate/how much learning
		//column i is saved in plot P(i(u)) = Point as the value of i(0) at a value of u
		
		curveGeom.verts.push_back(square[0]);
		curveGeom.cols.push_back(glm::vec3(0.0, 0.0, 1.0));
	}
}

std::vector<glm::vec3> innerBSpline(std::vector<glm::vec3> initial, std::vector<glm::vec3> next, int d) {
	for (float u = 0; u < 5; u++) {
		next[0] = initial[0];
		next[1] = 0.5f * initial[0] + 0.5f * initial[1];
		for (int i = 1; i < d - 2; i++) {
			next[2 * i] = 0.75f * initial[i] + 0.25f * initial[i + 1];
			next[2 * i + 1] = 0.25f * initial[i] + 0.75f * initial[i + 1];
		}
		next[2 * d - 4] = 0.5f * initial[d - 2] + 0.5f * initial[d - 1];
		next[2 * d - 3] = initial[d - 1];
		initial.clear();
		for (int j = 0; j < next.size(); j++) {
			initial.push_back(next[j]);
		}
		d = initial.size();
		next.resize(d * 2 - 2);
	}
	return initial;
}

void BSpline(CPU_Geometry& curveGeom, CPU_Geometry const& squareGeom) {

	std::vector<glm::vec3> square;

	for (auto iter = squareGeom.verts.begin(); iter != squareGeom.verts.end(); ++iter) {
		square.push_back(*iter);
	}
	int d = square.size();
	std::vector<glm::vec3> newSquare(d*2 - 2);

	square = innerBSpline(square, newSquare, d);

	for (int i = 0; i < square.size(); i++) { //copy whatever is in square
		curveGeom.verts.push_back(square[i]);
		curveGeom.cols.push_back(glm::vec3(0.0, 0.0, 1.0));
	}
}

glm::mat4 rotation(float angle) {
	return glm::rotate(glm::mat4(1.0f), glm::radians(angle), glm::vec3(0.0f, 1.0f, 0.0f));
}

int main() {
	Log::debug("Starting main");

	//std::cout << (*iter).x << "    " << (*iter).y << std::endl;

	// WINDOW
	glfwInit();
	Window window(800, 800, "CPSC 453"); // can set callbacks at construction if desired


	GLDebug::enable();

	// CALLBACKS
	auto a3 = std::make_shared<Assignment3>(800, 800);
	window.setCallbacks(a3);


	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	// The current CPU_Geometry and GPU_Geometry classes are defined in
	// Geometry.h/Geometry.cpp They will work for this assignment, but for some of
	// the bonuses you may have to modify them.
	std::vector<glm::vec3> originQuad;

	originQuad.push_back(glm::vec3{ -0.5, 0.5, 0 }); //top left
	originQuad.push_back(glm::vec3{ -0.5, -0.5, 0 }); //bottom left
	originQuad.push_back(glm::vec3{ 0.5, -0.5, 0 }); //bottom right
	originQuad.push_back(glm::vec3{ 0.5, 0.5, 0 }); //top right

	CPU_Geometry square;

	/*positiveXFace(originQuad, square);
	positiveYFace(originQuad, square);*/
	originalSquare(originQuad, square);
	/*
	negativeXFace(originQuad, square);
	negativeYFace(originQuad, square);
	negativeZFace(originQuad, square);*/

	CPU_Geometry ccurve;
	GPU_Geometry gcurve;

	bezierCurve(ccurve, square);
	updateGPUGeometry(gcurve, ccurve);

	CPU_Geometry revolvecCurve = ccurve;
	GPU_Geometry revolvegCurve;

	CPU_Geometry finalRevolvecCurve;
	GPU_Geometry finalRevolvegCurve;

	updateGPUGeometry(revolvegCurve, revolvecCurve);

	GPU_Geometry pointsGPUGeom;
	updateGPUGeometry(pointsGPUGeom, square);

	// Reset the colors to green
	square.cols.clear();
	square.cols.resize(square.verts.size(), glm::vec3{0.0, 1.0, 0.0});

	GPU_Geometry linesGPUGeom;
	updateGPUGeometry(linesGPUGeom, square);

	CPU_Geometry cSurfacePoints;
	GPU_Geometry gSurfacePoints;

	CPU_Geometry cSurface;
	GPU_Geometry gSurface;

	std::vector<glm::vec3> sur = std::vector{
		glm::vec3{-2.0f, 0.0f, -2.0f}, glm::vec3{-1.0f, 0.0f, -2.0f}, glm::vec3{0.0f, 0.0f, -2.0f}, glm::vec3{1.0f, 0.0f, -2.0f}, glm::vec3{2.0f, 0.0f, -2.0f},
		glm::vec3{-2.0f, 0.0f, -1.0f}, glm::vec3{-1.0f, 1.0f, -1.0f}, glm::vec3{0.0f, 1.0f, -1.0f}, glm::vec3{1.0f, 1.0f, -1.0f}, glm::vec3{2.0f, 0.0f, -1.0f},
		glm::vec3{-2.0f, 0.0f, 0.0f}, glm::vec3{-1.0f, 1.0f, 0.0f}, glm::vec3{0.0f, -1.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 0.0f}, glm::vec3{2.0f, 0.0f, 0.0f},
		glm::vec3{-2.0f, 0.0f, 1.0f}, glm::vec3{-1.0f, 1.0f, 1.0f}, glm::vec3{0.0f, 1.0f, 1.0f}, glm::vec3{1.0f, 1.0f, 1.0f}, glm::vec3{2.0f, 0.0f, 1.0f},
		glm::vec3{-2.0f, 0.0f, 2.0f}, glm::vec3{-1.0f, 0.0f, 2.0f}, glm::vec3{0.0f, 0.0f, 2.0f}, glm::vec3{1.0f, 0.0f, 2.0f}, glm::vec3{2.0f, 0.0f, 2.0f},
	};

	std::vector<glm::vec3> sur2 = std::vector{
		glm::vec3{-2.4f, 0.0f, -2.0f}, glm::vec3{-1.2f, -0.5f, -1.9f}, glm::vec3{0.0f, 1.0f, -1.8f}, glm::vec3{1.0f, 0.6f, -1.7f}, glm::vec3{1.6f, 0.3f, -1.6f},
		glm::vec3{-2.2f, 2.0f, -1.0f}, glm::vec3{-1.0f, -2.0f, -0.9f}, glm::vec3{-0.1f, 1.0f, -1.0f}, glm::vec3{1.1f, 0.8f, -1.2f}, glm::vec3{2.0f, 0.4f, -1.0f},
		glm::vec3{-1.8f, 0.9f, 0.0f}, glm::vec3{-1.1f, -1.5f, 0.0f}, glm::vec3{0.2f, -1.0f, 0.02f}, glm::vec3{1.0f, 1.0f, -0.1f}, glm::vec3{2.2f, 0.3f, 0.1f},
		glm::vec3{-1.9f, -0.5f, 1.0f}, glm::vec3{-1.3f, 0.3f, 1.2f}, glm::vec3{0.2f, 1.0f, 1.2f}, glm::vec3{0.9f, -2.4f, 1.0f}, glm::vec3{2.0f, -0.2f, 1.3f},
		glm::vec3{-2.0f, 0.0f, 2.0f}, glm::vec3{-1.0f, 0.2f, 2.1f}, glm::vec3{0.2f, 0.6f, 2.2f}, glm::vec3{0.8f, -1.0f, 2.3f}, glm::vec3{2.2f, -0.2f, 2.0f},
	};

	std::vector<glm::vec3> surface;
	for (int i = 0; i < sur.size(); i++) {
		surface.push_back(0.4f * sur[i]);
		cSurfacePoints.verts.push_back(surface[i]);
		cSurfacePoints.cols.push_back(glm::vec3{ 1.0f, 0.0f, 0.0f });
	}
	updateGPUGeometry(gSurfacePoints, cSurfacePoints);


	glPointSize(10.0f);

	// Note this call only work on some systems, unfortunately.
	// In order for them to work, you have to comment out line 60
	// If you're on a mac, you can't comment out line 60, so you
	// these will have no effect. :(
	// glLineWidth(5.0f);

	auto timeElapsed = glfwGetTime();
	float rotationAngle = 0.0f;
	float degreesPerSecond = 90.0f; //4 sec to go all the way around
	auto timeClicked = 0.0f;
	bool stuck = false;
	bool firstTime = true;
	bool bSpline = false;
	bool surf = false;
	bool firstTimeMove = true;
	float d = 0.5f;
	auto iter = checkClosestPoint(square, glm::vec3(a3->getCursorPosition(), 0.0f));
	glm::mat4 M = glm::mat4(1.0f);
	glm::mat4 V = glm::mat4(1.0f);
	glm::mat4 P = glm::mat4(1.0f);
	glm::vec3 cameraPosition = glm::vec3(0.0f, 0.0f, 0.5f);
	glm::vec3 centerPoint = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 cameraDirection = glm::normalize(cameraPosition - centerPoint);
	float cameraAngle = 0;
	int scene = 1;

	// RENDER LOOP
	while (!window.shouldClose()) {
		auto newTimeEleapsed = glfwGetTime();
		auto dt = newTimeEleapsed - timeElapsed;
		timeElapsed = newTimeEleapsed;

		glfwPollEvents();

		//resting other scene variables when scene is 1
		if (a3->getScene() == 1) {
			V = glm::mat4(1.0f);
			P = glm::mat4(1.0f);
			cameraPosition = glm::vec3(0.0f, 0.0f, 3.0f);
			centerPoint = glm::vec3(0.0f, 0.0f, 0.0f);
			scene = 1;
		}

		//changing from Bezier to BSpline
		if (a3->getScene() < 3 && bSpline != a3->getSpline()) {

			bSpline = a3->getSpline();
			ccurve.verts.clear();
			ccurve.cols.clear();
			if (square.verts.size() > 1) {
				if (!bSpline) {
					bezierCurve(ccurve, square);
				}
				else {
					BSpline(ccurve, square);
				}
			}
			updateGPUGeometry(gcurve, ccurve);

			revolvecCurve.verts.clear();
			revolvecCurve.cols.clear();

			revolvecCurve = ccurve;
			updateGPUGeometry(revolvegCurve, revolvecCurve);

		}

		//State of cursor on points & changing points
		if (a3->getState() > 0 && a3->getScene() == 1) {

			int pointState = a3->getState();
			if (pointState == 1 && square.verts.size() > 0) { //move point
				if (firstTime) { //sticks to point
					iter = checkClosestPoint(square, glm::vec3(a3->getCursorPosition(), 0.0f));
					d = getDistance(*iter, glm::vec3(a3->getCursorPosition(), 0.0f));
					firstTime = false;
				}
				timeClicked += dt;
			}
			else if (pointState == 2) {

				if (timeClicked > 0.25f) { //don't add point
					pointState = 0;
				}
				else if (square.verts.size() == 0) {
					square.verts.push_back(glm::vec3(a3->getCursorPosition(), 0.0f));
					square.cols.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
					pointState = 0;
				}

				timeClicked = 0.0f;
				stuck = false;
				firstTime = true;
				a3->updateState();
			}
			else if (pointState == 3) { //delete point
				a3->getState();
				iter = checkClosestPoint(square, glm::vec3(a3->getCursorPosition(), 0.0f));
				d = getDistance(*iter, glm::vec3(a3->getCursorPosition(), 0.0f));
				timeClicked = 0.0f;
				stuck = false;
				firstTime = true;
				a3->updateState();
			}
			else if (pointState == 4) { //delete point
				square.verts.clear();
				square.cols.clear();
				ccurve.verts.clear();
				ccurve.cols.clear();
				timeClicked = 0.0f;
				stuck = false;
				firstTime = true;
				a3->updateState();
			}


			square.cols.clear();
			square.cols.resize(square.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
			
			changePoints(pointState, timeClicked, square, a3->getCursorPosition(), stuck, iter, d);
			updateGPUGeometry(pointsGPUGeom, square);

			square.cols.clear();
			square.cols.resize(square.verts.size(), glm::vec3{ 0.0, 1.0, 0.0 });

			updateGPUGeometry(linesGPUGeom, square);

			ccurve.verts.clear();
			ccurve.cols.clear();

			if (square.verts.size() > 1) {
				if (!bSpline) {
					bezierCurve(ccurve, square);
				}
				else {
					BSpline(ccurve, square);
				}
			}
			updateGPUGeometry(gcurve, ccurve);

			revolvecCurve.verts.clear();
			revolvecCurve.cols.clear();

			revolvecCurve = ccurve;
			updateGPUGeometry(revolvegCurve, revolvecCurve);

			square.cols.clear();
			square.cols.resize(square.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
		}

		//reseting camera % curves when changing scene
		if (scene != a3->getScene() && a3->getScene() > 1) {
			scene = a3->getScene();

			//reseting camera when changing scene
			V = glm::lookAt(
				cameraPosition, //camera position
				centerPoint, //point to center at
				glm::vec3(0.0f, 1.0f, 0.0f));//up axis
			P = glm::perspective(glm::radians(45.0f), 800.0f / 800.0f, 0.1f, 10.0f);

			//doing 3D object once
			if (scene == 3) {
				finalRevolvecCurve.verts.clear();
				finalRevolvecCurve.cols.clear();
				if (!bSpline) {
					bSpline = true;

					ccurve.verts.clear();
					ccurve.cols.clear();
					square.cols.clear();
					square.cols.resize(square.verts.size(), glm::vec3{ 0.0, 1.0, 0.0 });

					BSpline(ccurve, square);
					updateGPUGeometry(gcurve, ccurve);

					revolvecCurve.verts.clear();
					revolvecCurve.cols.clear();
					revolvecCurve = ccurve;
					updateGPUGeometry(revolvegCurve, revolvecCurve);
				}

				std::vector<glm::vec3> line;
				for (auto iter = ccurve.verts.begin(); iter != ccurve.verts.end(); ++iter) {
					line.push_back(*iter);
				}
				std::vector<glm::vec3> newLine;

				for (int i = 0; i < 91; i++) {
					for (int j = 0; j < line.size(); j++) {
						revolvecCurve.verts.push_back(rotation(4.0f*i) * glm::vec4(line[j], 1.0f));
						revolvecCurve.cols.push_back(glm::vec3{ 0.0f, 0.0f, 1.0f });
					}
				}

				std::vector<int> h = { 0, 1, ((int)ccurve.verts.size()), 1, (1 + (int)ccurve.verts.size()), ((int)ccurve.verts.size()) };
				int pos = 0;
				for (int i = 0; i < 90; i++) {
					pos = pos + ccurve.verts.size();
					for (int j = 0; j < ccurve.verts.size() - 1; j++) {
						for (int k = 0; k < h.size(); k++) {
							finalRevolvecCurve.verts.push_back(revolvecCurve.verts[pos + h[k] + j]);
							finalRevolvecCurve.cols.push_back(glm::vec3{ 0.0f, 0.0f, 1.0f });
						}
					}
				}
				updateGPUGeometry(finalRevolvegCurve, finalRevolvecCurve);
			}

			if (scene == 4) {

				std::vector<glm::vec3> finalPoints;
				for (int i = 0; i < 5; i++) {
					std::vector<glm::vec3> points = { surface.begin() + 5*i, surface.begin() + 5*(i+1)};
					std::vector<glm::vec3> next(5 * 2 - 2);

					points = innerBSpline(points, next, 5);
					for (int j = 0; j < points.size(); j++) { //copy whatever is in square
						finalPoints.push_back(points[j]);
					}
				}

				std::vector<glm::vec3> resultPoints;
				int size = 0.2f*finalPoints.size();
				//int nextSize;
				for (int i = 0; i < size; i++) {
					std::vector<glm::vec3> points = { finalPoints[i], finalPoints[i + size], finalPoints[i + 2 * size], finalPoints[i + 3 * size], finalPoints[i + 4 * size] };
					std::vector<glm::vec3> next(5 * 2 - 2);
					points = innerBSpline(points, next, 5);

					for (int j = 0; j < size; j++) {
						resultPoints.push_back(points[j]);
					}

				}

				std::vector<int> h = { 0, 1, size, 1, (1 + size), (size) };
				int pos = 0;
				for (int i = 0; i < size - 1; i++) {
					for (int j = 0; j < size - 1; j++) {
						for (int k = 0; k < h.size(); k++) {
							cSurface.verts.push_back(resultPoints[pos + h[k] + j]);
							cSurface.cols.push_back(glm::vec3{ 0.0f, 0.0f, 1.0f });
						}
					}
					pos = pos + size;
				}
				
				updateGPUGeometry(gSurface, cSurface);
			}
		}

		//changing surface for scene
		if (surf != a3->getSurf() && a3->getScene() == 4) {
			surf = a3->getSurf();
			scene = 6;
			cSurfacePoints.verts.clear();
			cSurfacePoints.verts.clear();
			cSurface.verts.clear();
			cSurface.verts.clear();

			surface.clear();
			for (int i = 0; i < sur.size(); i++) {
				if (surf) {
					surface.push_back(0.4f * sur2[i]);
				}
				else {
					surface.push_back(0.4f * sur[i]);
				}
				cSurfacePoints.verts.push_back(surface[i]);
				cSurfacePoints.cols.push_back(glm::vec3{ 1.0f, 0.0f, 0.0f });
			}
			updateGPUGeometry(gSurfacePoints, cSurfacePoints);
		}

		//changing camera based on input
		if (a3->getScene() > 1 ) {
			int move = a3->getMovement();
			float tooCloseTooCenter = sqrtf(pow((cameraPosition[0]), 2) + pow((cameraPosition[2]), 2));
			if (move == 1 && tooCloseTooCenter > 0.05f) { //zoom in
				cameraPosition[0] = cameraPosition[0] - 0.05f * cameraDirection[0];
				cameraPosition[2] = cameraPosition[2] - 0.05f * cameraDirection[2];
				a3->updateMovement();
			}
			else if (move == 2) { //zoom out
				cameraPosition[0] = cameraPosition[0] + 0.05f * cameraDirection[0];
				cameraPosition[2] = cameraPosition[2] + 0.05f * cameraDirection[2];
				a3->updateMovement();
			}
			else if (move == 3) { //rotate left
				cameraPosition = rotation(-2.0f) * glm::vec4(cameraPosition, 1.0f);
				a3->updateMovement();
			}
			else if (move == 4) { //rotate right
				cameraPosition = rotation(2.0f) * glm::vec4(cameraPosition, 1.0f);
				a3->updateMovement();
			}
			else if (move == 5) { //move with cursor
				glm::vec2 initialPosition;
				glm::vec3 initialCameraPosition;
				glm::vec3 initialCenterPoint;
				if (firstTimeMove) {
					initialPosition = a3->getCursorPosition();
					initialCameraPosition = cameraPosition;
					initialCenterPoint = centerPoint;
					firstTimeMove = false;
				}
				else {
					glm::vec2 positionChange = a3->getCursorPosition() - initialPosition;
					glm::vec3 position;
					if (initialCameraPosition[2] > 0){
						position = glm::vec3(positionChange.x * sin(cameraAngle), -positionChange.y, positionChange.x * cos(cameraAngle));
					}
					else {
						position = glm::vec3(-positionChange.x * sin(cameraAngle), -positionChange.y, -positionChange.x * cos(cameraAngle));
					}
					if (cameraAngle < 0) {
						position[2] = -position[2];
					}
					else {
						position[0] = -position[0];
					}
					cameraPosition = glm::vec3(initialCameraPosition + position);
					centerPoint = glm::vec3(initialCenterPoint + position);
				}
			}
			else {
				firstTimeMove = true;
			}
			cameraDirection = glm::normalize(cameraPosition - centerPoint);
			cameraAngle = atan(cameraDirection[2] / cameraDirection[0]);
			V = glm::lookAt(
				cameraPosition, //camera position
				centerPoint, //point to center at
				glm::vec3(0.0f, 1.0f, 0.0f));//up axis

		}

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);

		if ((a3->getWireframe())) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		else {
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		shader.use();

		GLint myMLoc = glGetUniformLocation(shader, "M");
		glUniformMatrix4fv(myMLoc, 1, GL_FALSE, glm::value_ptr(M));
		GLint myVLoc = glGetUniformLocation(shader, "V");
		glUniformMatrix4fv(myVLoc, 1, GL_FALSE, glm::value_ptr(V));
		GLint myPLoc = glGetUniformLocation(shader, "P");
		glUniformMatrix4fv(myPLoc, 1, GL_FALSE, glm::value_ptr(P));

		if (scene == 1 || scene == 2) {
			linesGPUGeom.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(square.verts.size()));

			pointsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(square.verts.size()));
		}

		if (scene < 3) {
			gcurve.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(ccurve.verts.size()));
		}
		else if (scene == 3) {
			finalRevolvegCurve.bind();
			glDrawArrays(GL_TRIANGLES, 0, GLsizei(finalRevolvecCurve.verts.size()));
		}
		else {
			gSurfacePoints.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(cSurfacePoints.verts.size()));

			gSurface.bind();
			glDrawArrays(GL_TRIANGLES, 0, GLsizei(cSurface.verts.size()));
		}

		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}
