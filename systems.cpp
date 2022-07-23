#include<math.h>

//a basic linear scalar ODE ideal for testing

float linearScalarODE(float t, float y) {
	float ydot = 1.0 * y + sin(t);
	return ydot;
}