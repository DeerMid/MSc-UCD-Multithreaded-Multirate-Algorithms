#include<math.h>

//a basic linear scalar ODE ideal for testing

float linearScalerODE(float t, float y) {
	float ydot = 1.0 * y + sin(t);
	return ydot;
}