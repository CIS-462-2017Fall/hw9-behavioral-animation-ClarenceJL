// Spark.cpp: implementation of the ASpark class.
//
//////////////////////////////////////////////////////////////////////

#include "aSpark.h"
#include <math.h>

const float PI = 3.14159265358;

#ifndef GRAVITY
#define GRAVITY 9.8f
#endif

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ASpark::ASpark()
{
	//coefficients of restitution equals 0.25
	m_COR = 0.25f;
	m_mass = 1.0;
}

ASpark::ASpark(float* color): AParticle()
{
	for (int i = 0; i < 3; i++)
		m_color[i] = color[i];
 
	//coefficients of restitution equals 0.25
	m_COR = 0.25f;
}

ASpark::~ASpark()
{

}

//Set attractor position
void ASpark::setAttractor(vec3 position)
{
	m_attractorPos = position;
}

//Set repeller position
void ASpark::setRepeller(vec3 position)
{
	m_repellerPos = position;
}

void ASpark::setWind(vec3 wind)
{
	m_windForce = wind;
}

void ASpark::display()
{
	float fadeTime = 3.0;
	if (m_alive)
	{
		float alpha = 1.0;
		if (m_state[10] < fadeTime)
		{
			alpha = m_state[10] / 10.0f;
		}
		float scale = 1.0;

		glPushMatrix();
		glColor4f(m_color[0], m_color[1], m_color[2], alpha);
		glTranslatef(m_state[0], m_state[1], m_state[2]);
		glScalef(scale, scale, scale);
		glutSolidSphere(1.0, 10, 10);
		glPopMatrix();
	}

}
	


void ASpark::update(float deltaT, int extForceMode)
{
	m_deltaT = deltaT;
	if (m_state[10] <= 0.0)
	{
		m_alive = false;
		return;
	}

	if (!(extForceMode & EXT_SPARKFORCES_ACTIVE))
		extForceMode = 0;
	
	computeForces(extForceMode);
	
	updateState(deltaT, EULER);

	resolveCollisions();
	
}


 
void ASpark::computeForces(int extForceMode)
//	computes the forces applied to this spark
{
	// zero out all forces
	m_state[6] = 0.0;
	m_state[7] = 0.0;
	m_state[8] = 0.0;

	// gravity force
	addForce(m_mass*m_gravity);


	// wind force
	if (extForceMode & WIND_ACTIVE)
	{
		//TODO: Add your code here
		addForce(m_windForce);

	}

	if (extForceMode & DRAG_ACTIVE)
	{
		//TODO: Add your code here
		float c = 1.5; // constant for drag force
		addForce(- c*m_Vel);
	}


	// attractor force
	if (extForceMode & ATTRACTOR_ACTIVE)
	{
		//TODO: Add your code here
		vec3 dir = m_attractorPos - m_Pos;
		float r = dir.Length();
		if (r < 1e-10) {
			r = 1e-10;
			dir = { 0,0,0 };
		}
		else
			dir.Normalize();
		float magnitude = 1000.0/(r*r);
		addForce(m_mass * magnitude * dir);
	
	}

	// repeller force
	if (extForceMode & REPELLER_ACTIVE)
	{
		//TODO: Add your code here
		vec3 dir = m_Pos - m_repellerPos;
		float r = dir.Length();
		if (r < 1e-10) {
			r = 1e-10;
			dir = { 0,0,0 };
		}
		else
			dir.Normalize();
		float magnitude = 1000.0/(r*r);
		addForce(m_mass * magnitude * dir);
	}


	// random force
	if (extForceMode & RANDOM_ACTIVE)
	{
		
		//TODO: Add your code here
		vec3 rand_dir;
		float theta = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * PI * 2.0; //random force direction between 0~2*pi
		rand_dir[0] = cos(theta);
		rand_dir[1] = sin(theta);
		rand_dir[2] = 0.0;
		/*
		//3D random direction
		float z = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 2.0 - 1.0;
		rand_dir[0] = sqrt(1-z*z)*cos(theta);
		rand_dir[1] = sqrt(1-z*z)*sin(theta);
		rand_dir[2] = z;
		rand_dir.Normalize();
		*/
		float rand_magnitude = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 250; //random force magnitude between 0~250
		addForce(rand_magnitude * rand_dir);
		
	}

}

void ASpark::resolveCollisions()
// resolves collisions of the spark with the ground
{
	//TODO: Add  code here that reverses the y value of the spark velocity vector 
	//      when the y position value of the spark is < 0
	if(m_state[1] < 0){
		m_state[4] = - m_COR*m_state[4]; //turn the velocity in y direction
		m_Vel[1] = m_state[4];
		m_stateDot[1] = m_state[4];
	}

}
