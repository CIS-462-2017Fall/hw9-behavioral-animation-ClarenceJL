#include "aBehaviors.h"

#include <math.h>
#include "GL/glew.h"
#include "GL/glut.h"


#define PI 3.14159265358979f

// Base Behavior
///////////////////////////////////////////////////////////////////////////////
Behavior::Behavior()
{
}

Behavior::Behavior( char* name) 
{
	m_name = name;
	m_pTarget = NULL;
}

Behavior::Behavior( Behavior& orig) 
{
	m_name = orig.m_name;
	m_pTarget = NULL;
}

string& Behavior::GetName() 
{
    return m_name;
}

// Behaviors derived from Behavior
//----------------------------------------------------------------------------//
// Seek behavior
///////////////////////////////////////////////////////////////////////////////
// For the given the actor, return a desired velocity in world coordinates
// Seek returns a maximum velocity towards the target
// m_pTarget contains target world position
// actor.getPosition() returns Agent's world position


Seek::Seek( AJoint* target) 
{
	m_name = "seek";
	m_pTarget = target;

}

Seek::Seek( Seek& orig) 
{
	m_name = "seek";
	m_pTarget = orig.m_pTarget;
}


Seek::~Seek()
{
}

vec3 Seek::calcDesiredVel( BehaviorController* actor)
{
	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	vec3 targetPos = m_pTarget->getLocalTranslation();
	vec3 actorPos = actor->getPosition();

	// TODO: add your code here to compute Vdesired
	vec3 err = targetPos - actorPos;
	double thetad = atan2(err[0], err[2]);
	double speed = BehaviorController::gMaxSpeed;
	Vdesired[0] = speed * sin(thetad); Vdesired[2] = speed * cos(thetad);

	return Vdesired;
}


// Flee behavior
///////////////////////////////////////////////////////////////////////////////
// For the given the actor, return a desired velocity in world coordinates
// Flee calculates a a maximum velocity away from the target
// m_pTarget contains target world position
// actor.getPosition() returns Agent's world position

Flee::Flee( AJoint* target) 
{
	m_name = "flee";
	m_pTarget = target;
}

Flee::Flee( Flee& orig) 
{
	m_name = "flee";
	m_pTarget = orig.m_pTarget;
}

Flee::~Flee()
{
}

vec3 Flee::calcDesiredVel( BehaviorController* actor)
{
	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	vec3 targetPos = m_pTarget->getLocalTranslation();
	vec3 actorPos = actor->getPosition();

	// TODO: add your code here to compute Vdesired
	vec3 err = targetPos - actorPos;
	double thetad = atan2(-err[0], -err[2]);
	double speed = BehaviorController::gMaxSpeed;
	Vdesired[0] = speed * sin(thetad); Vdesired[2] = speed * cos(thetad);

	return Vdesired;

}

// Arrival behavior
///////////////////////////////////////////////////////////////////////////////
// Given the actor, return a desired velocity in world coordinates
// Arrival returns a desired velocity vector whose speed is proportional to
// the actors distance from the target
// m_pTarget contains target world position
// actor.getPosition() returns Agent's world position
//  Arrival strength is in BehavioralController::KArrival


Arrival::Arrival( AJoint* target) 
{
	m_name = "arrival";
	m_pTarget = target;
}

Arrival::Arrival( Arrival& orig) 
{
	m_name = "arrival";
	m_pTarget = orig.m_pTarget;
}

Arrival::~Arrival()
{
}

vec3 Arrival::calcDesiredVel( BehaviorController* actor)
{
	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	vec3 targetPos = m_pTarget->getLocalTranslation();
	vec3 actorPos = actor->getPosition();

	// TODO: add your code here to compute Vdesired
	vec3 err = targetPos - actorPos;
	Vdesired = BehaviorController::KArrival * err;


	return Vdesired;
}


// Departure behavior
///////////////////////////////////////////////////////////////////////////////
// Given the actor, return a desired velocity in world coordinates
// Arrival returns a desired velocity vector whose speed is proportional to
// 1/(actor distance) from the target
// m_pTarget contains target world position
// actor.getPosition() returns Agent's world position
//  Departure strength is in BehavioralController::KDeparture

Departure::Departure(AJoint* target) 
{
	m_name = "departure";
	m_pTarget = target;
}

Departure::Departure( Departure& orig) 
{
	m_name = "departure";
	m_pTarget = orig.m_pTarget;
}

Departure::~Departure()
{
}

vec3 Departure::calcDesiredVel( BehaviorController* actor)
{
	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	vec3 targetPos = m_pTarget->getLocalTranslation();
	vec3 actorPos = actor->getPosition();

	// TODO: add your code here to compute Vdesired
	vec3 err = targetPos - actorPos;
	Vdesired = -BehaviorController::KDeparture * err / (err.Length()*err.Length());

	return Vdesired;
}


// Avoid behavior
///////////////////////////////////////////////////////////////////////////////
//  For the given the actor, return a desired velocity in world coordinates
//  If an actor is near an obstacle, avoid adds a normal response velocity to the 
//  the desired velocity vector computed using arrival
//  Agent bounding sphere radius is in BehavioralController::radius
//  Avoidance parameters are  BehavioralController::TAvoid and BehavioralController::KAvoid

Avoid::Avoid(AJoint* target, vector<Obstacle>* obstacles) 
{
	m_name = "avoid";
	m_pTarget = target;
	mObstacles = obstacles;
}

Avoid::Avoid( Avoid& orig) 
{
	m_name = "avoid";
	m_pTarget = orig.m_pTarget;
	mObstacles = orig.mObstacles;
}

Avoid::~Avoid()
{
}

vec3 Avoid::calcDesiredVel( BehaviorController* actor)
{

	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	m_actorPos = actor->getPosition();
	m_actorVel = actor->getVelocity();

	//TODO: add your code here

	// Step 1. compute initial value for Vdesired = Varrival so agent moves toward target
	Arrival arr_obj(m_pTarget);
	vec3 Varrival = arr_obj.calcDesiredVel(actor);

	//TODO: add your code here to compute Vavoid 
	
	// Step 2. compute Lb
	//TODO: add your code here
	double Lb = BehaviorController::TAvoid * m_actorVel.Length();

	// Step 3. find closest obstacloe 
	//TODO: add your code here
	vector<Obstacle>& obstacles = *mObstacles;
	double minDist = 1e10; 
	int minIdx = -1;
	vec3 d;
	for(int i = 0; i < obstacles.size(); ++i){
		m_obstaclePos = obstacles[i].m_Center.getLocalTranslation();
		d = m_obstaclePos - m_actorPos;
		double dist = d.Length();
		if(Dot(d, m_actorVel) > 0 && dist < minDist){
			minDist = dist;
			minIdx = i;
		}
	}
	if (minIdx == -1)
		return Varrival;//no obstacle ahead of agent

	m_obstaclePos = obstacles[minIdx].m_Center.getLocalTranslation();
	d = m_obstaclePos - m_actorPos;

	// Step 4. determine whether agent will collide with closest obstacle (only consider obstacles in front of agent)
	//TODO: add your code here
	vec3 actorVelDir = m_actorVel; 
	actorVelDir.Normalize(); //agent's moving direction (unit length)
	double dz_len = d * actorVelDir, dx_len = sqrt(d.Length()*d.Length()-dz_len*dz_len);
	double rb = actor->gAgentRadius, ro = obstacles[minIdx].m_Radius;
	bool collide = false;
	if(dz_len <= Lb + rb + ro && dx_len <= rb + ro)
		collide = true;

	// Step 5.  if potential collision detected, compute Vavoid and set Vdesired = Varrival + Vavoid
	//TODO: add your code here
	vec3 Vavoid(0, 0, 0);
	if(collide){
		//compute normal avoidance velocity
		vec3 ey(0, 1, 0);
		Vavoid = - (actorVelDir ^ ey);//direction of Vavoid
		Vavoid *= BehaviorController::KAvoid * (ro+rb-dx_len) / (ro+rb); //magnitude of Vavoid
		/*
		//or compute radial avoidance velocity
		vec3 n = Lb * actorVelDir - d;
		double nlen = n.Length();
		n.Normalize();
		Vavoid = n;
		if (rb + ro > nlen)
			Vavoid *= BehaviorController::KAvoid * (ro + rb - nlen) / (ro + rb);
		else
			Vavoid *= 0;
		*/
	}

	Vdesired = Varrival + Vavoid;

	return Vdesired;
	
}

void Avoid::display( BehaviorController* actor)
{
	//  Draw Debug info
	vec3 angle = actor->getOrientation();
	vec3 vel = actor->getVelocity();
	vec3 dir = vec3(cos(angle[1]), 0, sin(angle[1]));
	vec3 probe = dir * (vel.Length()/BehaviorController::gMaxSpeed)*BehaviorController::TAvoid;
	
	glBegin(GL_LINES);
	glColor3f(0, 0, 1);
	glVertex3f(m_actorPos[0], m_actorPos[1], m_actorPos[2]);
	glVertex3f(m_obstaclePos[0], m_obstaclePos[1], m_obstaclePos[2]);
	glColor3f(0, 1, 1);
	glVertex3f(m_actorPos[0], m_actorPos[1], m_actorPos[2]);
	glVertex3f(m_actorPos[0] + probe[0], m_actorPos[1] + probe[1], m_actorPos[2] + probe[2]);
	glEnd();
}


// Wander Behavior
///////////////////////////////////////////////////////////////////////////////
// For the given the actor, return a desired velocity in world coordinates
// Wander returns a desired velocity vector whose direction changes at randomly from frame to frame
// Wander strength is in BehavioralController::KWander

Wander::Wander() 
{
	m_name = "wander";
	m_Wander = vec3(1.0, 0.0, 0.0);
}

Wander::Wander( Wander& orig) 
{
	m_name = "wander";
	m_Wander = orig.m_Wander;
}

Wander::~Wander()
{
}

vec3 Wander::calcDesiredVel( BehaviorController* actor)
{
	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	vec3 actorPos = actor->getPosition();


	// compute Vdesired = Vwander

	// Step. 1 find a random direction
	//TODO: add your code here
  	double theta_rand = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 2 * PI;
  	vec3 dir_rand(cos(theta_rand),0,sin(theta_rand));

	// Step2. scale it with a noise factor
	//TODO: add your code here
  	vec3 r_noise = BehaviorController::KNoise * dir_rand;

	// Step3. change the current Vwander to point to a random direction
	//TODO: add your code here

	vec3 Vwander = BehaviorController::KWander * m_Wander; //let m_Wander store the direction of the previous v_wander
  	Vwander = (Vwander + r_noise).Normalize();
  	m_Wander = Vwander;

	// Step4. scale the new wander velocity vector and add it to the nominal velocity
	//TODO: add your code here
	vec3 Vnominal(1, 0, 0);
  	Vwander *= BehaviorController::KWander;
  	Vdesired = Vnominal + Vwander;

	return Vdesired;
}


// Alignment behavior
///////////////////////////////////////////////////////////////////////////////
// For the given the actor, return a desired velocity vector in world coordinates
// Alignment returns the average velocity of all active agents in the neighborhood
// agents[i] gives the pointer to the ith agent in the environment
// Alignment parameters are in BehavioralController::RNeighborhood and BehavioralController::KAlign


Alignment::Alignment(AJoint* target, vector<AActor>* agents) 
{
	m_name = "alignment";
	m_pAgentList = agents;
	m_pTarget = target;
}



Alignment::Alignment( Alignment& orig) 
{
	m_name = orig.m_name;
	m_pAgentList = orig.m_pAgentList;
	m_pTarget = orig.m_pTarget;

}

Alignment::~Alignment()
{
}

vec3 Alignment::calcDesiredVel( BehaviorController* actor)
{
	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	vec3 targetPos = m_pTarget->getLocalTranslation();
	vec3 actorPos = actor->getPosition();
	vector<AActor>& agentList = *m_pAgentList;
	

	// compute Vdesired 
	
	// Step 1. compute value of Vdesired for fist agent (i.e. m_AgentList[0]) using an arrival behavior so it moves towards the target
	BehaviorController* leader = agentList[0].getBehaviorController(); // first agent is the leader
	//TODO: add your code here
	if (actor == leader) {
		Arrival arr_obj(m_pTarget);
		return arr_obj.calcDesiredVel(actor);
	}
	// Step 2. if not first agent compute Valign as usual
	//TODO: add your code here
	double weight_sum = 0;
	for (int i = 0; i < agentList.size(); ++i) {
		BehaviorController* neighbor = agentList[i].getBehaviorController();
		if(neighbor == actor)
			continue;
		vec3 d = actorPos - neighbor->getPosition();
		if(d.Length() < BehaviorController::gKNeighborhood){
			Vdesired += neighbor->getVelocity();
			weight_sum += 1.0;
		}
	}	
	Vdesired *= BehaviorController::KAlignment;
	if (weight_sum > 0)
		Vdesired /= weight_sum;
	return Vdesired;
}

// Separation behavior
///////////////////////////////////////////////////////////////////////////////
// For the given te actor, return a desired velocity vector in world coordinates
// Separation tries to maintain a constant distance between all agents
// within the neighborhood
// agents[i] gives the pointer to the ith agent in the environment
// Separation settings are in BehavioralController::RNeighborhood and BehavioralController::KSeperate

 

Separation::Separation( AJoint* target,  vector<AActor>* agents) 
{
	m_name = "separation";
	m_AgentList = agents;
	m_pTarget = target;
}

Separation::~Separation()
{
}

Separation::Separation( Separation& orig) 
{
	m_name = "separation";
	m_AgentList = orig.m_AgentList;
	m_pTarget = orig.m_pTarget;
}

vec3 Separation::calcDesiredVel( BehaviorController* actor)
{
	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	vec3 targetPos = m_pTarget->getLocalTranslation();
	vec3 actorPos = actor->getPosition();
	vector<AActor>& agentList = *m_AgentList;
	
	// compute Vdesired = Vseparate
	// TODO: add your code here to compute Vdesired 
	double weight_sum = 0.0;
	for (int i = 0; i < agentList.size(); ++i) {
		BehaviorController* neighbor = agentList[i].getBehaviorController();
		if(neighbor == actor)
			continue;
		
		vec3 d = actorPos - neighbor->getPosition();
		if (d.Length() < BehaviorController::gKNeighborhood) {
			Vdesired += d / (d.Length()*d.Length());
			weight_sum++;
		}
	}
	Vdesired *= BehaviorController::KSeparation;
	/*
	if (weight_sum > 0)
		Vdesired /= weight_sum;
	*/
	if (Vdesired.Length() < 5.0)
		Vdesired = 0.0;
	
	return Vdesired;
}


// Cohesion behavior
///////////////////////////////////////////////////////////////////////////////
// For the given actor, return a desired velocity vector in world coordinates
// Cohesion moves actors towards the center of the group of agents in the neighborhood
//  agents[i] gives the pointer to the ith agent in the environment
//  Cohesion parameters are in BehavioralController::RNeighborhood and BehavioralController::KCohesion


Cohesion::Cohesion( vector<AActor>* agents) 
{
	m_name = "cohesion";
	m_AgentList = agents;
}

Cohesion::Cohesion( Cohesion& orig) 
{
	m_name = "cohesion";
	m_AgentList = orig.m_AgentList;
}

Cohesion::~Cohesion()
{
}

vec3 Cohesion::calcDesiredVel( BehaviorController* actor)
{
	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	vec3 actorPos = actor->getPosition();
	vector<AActor>& agentList = *m_AgentList;
	
	// compute Vdesired = Vcohesion
	// TODO: add your code here 
	vec3 Xcm(0, 0, 0);
	double weight_sum = 0.0;
	for (int i = 0; i < agentList.size(); ++i) {
		BehaviorController* neighbor = agentList[i].getBehaviorController();
		if (neighbor == actor)
			continue;
		vec3 neiPos = neighbor->getPosition();
		vec3 d = actorPos - neiPos;
		if(d.Length() < BehaviorController::gKNeighborhood){
			Xcm += neiPos;
			weight_sum += 1;
		}
	}
	if(weight_sum > 0)
		Xcm /= weight_sum;
	Vdesired = BehaviorController::KCohesion * (Xcm - actorPos);

	return Vdesired;
}

// Flocking behavior
///////////////////////////////////////////////////////////////////////////////
// For the given actor, return a desired velocity vector  in world coordinates
// Flocking combines separation, cohesion, and alignment behaviors
//  Utilize the Separation, Cohesion and Alignment behaviors to determine the desired velocity vector


Flocking::Flocking( AJoint* target,  vector<AActor>* agents) 
{
	m_name = "flocking";
	m_AgentList = agents;
	m_pTarget = target;
}

Flocking::Flocking( Flocking& orig) 
{
	m_name = "flocking";
	m_AgentList = orig.m_AgentList;
	m_pTarget = orig.m_pTarget;
}

Flocking::~Flocking()
{
}

vec3 Flocking::calcDesiredVel( BehaviorController* actor)
{
	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	vec3 actorPos = actor->getPosition();
	vector<AActor>& agentList = *m_AgentList;

	// compute Vdesired = Vflocking
	// TODO: add your code here 
	Separation sep_obj(m_pTarget, m_AgentList);
	vec3 Vseparate = sep_obj.calcDesiredVel(actor);
	Alignment ali_obj(m_pTarget, m_AgentList);
	vec3 Valign = ali_obj.calcDesiredVel(actor);
	Cohesion coh_obj(m_AgentList);
	vec3 Vcohesion = coh_obj.calcDesiredVel(actor);
	//double Cseparation = 0.3, Calignment = 0.4, Ccohesion = 0.3;      //define in behaviorControl and tune in GUI
	Vdesired = BehaviorController::Cseparation * Vseparate + BehaviorController::Calignment * Valign + BehaviorController::Ccohesion * Vcohesion;

	return Vdesired;
}

//	Leader behavior
///////////////////////////////////////////////////////////////////////////////
// For the given actor, return a desired velocity vector in world coordinates
// If the agent is the leader, move towards the target; otherwise, 
// follow the leader at a set distance behind the leader without getting to close together
//  Utilize Separation and Arrival behaviors to determine the desired velocity vector
//  You need to find the leader, who is always agents[0]

Leader::Leader( AJoint* target, vector<AActor>* agents) 
{
	m_name = "leader";
	m_AgentList = agents;
	m_pTarget = target;
}

Leader::Leader( Leader& orig) 
{
	m_name = "leader";
	m_AgentList = orig.m_AgentList;
	m_pTarget = orig.m_pTarget;
}

Leader::~Leader()
{
}

vec3 Leader::calcDesiredVel(BehaviorController* actor)
{
	
	vec3 Vdesired = vec3(0.0, 0.0, 0.0);
	vec3 actorPos = actor->getPosition();
	vector<AActor>& agentList = *m_AgentList;

	// TODO: compute Vdesired  = Vleader
	// followers should stay directly behind leader at a distance of -200 along the local z-axis

	//float CSeparation = 2.0;  float CArrival = 2.0;
	
	BehaviorController* leader = agentList[0].getBehaviorController(); // first agent is the leader
	mat3 Rmat = leader->getGuide().getLocalRotation();  // is rotattion matrix of lead agent

	//Vseparation
	Separation sep_obj(m_pTarget, m_AgentList);
	vec3 Vseparation = sep_obj.calcDesiredVel(actor);

	//Varrival
	vec3 targetPos = m_pTarget->getLocalTranslation();
	if (actor != leader)
		targetPos += Rmat * vec3(0, 0, -200);
	vec3 err = targetPos - actorPos;
	vec3 Varrival = BehaviorController::KArrival * err;

	if (actor == leader)
		return Varrival;
	else
		Vdesired = BehaviorController::Cseparation * Vseparation + BehaviorController::Carrival *Varrival;

	
	return Vdesired;
}

///////////////////////////////////////////////////////////////////////////////

