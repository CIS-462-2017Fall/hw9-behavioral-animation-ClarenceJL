#include "aIKController.h"
#include "GL/glut.h"

#include "aActor.h"

#include <iostream>      
#include <exception>

#include "aMatrix.h"

#pragma warning (disable : 4018)

int IKController::gIKmaxIterations = 5;
double IKController::gIKEpsilon = 0.1;

const double PI = 3.1415926535897;

// AIKchain class functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////
AIKchain::AIKchain()
{
	mWeight0 = 0.1;
}

AIKchain::~AIKchain()
{

}

AJoint* AIKchain::getJoint(int index) 
{ 
	return mChain[index]; 
}

void AIKchain::setJoint(int index, AJoint* pJoint) 
{ 
	mChain[index] = pJoint; 
}

double AIKchain::getWeight(int index) 
{ 
	return mWeights[index]; 
}

void AIKchain::setWeight(int index, double weight) 
{ 
	mWeights[index] = weight; 
}

int AIKchain::getSize() 
{ 
	return mChain.size(); 
}

std::vector<AJoint*>& AIKchain::getChain() 
{ 
	return mChain; 
}

std::vector<double>& AIKchain::getWeights() 
{ 
	return mWeights; 
}

void AIKchain::setChain(std::vector<AJoint*> chain) 
{
	mChain = chain; 
}

void AIKchain::setWeights(std::vector<double> weights) 
{ 
	mWeights = weights; 
}

// AIKController class functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////

IKController::IKController()
{
	m_pActor = NULL;
	m_pSkeleton = NULL;
	mvalidLimbIKchains = false;
	mvalidCCDIKchains = false;
	mvalidPseudoInvIKchains = false;
	mvalidOtherIKchains = false;

	// Limb IK
	m_pEndJoint = NULL;
	m_pMiddleJoint = NULL;
	m_pBaseJoint = NULL;
	m_rotationAxis = vec3(0.0, 1.0, 0.0);

	ATransform desiredTarget = ATransform();
	mTarget0.setLocal2Parent(desiredTarget);  // target associated with end joint
	mTarget1.setLocal2Parent(desiredTarget);  // optional target associated with middle joint - used to specify rotation of middle joint about end/base axis
	mTarget0.setLocal2Global(desiredTarget);
	mTarget1.setLocal2Global(desiredTarget);

	//CCD IK
	mWeight0 = 0.1;  // default joint rotation weight value

}

IKController::~IKController()
{
}

ASkeleton* IKController::getSkeleton()
{
	return m_pSkeleton;
}

const ASkeleton* IKController::getSkeleton() const
{
	return m_pSkeleton;
}

ASkeleton* IKController::getIKSkeleton()
{
	return &mIKSkeleton;
}

const ASkeleton* IKController::getIKSkeleton() const
{
	return &mIKSkeleton;
}

AActor* IKController::getActor()
{
	return m_pActor;
}

void IKController::setActor(AActor* actor)

{
	m_pActor = actor;
	m_pSkeleton = m_pActor->getSkeleton();
}


AIKchain IKController::createIKchain(int endJointID, int desiredChainSize, ASkeleton* pSkeleton)
{
	// TODO: given the end joint ID and the desired size (i.e. length) of the IK chain, 
	// 1. add the corresponding skeleton joint pointers to the AIKChain "chain" vector data member starting with the end joint
	// 2. also add weight values to the associated AIKChain "weights" vector data member for use in the CCD IK implemention
	// Note: desiredChainSize = -1 should create an IK chain of maximum length (i.e. where the last chain joint is the joint before the root joint)
	bool getMaxSize = false;

	int EndJointID = endJointID;
	std::vector<AJoint*> chain;
	std::vector<double> weights;

	chain.clear();
	weights.clear();
	if (desiredChainSize == -1)
		getMaxSize = true;

	if ((EndJointID >= 0) && (EndJointID < pSkeleton->getNumJoints()))
	{
		AJoint* pJoint = pSkeleton->getJointByID(endJointID);

		// TODO: add code here to generate chain of desired size or terminate at the joint before root joint, so that root will not change during IK	
		// also add weight values to corresponding weights vector  (default value = 0.1)
		if (getMaxSize) {
			while (pJoint != pSkeleton->getRootNode()) {
				chain.push_back(pJoint);
				weights.push_back(0.1); // you may need to change this later
				pJoint = pJoint->getParent();
			}
		}
		else {
			while (desiredChainSize--) {
				chain.push_back(pJoint);
				weights.push_back(mWeight0); // you may need to change this later
				pJoint = pJoint->getParent();
			}
		}
	}
	AIKchain result;
	result.setChain(chain);
	result.setWeights(weights);

	return result;
}

bool IKController::IKSolver_Limb(int endJointID, const ATarget& target)
{
	// Implements the analytic/geometric IK method assuming a three joint limb  

	if (!mvalidLimbIKchains)
	{
		mvalidLimbIKchains = createLimbIKchains();
		//assert(mvalidLimbIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;

	switch (endJointID)
	{
	case mLhandID:
		mLhandTarget = target;
		computeLimbIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
		break;
	case mRhandID:
		mRhandTarget = target;
		computeLimbIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
		break;
	case mLfootID:
		mLfootTarget = target;
		computeLimbIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
		break;
	case mRfootID:
		mRfootTarget = target;
		computeLimbIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
		break;
	case mRootID:
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeLimbIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
		computeLimbIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
		computeLimbIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
		computeLimbIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
		break;
	default:
		mIKchain = createIKchain(endJointID, 3, &mIKSkeleton);
		computeLimbIK(target, mIKchain, axisY, &mIKSkeleton);
		break;
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}


int IKController::createLimbIKchains()
{
	bool validChains = false;
	int desiredChainSize = 3;

	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);
	
	if (mLhandIKchain.getSize() == 3 && mRhandIKchain.getSize() == 3 && mLfootIKchain.getSize() == 3 && mRfootIKchain.getSize() == 3)
	{
		validChains = true;
		
		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}

int IKController::computeLimbIK(ATarget target, AIKchain& IKchain, const vec3 midJointAxis, ASkeleton* pIKSkeleton)
{
	// TODO: Implement the analytic/geometric IK method assuming a three joint limb  
	// The actual position of the end joint should match the target position within some episilon error 
	// the variable "midJointAxis" contains the rotation axis for the middle joint
	
	bool result = false;
	int endJointID;
	mTarget0 = target;

	if (IKchain.getSize() > 0)
		 endJointID = IKchain.getJoint(0)->getID();
	else endJointID = -1;

	if ((endJointID >= 0) && (endJointID < pIKSkeleton->getNumJoints()))
	{
		m_pEndJoint = IKchain.getJoint(0);
		m_pMiddleJoint = IKchain.getJoint(1);
		m_pBaseJoint = IKchain.getJoint(2);

		//TODO:
		// 1. compute error vector between target and end joint
		vec3 e = mTarget0.getGlobalTranslation() - m_pEndJoint->getGlobalTranslation();
		double le = e.Length();

		// 2. compute vector between end Joint and base joint
		vec3 r = m_pEndJoint->getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
		double lr = r.Length();

		// 3. compute vector between target and base joint
		vec3 rd = mTarget0.getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();

		// 4. Compute desired angle for middle joint 
		double lrd = rd.Length();
		double l1 = m_pMiddleJoint->getLocalTranslation().Length();
		double l2 = m_pEndJoint->getLocalTranslation().Length();
		//if the target lies beyond the configuration space:
		/*
		if (lrd < fabs(l1 - l2)) {
			rd = rd * fabs(l1 - l2) / lrd;
			lrd = fabs(l1 - l2);
		}
		if (lrd > l1 + l2) {
			rd = rd * (l1 + l2) / lrd;
			lrd = l1 + l2;
		}
		*/
		double cM = (l1*l1 + l2*l2 - lrd*lrd) / (2*l1*l2);
		if (cM < -1) cM = -1;
		if (cM > 1) cM = 1;
		double MiddleJointAngle = acos(cM);

		// 5. given desired angle and midJointAxis, compute new local middle joint rotation matrix and update joint transform
		mat3 MiddelJointRot; 
		MiddelJointRot.FromAxisAngle(midJointAxis, PI-MiddleJointAngle);
		m_pMiddleJoint->setLocalRotation(MiddelJointRot);
		m_pMiddleJoint->updateTransform();

		// 6. compute vector between target and base joint
		rd = mTarget0.getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
		lrd = rd.Length();

		// 7. Compute base joint rotation axis (in global coords) and desired angle
		vec3 baseJointAxis = r.Cross(rd); 
		baseJointAxis.Normalize();
		double cB = (lrd*lrd + lr*lr - le*le) / (2*lrd*lr);
		if (cB < -1)
			cB = -1;
		if (cB > 1)
			cB = 1;
		double baseJointAngle = acos(cB);

		// 8. transform base joint rotation axis to local coordinates
		baseJointAxis = m_pBaseJoint->getGlobalRotation().Transpose() * baseJointAxis;

		// 9. given desired angle and local rotation axis, compute new local rotation matrix and update base joint transform
		mat3 BaseJointRot;
		BaseJointRot.FromAxisAngle(baseJointAxis, baseJointAngle);
		BaseJointRot = m_pBaseJoint->getLocalRotation() * BaseJointRot;
		m_pBaseJoint->setLocalRotation(BaseJointRot);
		m_pBaseJoint->updateTransform();

	}

	result = true;
	return result;

}



bool IKController::IKSolver_CCD(int endJointID, const ATarget& target)
{
	// Implements the CCD IK method assuming a three joint limb 

	bool validChains = false;

	if (!mvalidCCDIKchains)
	{
		mvalidCCDIKchains = createCCDIKchains();
		//assert(mvalidCCDIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;

	switch (endJointID)
	{
	case mLhandID:
		mLhandTarget = target;
		computeCCDIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		break;
	case mRhandID:
		mRhandTarget = target;
		computeCCDIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		break;
	case mLfootID:
		mLfootTarget = target;
		computeCCDIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		break;
	case mRfootID:
		mRfootTarget = target;
		computeCCDIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
		break;
	case mRootID:
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeCCDIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		computeCCDIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		computeCCDIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		computeCCDIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
		break;
	default:
		mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
		computeCCDIK(target, mIKchain, &mIKSkeleton);
		break;
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}

int IKController::createCCDIKchains()
{
	bool validChains = false;

	int desiredChainSize = -1;  // default of -1 creates IK chain of maximum length from end joint to child joint of root


	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);

	if (mLhandIKchain.getSize() > 1 && mRhandIKchain.getSize() > 1 && mLfootIKchain.getSize() > 1 && mRfootIKchain.getSize() > 1)
	{
		validChains = true;

		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}

int IKController::computeCCDIK(ATarget target, AIKchain& IKchain, ASkeleton* pIKSkeleton)
{

	// TODO: Implement CCD IK  
	// The actual position of the end joint should match the desiredEndPos within some episilon error 

	bool result = false;

	mTarget0 = target;
	vec3 desiredEndPos = mTarget0.getGlobalTranslation();  // Get desired position of EndJoint

	int chainSize = IKchain.getSize();
	if (chainSize == 0) // There are no joints in the IK chain for manipulation
		return false;

	double epsilon = gIKEpsilon;
	int maxIterations = gIKmaxIterations;
	int numIterations = 0;

	m_pEndJoint = IKchain.getJoint(0);
	int endJointID = m_pEndJoint->getID();
	m_pBaseJoint = IKchain.getJoint(chainSize - 1);

	pIKSkeleton->copyTransforms(m_pSkeleton);

	if ((endJointID >= 0) && (endJointID < pIKSkeleton->getNumJoints()))
	{
		//TODO:
		for (; numIterations < maxIterations; ++numIterations) {
			vec3 e, r_i2end;
			AJoint * curr;
			for (int i = 1; i < chainSize; ++i) {
				curr = IKchain.getJoint(i);
				// 1. compute axis and angle for each joint in the IK chain (distal to proximal) in global coordinates
				e = desiredEndPos - m_pEndJoint->getGlobalTranslation();
				r_i2end = m_pEndJoint->getGlobalTranslation() - curr->getGlobalTranslation();
				double angle = (r_i2end ^ e).Length() / (r_i2end*r_i2end + e*e);
				vec3 axis = (r_i2end ^ e).Normalize();
				// 2. once you have the desired axis and angle, convert axis to local joint coords 
				axis = curr->getGlobalRotation().Transpose() * axis;
				// 3. multiply angle by corresponding joint weight value
				angle = IKchain.getWeight(i) * angle;
				// 4. compute new local joint rotation matrix
				mat3 rot;
				rot.FromAxisAngle(axis, angle);
				curr->setLocalRotation(curr->getLocalRotation() * rot);
				// 5. update joint transform
				curr->updateTransform();

			}
			//if the distance between end joint and target are smaller than threshold
			e = desiredEndPos - m_pEndJoint->getGlobalTranslation();
			if (e.Length() < epsilon)
				break;
		}// 6. repeat same operations above for each joint in the IKchain from end to base joint

	}
	result = true;
	return result;

}


/*
 * Solve Inverse Kinematics using Pseudo-Inverse of Jacobian 
 */

bool IKController::IKSolver_PseudoInv(int endJointID, const ATarget& target)
{
	bool result = false;

	// TODO: Implement Pseudo Inverse-based IK  
	// The actual position of the end joint should match the target position after the skeleton is updated with the new joint angles

	if (!mvalidPseudoInvIKchains){
		mvalidPseudoInvIKchains = createCCDIKchains();
		//assert(mvalidLimbIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;

	switch (endJointID)
	{
	case mLhandID:
		mLhandTarget = target;
		computePsudoInvIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		break;
	case mRhandID:
		mRhandTarget = target;
		computePsudoInvIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		break;
	case mLfootID:
		mLfootTarget = target;
		computePsudoInvIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		break;
	case mRfootID:
		mRfootTarget = target;
		computePsudoInvIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
		break;
	case mRootID:
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computePsudoInvIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		computePsudoInvIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		computePsudoInvIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		computePsudoInvIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
		break;
	default:
		mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
		computePsudoInvIK(target, mIKchain,  &mIKSkeleton);
		break;
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);


	return result;
}



int IKController::computePsudoInvIK(ATarget target, AIKchain& IKchain, ASkeleton* pIKSkeleton) {

	bool result = false;
	
	mTarget0 = target;
	vec3 desiredEndPos = mTarget0.getGlobalTranslation();  // Get desired position of EndJoint

	int chainSize = IKchain.getSize();
	if (chainSize == 0) // There are no joints in the IK chain for manipulation
		return false;

	double epsilon = 0.01;
	int maxIterations = gIKmaxIterations;
	int numIterations = 0;

	m_pEndJoint = IKchain.getJoint(0);
	int endJointID = m_pEndJoint->getID();
	m_pBaseJoint = IKchain.getJoint(chainSize - 1);

	
	AJoint* currJoint;

	matrix<double> Jacobian(3, 3 * (chainSize-1));
	matrix<double> Theta(3 * (chainSize - 1), 1);

	double maxReach = 0;
	
	for (int j = 1; j < chainSize; ++j) {
		maxReach += IKchain.getJoint(j - 1)->getLocalTranslation().Length();
		currJoint = IKchain.getJoint(j);
		//compute Bj
		mat3 Rj0 = currJoint->getGlobalRotation();
		vec3 rjn = m_pEndJoint->getGlobalTranslation() - currJoint->getGlobalTranslation();
		mat3 Bj;
		for (int col = 0; col < 3; ++col) {
			vec3 bc = Rj0.GetCol(col) ^ rjn;
			Bj.SetCol(col, bc);
		}
		//compute Lj
		mat3 Rj = currJoint->getLocalRotation();
		vec3 eulerAngle;
		Rj.ToEulerAngles(mat3::RotOrder::XYZ, eulerAngle);
		mat3 Lj;//default x-y-z order of rotation
		Lj[0][0] = cos(eulerAngle[1])*cos(eulerAngle[2]); Lj[0][1] = sin(eulerAngle[2]); Lj[0][2] = 0;
		Lj[1][0] = -cos(eulerAngle[1])*sin(eulerAngle[2]); Lj[1][1] = cos(eulerAngle[2]); Lj[1][2] = 0;
		Lj[2][0] = sin(eulerAngle[1]); Lj[2][1] = 0; Lj[2][2] = 1;
		//compute Jj
		mat3 Jj = Bj*Lj; 
		//copy to the full Jacobian matrix
		for (int r = 0; r < 3; ++r) {
			Theta(3 * (j - 1) + r, 0) = eulerAngle[r];
			for (int c = 0; c < 3; ++c)
				Jacobian(r, 3 * (j-1) + c) = Jj[r][c];
		}
	}

	vec3 displacement = desiredEndPos - m_pEndJoint->getGlobalTranslation();
	//handle the case when the target is out of reach:
	vec3 base2target = mTarget0.getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
	if (base2target.Length() > 0 && maxReach - epsilon < base2target.Length() ) {
		displacement -= (1 - (maxReach - epsilon ) / base2target.Length())*base2target;
	}
	
	//position shift
	matrix<double> dx(3, 1);
	dx(0, 0) = displacement[0]; dx(1, 0) = displacement[1]; dx(2, 0) = displacement[2];

	//right pseudo-inverse
	matrix<double> dTheta = (~Jacobian)*(!(Jacobian * (~Jacobian))) * dx;

	Theta += dTheta;
	
	for (int j = 1; j < chainSize; ++j) {
		vec3 eulerAngle(Theta(3*j-3,0), Theta(3*j-2,0), Theta(3*j-1,0));
		mat3 Rnew; 
		Rnew = Rnew.FromEulerAngles(mat3::RotOrder::XYZ, eulerAngle);
		currJoint = IKchain.getJoint(j);
		currJoint->setLocalRotation(Rnew);
		currJoint->updateTransform();
	}

	result = true;
	return result;
}

/*
 * An enhanced method combining Limb IK and CCD 
 */

bool IKController::IKSolver_Other(int endJointID, const ATarget& target)
{
	
	bool result = false;
	
	if (!mvalidOtherIKchains){
		mvalidOtherIKchains = createCCDIKchains();//create CCD chains (end joint to root) for this method 
		//assert(mvalidCCDIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;

	switch (endJointID){
	case mLhandID:
		mLhandTarget = target;
		computeOtherIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
		break;
	case mRhandID:
		mRhandTarget = target;
		computeOtherIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
		break;
	case mLfootID:
		mLfootTarget = target;
		computeOtherIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
		break;
	case mRfootID:
		mRfootTarget = target;
		computeOtherIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
		break;
	case mRootID:
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeOtherIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
		computeOtherIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
		computeOtherIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
		computeOtherIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
		break;
	default:
		mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
		computeOtherIK(target, mIKchain, axisY, &mIKSkeleton);
		break;
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	result = true;
	return result;
}


int IKController::computeOtherIK(ATarget target, AIKchain& IKchain, const vec3 midJointAxis, ASkeleton* pIKSkeleton) {
	bool result = false;
	int endJointID;
	mTarget0 = target;

	if (IKchain.getSize() > 0)
		endJointID = IKchain.getJoint(0)->getID();
	else endJointID = -1;

	if ((endJointID < 0) || (endJointID >= pIKSkeleton->getNumJoints()))
		return result;

	m_pEndJoint = IKchain.getJoint(0);
	m_pMiddleJoint = IKchain.getJoint(1);
	m_pBaseJoint = IKchain.getJoint(2);

	vec3 l1 = m_pEndJoint->getLocalTranslation(), l2 = m_pMiddleJoint->getLocalTranslation();
	vec3 dist = mTarget0.getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
	if (dist.Length() > l1.Length() + l2.Length()) 
		result = computeCCDIK(target, IKchain, pIKSkeleton);
	else {
		vector<AJoint*> shortChain({ m_pEndJoint, m_pMiddleJoint, m_pBaseJoint });
		vector<double> weights(3, 0.1);
		mvalidOtherIKchains = false;//re-compute IK-chain next time
		IKchain.setChain(shortChain);
		IKchain.setWeights(weights);
		result = computeLimbIK(target, IKchain, midJointAxis, pIKSkeleton);
	}
	result = true;
	return result;

};