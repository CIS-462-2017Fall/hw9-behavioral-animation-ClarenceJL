#include "aSplineVec3.h"
#include <algorithm>
#include <Eigen/Dense>
#include <iostream>

#pragma warning(disable:4018)
#pragma warning(disable:4244)

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

ASplineVec3::ASplineVec3() : mInterpolator(new ALinearInterpolatorVec3())
{
}

ASplineVec3::~ASplineVec3()
{
    delete mInterpolator;
}

void ASplineVec3::setFramerate(double fps)
{
    mInterpolator->setFramerate(fps);
}

double ASplineVec3::getFramerate() const
{
    return mInterpolator->getFramerate();
}

void ASplineVec3::setLooping(bool loop)
{
    mLooping = loop;
}

bool ASplineVec3::getLooping() const
{
    return mLooping;
}

void ASplineVec3::setInterpolationType(ASplineVec3::InterpolationType type)
{
    double fps = getFramerate();

    delete mInterpolator;
    switch (type)
    {
    case LINEAR: mInterpolator = new ALinearInterpolatorVec3(); break;
    case CUBIC_BERNSTEIN: mInterpolator = new ABernsteinInterpolatorVec3(); break;
    case CUBIC_CASTELJAU: mInterpolator = new ACasteljauInterpolatorVec3(); break;
	case CUBIC_MATRIX: mInterpolator = new AMatrixInterpolatorVec3(); break; 
	case CUBIC_HERMITE: mInterpolator = new AHermiteInterpolatorVec3(); break;
	case CUBIC_BSPLINE: mInterpolator = new ABSplineInterpolatorVec3(); break;
    };
    
    mInterpolator->setFramerate(fps);
    computeControlPoints();
    cacheCurve();
}

ASplineVec3::InterpolationType ASplineVec3::getInterpolationType() const
{
    return mInterpolator->getType();
}

void ASplineVec3::editKey(int keyID, const vec3& value)
{
    assert(keyID >= 0 && keyID < mKeys.size());
    mKeys[keyID].second = value;
    computeControlPoints();
    cacheCurve();
}

void ASplineVec3::editControlPoint(int ID, const vec3& value)
{
    assert(ID >= 0 && ID < mCtrlPoints.size()+2);
    if (ID == 0)
    {
        mStartPoint = value;
        computeControlPoints();
    }
    else if (ID == mCtrlPoints.size() + 1)
    {
        mEndPoint = value;
        computeControlPoints();
    }
    else mCtrlPoints[ID-1] = value;
    cacheCurve();
}

void ASplineVec3::appendKey(double time, const vec3& value, bool updateCurve)
{
    mKeys.push_back(Key(time, value));

    if (mKeys.size() >= 2)
    {
        int totalPoints = mKeys.size();

        //If there are more than 1 interpolation point, set up the 2 end points to help determine the curve.
        //They lie on the tangent of the first and last interpolation points.
        vec3 tmp = mKeys[0].second - mKeys[1].second;
        double n = tmp.Length();
        mStartPoint = mKeys[0].second + (tmp / n) * n * 0.25; // distance to endpoint is 25% of distance between first 2 points

        tmp = mKeys[totalPoints - 1].second - mKeys[totalPoints - 2].second;
        n = tmp.Length();
        mEndPoint = mKeys[totalPoints - 1].second + (tmp / n) * n * 0.25;
    }

    if (updateCurve)
    {
        computeControlPoints();
        cacheCurve();
    }
}

void ASplineVec3::appendKey(const vec3& value, bool updateCurve)
{
    if (mKeys.size() == 0)
    {
        appendKey(0, value, updateCurve);
    }
    else
    {
        double lastT = mKeys[mKeys.size() - 1].first;
        appendKey(lastT + 1, value, updateCurve);
    }
}

void ASplineVec3::deleteKey(int keyID)
{
    assert(keyID >= 0 && keyID < mKeys.size());
    mKeys.erase(mKeys.begin() + keyID);
    computeControlPoints();
    cacheCurve();
}

vec3 ASplineVec3::getKey(int keyID)
{
    assert(keyID >= 0 && keyID < mKeys.size());
    return mKeys[keyID].second;
}

int ASplineVec3::getNumKeys() const
{
    return mKeys.size();
}

vec3 ASplineVec3::getControlPoint(int ID)
{
    assert(ID >= 0 && ID < mCtrlPoints.size()+2);
    if (ID == 0) return mStartPoint;
    else if (ID == mCtrlPoints.size() + 1) return mEndPoint;
    else return mCtrlPoints[ID-1];
}

int ASplineVec3::getNumControlPoints() const
{
    return mCtrlPoints.size()+2; // include endpoints
}

void ASplineVec3::clear()
{
    mKeys.clear();
}

double ASplineVec3::getDuration() const 
{
    return mKeys[mKeys.size()-1].first;
}

double ASplineVec3::getNormalizedTime(double t) const 
{
    return (t / getDuration());
}

vec3 ASplineVec3::getValue(double t)
{
    if (mCachedCurve.size() == 0) return vec3();

    double dt = mInterpolator->getDeltaTime();
    int rawi = (int)(t / dt); // assumes uniform spacing
    int i = rawi % mCachedCurve.size();
    double frac = t - rawi*dt;
    int inext = i + 1;
    if (!mLooping) inext = std::min<int>(inext, mCachedCurve.size() - 1);
    else inext = inext % mCachedCurve.size();

    vec3 v1 = mCachedCurve[i];
    vec3 v2 = mCachedCurve[inext];
    vec3 v = v1*(1 - frac) + v2 * frac;
    return v;
}

void ASplineVec3::cacheCurve()
{
    mInterpolator->interpolate(mKeys, mCtrlPoints, mCachedCurve);
}

void ASplineVec3::computeControlPoints()
{
    mInterpolator->computeControlPoints(mKeys, mCtrlPoints, mStartPoint, mEndPoint);
}

int ASplineVec3::getNumCurveSegments() const
{
    return mCachedCurve.size();
}

vec3 ASplineVec3::getCurvePoint(int i) const
{
    return mCachedCurve[i];
}

//---------------------------------------------------------------------
AInterpolatorVec3::AInterpolatorVec3(ASplineVec3::InterpolationType t) : mDt(1.0 / 120.0), mType(t)
{
}

void AInterpolatorVec3::setFramerate(double fps)
{
    mDt = 1.0 / fps;
}

double AInterpolatorVec3::getFramerate() const
{
    return 1.0 / mDt;
}

double AInterpolatorVec3::getDeltaTime() const
{
    return mDt;
}

void AInterpolatorVec3::interpolate(const std::vector<ASplineVec3::Key>& keys, 
    const std::vector<vec3>& ctrlPoints, std::vector<vec3>& curve)
{
	vec3 val = 0.0;
	double u = 0.0; 

	curve.clear();
	
	int numSegments = keys.size() - 1;
    for (int segment = 0; segment < numSegments; segment++)
    {
        for (double t = keys[segment].first; t < keys[segment+1].first - FLT_EPSILON; t += mDt)
        {
            
			// TODO: Compute u, fraction of duration between segment and segmentnext, for example,
            // u = 0.0 when t = keys[segment].first  
            // u = 1.0 when t = keys[segment+1].first
        	u = (t - keys[segment].first) / (keys[segment+1].first - keys[segment].first);

            val = interpolateSegment(keys, ctrlPoints, segment, u);
            curve.push_back(val);
        }
    }
	// add last point
	if (keys.size() > 1)
	{
		u = 1.0;
		val = interpolateSegment(keys, ctrlPoints, numSegments-1, u);
		curve.push_back(val);
	}
	
    
}

/**************** interpolator functions ********************/

vec3 ALinearInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
   
	vec3 curveValue(0, 0, 0);
	vec3 key0 = keys[segment].second;
    vec3 key1 = keys[segment+1].second;

    // TODO: 
	//Step 1: Create a Lerp helper function
	//Step 2: Linear interpolate between key0 and key1 so that u = 0 returns key0 and u = 1 returns key1
    curveValue = (1-u)*key0 + u*(key1);

	return curveValue;
}

vec3 ABernsteinInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 b0; 
	vec3 b1;
	vec3 b2; 
	vec3 b3;
    vec3 curveValue(0,0,0);
    // TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	b0 = ctrlPoints[4*segment];
	b1 = ctrlPoints[4*segment+1];
	b2 = ctrlPoints[4*segment+2];
	b3 = ctrlPoints[4*segment+3];
	// Step2: Compute the interpolated value f(u) point using Bernstein polynomials
    curveValue = (1-u)*(1-u)*(1-u)*b0 + 3*u*(1-u)*(1-u)*b1 + 3*u*u*(1-u)*b2 + u*u*u*b3;

	return curveValue;

}


vec3 ACasteljauInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);
	
	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	b0 = ctrlPoints[4*segment];
	b1 = ctrlPoints[4*segment+1];
	b2 = ctrlPoints[4*segment+2];
	b3 = ctrlPoints[4*segment+3];
	// Step2: Compute the interpolated value f(u) point using  deCsteljau alogithm
	vec3 b01, b11, b21, b02, b12;
	b01 = (1-u)*b0 + u*b1;
	b11 = (1-u)*b1 + u*b2;
	b21 = (1-u)*b2 + u*b3;
	b02 = (1-u)*b01 + u*b11;
	b12 = (1-u)*b11 + u*b21;
	curveValue = (1-u)*b02 + u*b12;

	return curveValue;
}

vec3 AMatrixInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);

	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	b0 = ctrlPoints[4*segment];
	b1 = ctrlPoints[4*segment+1];
	b2 = ctrlPoints[4*segment+2];
	b3 = ctrlPoints[4*segment+3];
	// Step2: Compute the interpolated value f(u) point using  matrix method f(u) = GMU
	// Hint: Using Eigen::MatrixXd data representations for a matrix operations
	MatrixXd G(3, 4), M(4, 4);
	VectorXd U(4), f(3);
	M << 1, -3, 3, -1,
		 0, 3, -6, 3, 
		 0, 0, 3, -3, 
		 0, 0, 0, 1;
	G << b0[0], b1[0], b2[0], b3[0],
		b0[1], b1[1], b2[1], b3[1],
		b0[2], b1[2], b2[2], b3[2];
	U << 1,
		 u,
		 u*u,
		 u*u*u;
	f = G*M*U;
	curveValue[0] = f(0); curveValue[1] = f(1); curveValue[2] = f(2);
	return curveValue;
}

vec3 AHermiteInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{

    vec3 p0 = keys[segment].second;
    vec3 p1 = keys[segment + 1].second;
    vec3 q0 = ctrlPoints[segment]; // slope at p0
    vec3 q1 = ctrlPoints[segment + 1]; // slope at p1
	vec3 curveValue(0, 0, 0);

    // TODO: Compute the interpolated value h(u) using a cubic Hermite polynomial  
	double H0, H1, H2, H3;

	H0 = 1 - 3 * u*u + 2 * u*u*u;
	H3 = 3 * u*u - 2 * u*u*u;
	H1 = u - 2 * u*u + u*u*u;
	H2 = -u*u + u*u*u; 
	curveValue = H0*p0 + H3*p1 + H1*q0 + H2*q1;

    return curveValue;
}



vec3 ABSplineInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 curveValue(0, 0, 0);

	// Hint: Create a recursive helper function N(knots,n,j,t) to calculate BSpline basis function values at t, where
	//     knots = knot array
	//	   n = degree of the spline curves (n = 3 for cubic)
	//     j = curve interval on knot vector in which to interpolate
	//     t = time value

	// Step 0: recover t from segment and u
	double t = u*(keys[segment + 1].first - keys[segment].first) + keys[segment].first;

	// Step 1: determine the index j
	int j = segment + 3;

	// Step 2: compute the n nonzero Bspline Basis functions N given j
	double N0 = N(3, j - 3, t), N1 = N( 3, j - 2, t), N2 = N( 3, j - 1, t), N3 = N(3, j, t);

	// Step 3: get the corresponding control points from the ctrlPoints vector
	vec3 c0 = ctrlPoints[j - 3], c1 = ctrlPoints[j - 2], c2 = ctrlPoints[j - 1], c3 = ctrlPoints[j];

	// Step 4: compute the Bspline curveValue at time t
	curveValue = c0*N0 + c1*N1 + c2*N2 + c3*N3 ;

	//matrix method
	/*
	MatrixXd G(3, 4), M(4, 4);
	VectorXd U(4), f(3);
	M << 1.0/6, -1.0/2, 1.0/2, -1.0/6,
		 2.0/3, 0, -1.0, 1.0/2,
		 1.0/6, 1.0/2, 1.0/2, -1.0/2,
		 0, 0, 0, 1.0/6;
	G << c0[0], c1[0], c2[0], c3[0],
		 c0[1], c1[1], c2[1], c3[1],
		 c0[2], c1[2], c2[2], c3[2];
	U << 1,
		 u,
		 u*u,
		 u*u*u;
	f = G*M*U;
	curveValue[0] = f(0); curveValue[1] = f(1); curveValue[2] = f(2);
	*/

	return curveValue;
}


/**************** compute control points functions ********************/

void ACubicInterpolatorVec3::computeControlPoints(
    const std::vector<ASplineVec3::Key>& keys, 
    std::vector<vec3>& ctrlPoints, 
    vec3& startPoint, vec3& endPoint)
{
    ctrlPoints.clear();
    if (keys.size() <= 1) return;

    for (int i = 1; i < keys.size(); i++)
    {
        vec3 b0, b1, b2, b3;

        // TODO: compute b0, b1, b2, b3

        //Question: startPoint, endPoint
        vec3 s0 = (i > 1) ? (keys[i].second - keys[i-2].second)/(keys[i].first - keys[i-2].first): (keys[i].second - startPoint)/2;
        vec3 s1 = (i < keys.size() - 1) ? (keys[i+1].second - keys[i-1].second)/(keys[i+1].first - keys[i-1].first): (endPoint - keys[i-1].second)/2;

        b0 = keys[i-1].second;
        b3 = keys[i].second;
        b1 = b0 + s0/3;
        b2 = b3 - s1/3;

        ctrlPoints.push_back(b0);
        ctrlPoints.push_back(b1);
        ctrlPoints.push_back(b2);
        ctrlPoints.push_back(b3);
    }
}

void AHermiteInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints,
	vec3& startPoint, vec3& endPoint)
{
	ctrlPoints.clear();
	if (keys.size() <= 1) return;

	int numKeys = keys.size();


	// TODO: 
	// For each key point pi, compute the corresonding value of the slope pi_prime.
	// Hints: Using Eigen::MatrixXd for a matrix data structures, 
	// this can be accomplished by solving the system of equations AC=D for C.
	// Don't forget to save the values computed for C in ctrlPoints
	// For clamped endpoint conditions, set 1st derivative at first and last points (p0 and pm) to s0 and s1, respectively
	// For natural endpoints, set 2nd derivative at first and last points (p0 and pm) equal to 0

	MatrixXd A = MatrixXd::Zero(numKeys, numKeys), D = MatrixXd::Zero(numKeys, 3), C(numKeys, 3);

	// Step 1: Initialize A
	A(0, 0) = 2; A(0, 1) = 1;
	A(numKeys - 1, numKeys - 2) = 1; A(numKeys - 1, numKeys - 1) = 2;
	for (int row = 1; row < numKeys - 1; ++row) {
		A(row, row - 1) = 1;
		A(row, row) = 4;
		A(row, row + 1) = 1;
	}

	// Step 2: Initialize D
	for (int i = 0; i < 3; ++i) {
		D(0, i) = 3 * (keys[1].second[i] - keys[0].second[i]);
		D(numKeys - 1, i) = 3 * (keys[numKeys - 1].second[i] - keys[numKeys - 2].second[i]);
		for (int row = 1; row < numKeys - 1; ++row)
			D(row, i) = 3 * (keys[row + 1].second[i] - keys[row - 1].second[i]);
	}

	// Step 3: Solve AC=D for C
	C = A.colPivHouseholderQr().solve(D);
	
	// Step 4: Save control points in ctrlPoints
	for (int row = 0; row < numKeys; ++row)
		ctrlPoints.emplace_back(vec3(C(row, 0), C(row, 1), C(row, 2)));

}



void ABSplineInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints,
	vec3& startPt, vec3& endPt)
{
	ctrlPoints.clear();
	if (keys.size() <= 1) return;

	// TODO: c
	// Hints: 
	// 1. use Eigen::MatrixXd to calculate the control points by solving the system of equations AC=D for C

	// 2. Create a recursive helper function dN(knots,n,t,l) to calculate derivative BSpline values at t, where
	//     knots = knot array
	//	   n = degree of the spline curves (n = 3 for cubic)
	//     j = interval on knot vector in which to interpolate
	//     t = time value
	//     l = derivative (l = 1 => 1st derivative)

	int numKeys = keys.size();
	// Step 1: Calculate knot vector using a uniform BSpline
	//         (assune knots are evenly spaced 1 apart and the start knot is at time = 0.0)
	mKnots.clear();
	mKnots = std::vector<double>(numKeys + 6);
	for (int i = 0; i < mKnots.size(); ++i) {
		if (i<3)
			mKnots[i] = keys[0].first - (3 - i);
		else if (i >= numKeys + 3)
			mKnots[i] = keys.back().first + (i - numKeys - 2);
		else
			mKnots[i] = keys[i - 3].first;
	}

	// Step 2: Calculate A matrix  for a natural BSpline
	//         (Set 2nd derivative at t0 and tm to zero, where tm is the last point knot; m = #segments)
	MatrixXd A = MatrixXd::Zero(numKeys + 2, numKeys + 2);
	A(0, 0) = dN(3, 0, keys[0].first, 2); A(0, 1) = dN(3, 1, keys[0].first, 2);
	A(0, 2) = dN( 3, 2, keys[0].first, 2); A(0, 3) = dN(3, 3, keys[0].first, 2);
	for (int r = 1; r < numKeys; ++r) {
		A(r, r - 1) = N(3, r - 1, keys[r - 1].first);
		A(r, r) = N(3, r, keys[r - 1].first);
		A(r, r+1) = N(3, r + 1, keys[r - 1].first);
		A(r, r+2) = N( 3, r + 2, keys[r - 1].first);
	}
	A(numKeys, numKeys - 2) = N(3, numKeys - 2, keys[numKeys - 1].first);
	A(numKeys, numKeys - 1) = N( 3, numKeys - 1, keys[numKeys - 1].first);
	A(numKeys, numKeys) = N( 3, numKeys, keys[numKeys - 1].first);
	A(numKeys, numKeys + 1) = N( 3, numKeys + 1, keys[numKeys - 1].first);

	A(numKeys + 1, numKeys - 2) = dN( 3, numKeys - 2, keys[numKeys - 1].first, 2);
	A(numKeys + 1, numKeys - 1) = dN( 3, numKeys - 1, keys[numKeys - 1].first, 2);
	A(numKeys + 1, numKeys) = dN( 3, numKeys, keys[numKeys - 1].first, 2);
	A(numKeys + 1, numKeys + 1) = dN( 3, numKeys + 1, keys[numKeys - 1].first, 2);

	// Step 3: Calculate  D matrix composed of our target points to interpolate
	MatrixXd D = MatrixXd::Zero(numKeys + 2, 3);
	for (int r = 1; r <= numKeys; ++r) {
		D(r, 0) = keys[r - 1].second[0];
		D(r, 1) = keys[r - 1].second[1];
		D(r, 2) = keys[r - 1].second[2];
	}

	// Step 4: Solve AC=D for C 
	MatrixXd C = A.colPivHouseholderQr().solve(D);

	// Step 5: save control points in ctrlPoints
	for (int r = 0; r < numKeys + 2; ++r)
		ctrlPoints.emplace_back(vec3(C(r, 0), C(r, 1), C(r, 2)));

}


/**************** B-spline helper functions ********************/

double ABSplineInterpolatorVec3::N(int n, int j, double t) {
	if (n == 0)
		return (t >= mKnots[j] && t < mKnots[j + 1]) ? 1 : 0;
	if (t < mKnots[j] || t > mKnots[j + n + 1])
		return 0;
	return ((t - mKnots[j]) / (mKnots[j + n] - mKnots[j]))*N( n - 1, j, t) +
		((mKnots[j + n + 1] - t) / (mKnots[j + n + 1] - mKnots[j + 1]))*N( n - 1, j + 1, t);
}

double ABSplineInterpolatorVec3::dN(int n, int j, double t, int l) {
	if (l == 0) return N(n, j, t);
	return n*(dN(n - 1, j, t, l - 1) / (mKnots[j + n] - mKnots[j]) - dN(n - 1, j + 1, t, l - 1) / (mKnots[j + n + 1] - mKnots[j + 1]));
}



/**************** Quintic interpolator functions ******************

vec3 AQuinticBernsteinInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u) {
	std::vector<vec3> b(ctrlPoints.begin()+4*segment,ctrlPoints.begin()+4*segment+6);
	vec3 curveValue(0, 0, 0);

	for (int i = 5; i >= 0; i--)
		curveValue += pow(1 - u, i)*pow(u, 5 - i)*b[5 - i];

	return curveValue;
}


void AQuinticInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints,
	vec3& startPoint, vec3& endPoint) {

	ctrlPoints.clear();
	if (keys.size() <= 1) return;

	for (int i = 1; i < keys.size(); i++)
	{
		vec3 b0, b1, b2, b3, b4, b5;

		vec3 s0 = (i > 1) ? (keys[i].second - keys[i - 2].second) / (keys[i].first - keys[i - 2].first) : (keys[i].second - startPoint) / 2;
		vec3 s1 = (i < keys.size() - 1) ? (keys[i + 1].second - keys[i - 1].second) / (keys[i + 1].first - keys[i - 1].first) : (endPoint - keys[i - 1].second) / 2;

		b0 = keys[i - 1].second;
		b3 = keys[i].second;
		b1 = b0 + s0 / 3;
		b2 = b3 - s1 / 3;
		// 2nd order continuity
		// 

		ctrlPoints.push_back(b0);
		ctrlPoints.push_back(b1);
		ctrlPoints.push_back(b2);
		ctrlPoints.push_back(b3);
		ctrlPoints.push_back(b4);
		ctrlPoints.push_back(b5);
	}
}
*/
