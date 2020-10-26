#include "aActor.h"

#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor() 
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	 delete m_IKController;
	 delete m_BVHController;
	 delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();  

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode() )
		 return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	if (!m_pSkeleton->getRootNode()) { return; }

	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	// 2.	Set the y component of the guide position to 0
	// 3.	Set the global rotation of the guide joint towards the guideTarget
	m_Guide.setLocal2Global(m_Guide.getLocal2Global() * m_pSkeleton->getRootNode()->getLocal2Global());

	vec3 guidePos = m_Guide.getGlobalTranslation();
	guidePos[1] = 0.0;
	m_Guide.setGlobalTranslation(guidePos);

	mat3 rotMat;
	vec3 newDir = (guideTargetPos - m_Guide.getGlobalTranslation()).Normalize();
	rotMat.SetCol(0, vec3(0.0, 1.0, 0.0).Cross(newDir));
	rotMat.SetCol(1, vec3(0.0, 1.0, 0.0));
	rotMat.SetCol(2, newDir);
	m_Guide.setGlobalRotation(rotMat);
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space

	// 1.	Update the local translation of the root based on the left height and the right height
	AJoint* root = m_pSkeleton->getRootNode();
	vec3 rootPosLocal = root->getLocalTranslation();
	rootPosLocal[1] += std::max(leftHeight, rightHeight);
	root->setLocalTranslation(rootPosLocal);

	m_pSkeleton->update();
	
	// 2.	Update the character with Limb-based IK 
	
	// Rotate Foot
	if (rotateLeft)
	{
		// Update the local orientation of the left foot based on the left normal
		ATarget target;
		mat3 rotMat;
		rotMat.SetCol(0, leftNormal.Cross(m_Guide.getGlobalRotation()[2]));
		rotMat.SetCol(1, leftNormal);
		rotMat.SetCol(2, m_Guide.getGlobalRotation()[2]);
		m_IKController->IKSolver_Limb(leftFoot->getID(), target);
	}
	if (rotateRight)
	{
		// Update the local orientation of the right foot based on the right normal
		ATarget target;
		mat3 rotMat;
		rotMat.SetCol(0, rightNormal.Cross(m_Guide.getGlobalRotation()[2]));
		rotMat.SetCol(1, rightNormal);
		rotMat.SetCol(2, m_Guide.getGlobalRotation()[2]);
		m_IKController->IKSolver_Limb(rightFoot->getID(), target);
	}
	m_pSkeleton->update();
}
