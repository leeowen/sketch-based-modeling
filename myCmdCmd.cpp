//
// Copyright (C) Xiaosong Yang
// 
// File: ExtractCSsCmd.cpp
//
// MEL Command: ExtractCSs
//
// How to use it
// ExtractCSs curveShapeNodeName meshShapeNodeName NumberofCrossSections numberOfProfiles
// for example: 
//		ExtractCSs curveShape1 pCylinderShape1 4 10
//		ExtractCSs curveShape1 skin_patchShape 40 50
//
// Condition before use: 
//		1. because this program we use kObject space to get the coordinate, so before use, make sure to freeze the transformation of curve and mesh shape
//		2. make sure the cross section plane does intersect with the mesh in all 360 directions

/*
#include <maya/MSimple.h>
#include <maya/MSelectionList.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MBoundingBox.h>
#include <maya/MFnMesh.h>
#include <maya/MGlobal.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatVector.h>
#include <maya/MVector.h>
#include <maya/MQuaternion.h>
#include <math.h>
*/
//#include <C:\lyou\D_Drive\include\maya\MTypes.h>
//#include <C:\lyou\D_Drive\include\maya\MString.h>
#include <maya\MSimple.h>
#include <maya\MSelectionList.h>
#include <maya\MFnNurbsCurve.h>
#include <maya\MBoundingBox.h>
#include <maya\MFnMesh.h>
#include <maya\MGlobal.h>
#include <maya\MFloatPoint.h>
#include <maya\MFloatVector.h>
#include <maya\MVector.h>
#include <maya\MQuaternion.h>
#include <math.h>

DeclareSimpleCommand( ExtractCSs, "Xiaosong Yang", "7.0");

// second version, don't care some ray doesn't interest with the skin surface, just count how many intersection per-circle of ray
MStatus ExtractCSs::doIt( const MArgList& args )
//
//	Description:
//		implements the MEL ExtractCSs command.
//
//	Arguments:
//		args - the argument list that was passes to the command from MEL
//
//	Return Value:
//		MS::kSuccess - command succeeded
//		MS::kFailure - command failed (returning this value will cause the 
//                     MEL script that is being run to terminate unless the
//                     error is caught using a "catch" statement.
//
{
	MStatus stat = MS::kSuccess;

	// get the argument list: 1. curveNodeName; 2. Polygon Mesh Node name; 3. the number of CSs; 4. number of profiles
	int numArgs = args.length(&stat);
	if(numArgs != 4)
		return (MS::kFailure);

	MString curveNodeName, meshNodeName;
	int numCSs, numProfiles;
	stat = args.get(0, curveNodeName);
	stat = args.get(1, meshNodeName);
	stat = args.get(2, numCSs);
	stat = args.get(3, numProfiles);

	// find the obj of the curve
	MSelectionList tmpCurveSelect;
	tmpCurveSelect.add(curveNodeName);
	MObject curveNodeObj;
	stat = tmpCurveSelect.getDependNode(0, curveNodeObj);
	MFnNurbsCurve fnCurve(curveNodeObj);
	// find the obj of the mesh
	MSelectionList tmpMeshSelect;
	tmpMeshSelect.add(meshNodeName);
	MObject meshNodeObj;
	stat = tmpMeshSelect.getDependNode(0, meshNodeObj);
	MFnMesh fnMesh(meshNodeObj);
	MBoundingBox meshBox = fnMesh.boundingBox(&stat);
	double w = meshBox.width();
	double h = meshBox.height();
	double d = meshBox.depth();
	double biggestMeasure = sqrt(w*w + h*h + d*d);

	// settle down the curve local frame
	double curveLength = fnCurve.length();
	double step = curveLength*(1-0.00002)/numCSs; // leave curvelength*0.00001 at two ending
	double stepAngle = 3.1415926 *2.0/numProfiles;
	double currentAngle;
	double dist = curveLength*0.00001; // start from a little offset from Beginning point
	double para;
	double startPara =fnCurve.findParamFromLength(curveLength*0.00001, &stat);
	double endPara = fnCurve.findParamFromLength(curveLength*(1-0.00002), &stat);
	double stepPara = (endPara-startPara)/numCSs;
	MPoint currentPoint;
	MVector xAxis, yAxis, zAxis;
	MVector ray;
	MString buffer;
	buffer += numCSs;
	buffer += " ";
	buffer += numProfiles;
	buffer += "\n";
	MGlobal::displayInfo(buffer);
	buffer.clear();
	MFloatPoint hitPoint;
	bool succeed;
	MVector lastZAxis;
	para = startPara;
	xAxis = fnCurve.tangent(para, MSpace::kObject, &stat);
	xAxis.normalize();
	yAxis = fnCurve.normal(para, MSpace::kObject, &stat);
	stat= yAxis.normalize();
	lastZAxis = xAxis^yAxis;
	int numIntsect;

	for(int i=0;i<=numCSs;i++)
	{
		fnCurve.getPointAtParam(para, currentPoint, MSpace::kWorld);
		xAxis = fnCurve.tangent(para, MSpace::kObject, &stat);
		xAxis.normalize();
		yAxis = lastZAxis^xAxis;
		stat= yAxis.normalize();
		zAxis = xAxis^yAxis;
		lastZAxis = zAxis;
		currentAngle = 0.0;
		numIntsect = 0;
		for(int j=0;j<numProfiles;j++)
		{
			if(j==0)
				ray = yAxis;
			else
				ray = yAxis.rotateBy(MQuaternion(currentAngle, xAxis));
			// find the intersection point
			MFloatPoint tmpSource(currentPoint.x, currentPoint.y, currentPoint.z);
			MFloatVector tmpDirection(ray.x, ray.y, ray.z);
			succeed = fnMesh.closestIntersection(tmpSource, tmpDirection, NULL, NULL,false, MSpace::kObject, 100.0*biggestMeasure, false, NULL, hitPoint, NULL, NULL, NULL, NULL, NULL, 0.001, &stat);
			if(!succeed)
			{
			}
			else
			{
				buffer += MString(" ") + hitPoint.x + MString(" ") + hitPoint.y + MString(" ") + hitPoint.z;
				numIntsect++;
			}
			currentAngle += stepAngle;
		}
		buffer += "\n";
		MGlobal::displayInfo(MString(" ") + numIntsect + buffer);
		buffer.clear();

		para += stepPara;
	}

	// Since this class is derived off of MPxCommand, you can use the 
	// inherited methods to return values and set error messages
	//
	setResult( "ExtractCSs command executed!\n" );

	return stat;
}
