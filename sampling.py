import maya.cmds as cmds
'''Name the body parts --> sample the contours --> create new contour curves based on the sampling points'''
	
def name_the_part(part):
	for i in range(2):
		cmds.rename('curve'+str(i+1),part+'_curve'+str(i+1),ignoreShape=False)


def sample_the_contour(sampleNum=20,part='leg',componentPoints):
	
	name_the_part(part)
	
	interval=1.0/sampleNum
	for i in range(2):
		points=[]
		for j in range(sampleNum+1):
			pos=cmds.pointOnCurve('leg_curve'+str(i+1),turnOnPercentage=True,parameter=interval*j) 
			tuple(pos)
			points.append(pos)
		componetPoints.append(points)

	#create new contour curves after the sampling
	for i in range(2):
		curveAfterSample=cmds.curve(point=componetPoints[i],degree=3)
		cmds.rename(curveAfterSample,part+'_curve_afterSampling'+str(i+1))
		
		
	

	
