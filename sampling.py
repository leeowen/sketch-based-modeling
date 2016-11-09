import maya.cmds as cmds
	
for i in range(2):
	cmds.rename('curve'+str(i+1),'leg_curve'+str(i+1),ignoreShape=False)

sampleNum=20
legPoints=[] 
interval=1.0/sampleNum
for i in range(2):
	points=[]
	for j in range(sampleNum+1):
		pos=cmds.pointOnCurve('leg_curve'+str(i+1),turnOnPercentage=True,parameter=interval*j) 
		tuple(pos)
		points.append(pos)
	legPoints.append(points)

for i in range(2):
	curveAfterSample=cmds.curve(point=legPoints[i],degree=3)
	cmds.rename(curveAfterSample,'leg_curve_afterSampling'+str(i+1))
	

