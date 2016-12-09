from os import listdir
from os.path import isfile, join
from svgpathtools import Path, Line, QuadraticBezier, CubicBezier, Arc, svg2paths, wsvg, disvg
import math 
import cmath 
import numpy as np
#from emd import emd
from fitCurves import *
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d
from scipy.ndimage.filters import gaussian_filter
import time 
import pickle
dataPath = "./simpleData"
#dataPath = "./smallData"
queryPath = "querySVG.svg"
sampleLen = 10 
maxError = 5 
binNum = 10
def rOfCurve(curve, point):
	p = curve[0]*3*(-1)*(1-point)*(1-point) + \
	    curve[1]*3*((1-point)*(1-point) + point*2*(-1)*(1-point)) +\
	    curve[2]*3*((1-point)*2*point + point*point*(-1)) +\
	    curve[3]*3*point*point
	degree =  math.degrees(cmath.phase(p))
        if degree == 180:
                return 0
        if degree < 0:
        	if degree+ 180 == 180:
                	return 0 
                return degree + 180 
        if(degree >= 0)&(degree < 180):
		return degree 
def complexDist(c1, c2):
	return math.sqrt((c1.real-c2.real)**2 + (c1.imag-c2.imag)**2)
def lenOfCurve(curve, s, e):
	sample = int(math.floor(curve.length()/0.9 * (e-s)))
	length = 0
	if(sample < 1):
		return 0
	d = 1.0*(e - s)/sample
	for i in range(sample):
		length += complexDist(curve.point(s + i*d), curve.point(s+ (i+1)*d))
	return length
def radiusWithinRange(curve, start, rRange):
	if(1 - start < 0.001):	  
#		print "quick return"
#		print start
		return [start , 0]
	rNow = rOfCurve(curve, start) 
	rEnd = rOfCurve(curve, 1) 
	rMid = rOfCurve(curve, (1.0+start)/2) 
        if(rNow > rEnd):
		k = -1
	else:
		k = 1 
#	print "rNow = {}, rMid = {}, rEnd = {}".format(rNow,rMid,rEnd)
        if( not (k*rNow<= k*rMid<=k*rEnd)):
#		print "reverse"
		k = k * (-1)
#	print "rNow = {}".format(rNow)
	if(k == 1):
#		print "Up"
		rLevel =  math.floor(rNow/rRange)*rRange + rRange
	else:
#		print "Down"
		rLevel =  math.floor(rNow/rRange)*rRange
		if(rLevel == rNow):
			rLevel -= rRange
#	print "rLevel = {}".format(rLevel)
	if(rLevel > 180):
		rLevel -= 180
	if(rLevel < 0):
		rLevel += 180
#     179 => 2 UP
#     2 => 179 DOWN
	bStart = start
	bEnd = 1
	rEnd = rOfCurve(curve, 1)
	rMid = rOfCurve(curve, (start + 1)/2.0)
#	if((rEnd*k < rLevel*k) & (rMid*k < rLevel*k)):
#		return [1, lenOfCurve(curve,start,1)]
	r = 0
#	print "rLevel = {}".format(rLevel)
#        while((bEnd - bStart)> 0.005):
        while((bEnd - bStart)> 0.00005):
		bNow = (bStart + bEnd)/2.0
#		print "bEnd {}, r {}".format(bEnd, r)
		r = rOfCurve(curve, bNow)
		if( k == 1):
			if(rNow<= r <=rEnd):
				if(r > rLevel):
					bEnd = bNow
				else:
					bStart = bNow	
			else:		
				if(r < rEnd):
					bEnd = bNow
				else:
					if(r == 0):
						r = 180.0
					if(r > rLevel):
						bEnd = bNow
					else:
						bStart = bNow	
		else:
			if(rEnd<= r <=rNow):
				if(r < rLevel):
					bEnd = bNow
				else:
					bStart = bNow	
			else:		
				if(r > rNow):
					bEnd = bNow
				else:
					if(r < rLevel):
						bEnd = bNow
					else:
						bStart = bNow	
#	if(lenOfCurve(curve,start, bEnd)>100):
#		print "error??"
#		raw_input();
#	print "length is {}".format(lenOfCurve(curve,start, bEnd))
	return bEnd, lenOfCurve(curve,start, bEnd)
#Cut the line into several line segments
# input : a path 
# output : a path
def lineSegment(path, rRange=3.0):
        newPath = Path()
        for segment in path:
                if type(segment) is  Line :
                        newPath.append(segment)
                if type(segment) is  CubicBezier :
#                        n = int(math.floor(segment.length()/sampleLen)) + 1
#                        for i in range(n) :
#                                newPath.append(Line(segment.point(1.0*i/n), segment.point(1.0*(i+1)/n)))
#			print "Cut Bezier Curve"
			startPoint = segment[0]
			endPoint = segment[0]
			start = 0
			end = 0
#			while(cmath.polar(endPoint - segment.point(1))[0] > 0.001):
			while(start <  0.999):	  
				end, length = radiusWithinRange(segment, start, rRange)
				endPoint = segment.point(end)
                                l = cmath.polar(endPoint - startPoint)[0]
				if(l != 0):	
#					print "l != 0"
#                           		print cmath.polar(endPoint - startPoint)
#                                	print Line(startPoint, startPoint + length*(startPoint - endPoint)/cmath.polar(endPoint - startPoint)[0])
                        	        newPath.append(Line(startPoint, startPoint + length*(startPoint - endPoint)/cmath.polar(endPoint - startPoint)[0]))
					startPoint = endPoint
					start = end
        return newPath  

#Fit the line by Curve 
# input : a path with bezier curve
# output : a path only with line
def pathFitCurve(path, maxError, smooth = False):
	if(smooth == False):
		return path
	path2 = Path()
        for segment in path:
                if segment.length() > 0: 
       			path2.append(segment)
        for segment in path2:
                if type(segment) is  CubicBezier :
        		return path2  
        if len(path2) == 0 :
		return path2  
	points = [np.array([path2[0].start.real, path2[0].start.imag])]
	for segment in path2:
		points.append(np.array([segment.end.real, segment.end.imag]))
	curves = fitCurve(points, maxError)
	newPath = Path()
	for curve in curves:
		if len(curve) > 2:
			newPath.append(CubicBezier(complex(curve[0][0], curve[0][1]),
					  	   complex(curve[1][0], curve[1][1]),
						   complex(curve[2][0], curve[2][1]),
						   complex(curve[3][0], curve[3][1])))
		if len(curve) == 2:
			newPath.append(Line(complex(curve[0][0], curve[0][1]),
					    complex(curve[1][0], curve[1][1])))
	return newPath

# preprocess all the path in the list
# input : a list of path in the image
# output : only a list of line segment
def preprocess(pathList, smooth = False, rRange=3):
        newPath = Path()
        for path in pathList:
		path = pathFitCurve(path, maxError, smooth)
                for segment in path:
                        newPath.append(segment)
#	print "rRange = {}".format(rRange)
        newPath = lineSegment(newPath, rRange)
	return newPath
# calculate the angle of a line
# input : a line
# output : angle
def calcRadians(line):
        diff = line.start-line.end
        degree =  math.degrees(cmath.phase(diff))
        if degree == 180:
                return 0
        if degree < 0:
        	if degree+ 180 == 180:
                	return 0 
                return degree + 180 
        if(degree >= 0)&(degree < 180):
                return degree

def calcRadiansFeature(path):
	radiansFeature = []
	for line in path:
		radiansFeature.append(calcRadians(line))
	return np.array(radiansFeature)

def calcWeightFeature(path):
        if len(path) == 0: 
		return np.array([])
	lineWeight = [ segment.length() for segment in path]
	totalWeight = sum(lineWeight)
        if totalWeight == 0: 
		print path
		return np.array([])
		
	weightFeature = [ w / totalWeight for w in lineWeight]
	return np.array(weightFeature)

###WHAT is THE BUGGGG???
def calcMidFeature(path):
	'''
	lineWeight = [ segment.length() for segment in path]
	totalWeight = sum(lineWeight)
        print totalWeight
	midFeature = [ (segment.start + segment.end)/2 for segment in path]
        pointSum = complex(0,0)
	for i in range(len(path)):
		pointSum += lineWeight[i] * midFeature[i]
        
        if(totalWeight == 0):
		return np.array(midFeature)

	pointCenter = pointSum/totalWeight
	midFeature = [ p - pointCenter for p in midFeature]
	'''
	if len(path) == 0:
		return np.array([])	
	x = [ segment.start.real for segment in path] + [ segment.end.real for segment in path]
	max_x = max(x)
	min_x = min(x)
	y = [ segment.start.imag for segment in path] + [ segment.end.imag for segment in path]
	max_y = max(y)
	min_y = min(y)
	if(max_x == min_x):
		return []
	if(max_y == min_y):
		return []
	midFeature = []
        dx = max_x-min_x
	dy = max_y-min_y
	r = math.sqrt(dx*dx + dy*dy)
	for segment in path:
		mid = (segment.start + segment.end)/2
		mid = complex((mid.real - min_x)*100/r, (mid.imag - min_y)*100/r)
		midFeature.append(mid)
	return [np.array(midFeature), r]

def calcFeature(rFeature, wFeature, binNum):
	feature = [0] * binNum 
	for i in range(len(rFeature)):
		r = rFeature[i]
		w = wFeature[i]
		feature[int(math.floor(r*binNum/180.0))] += w
	return np.array(feature)

def calcFeature2(rFeature, wFeature, mFeature, binNum):
	feature = []
	for i in range(len(rFeature)):
		r = rFeature[i]
		w = wFeature[i]
		p = mFeature[0][i]
		rad = mFeature[1]
		dx = p.real
		dy = p.imag
		dl1 = dx * math.sin(r*pi/180.0) - dy * math.cos(r**pi/180.0)
		dl2 = -dx * math.cos(r*pi/180.0) - dy * math.sin(r**pi/180.0)
#		feature[int(math.floor(r*binNum/180.0))] += w
		feature.append([ r, dl1/rad, dl2/rad, w])
	return np.array(feature)

def calcVariance(f2, binNum):
	feature =[]
	tmp = [] 
	for i in range(binNum):
		tmp.append([])	
	for i in range(len(f2)):
		r = f2[i][0]
		x = f2[i][1]
		y = f2[i][2]
		w = f2[i][3]
		tmp[int(math.floor(r*binNum/180.0))].append([x, y, w])
	for i in range(binNum):
		angleList = tmp[i]
		sum_x = 0
		sum_y = 0
		mean_x = 0
		mean_y = 0

		for j in range(len(angleList)):
			sum_x += angleList[j][1] 
			sum_y += angleList[j][2] 

		mean_x = sum_x*1.0 / len(tmp[i]) 
		mean_y = sum_y*1.0 / len(tmp[i])
		
		rfeature = [0,0]
		for j in range(len(angleList)):
			rfeature[0] += angleList[j][3]*((angleList[j][1] - mean_x)**2)
			rfeature[1] += angleList[j][3]*((angleList[j][2] - mean_y)**2)
		feature.append(rfeature)
	return feature
def calcFeature3(f2, binNum):
	#50 * 50
	size = 5
	featureTmp = [[]] * binNum
	for i in range(len(f2)):
		featureTmp[int(math.floor(f2[i][0]*binNum/180.0))].append(f2[i]) 
	feature = [[]] * binNum
	for i in range(binNum):
		feature[i] = np.zeros([size,size])
		for p in featureTmp[i]:
			feature[i][int(math.floor(p[1]*4.9))][int(math.floor(p[2]*4.9))] += p[3]
	return np.array(feature)
	
def calcDistance0(f1, f2):
	dist = 0
	for i in range(len(f1)):
		diff = f1[i] - f2[i]
		dist += diff * diff
	return dist

def calcDistance(f1, f2, k):
	dist = 0
	l = len(f1)
	for i in range(l):
		diff = f1[(i+k)%l] - f2[i]
		dist += diff * diff
	return dist


def calcDistance(f1, f2, k, m):
	a = 0.5
	b = 1-a
	dist = 0
	v1 = f1[1]
	v2 = f2[1]
	f1 = f1[m]
	f2 = f2[m]
	l = len(f1)
	for i in range(l):
		diff = f1[(i+k)%l] - f2[i]
		diff2 = v1[(i+k)%l][0] - v2[i][0]
		diff3 = v1[(i+k)%l][1] - v2[i][1]
		
		dist += a*diff * diff
		dist += b*diff2 * diff2
		dist += b*diff3 * diff3
	return dist


def calcDistance3(f1, f2):
	dist = 0
	for i in range(len(f1)):
		diff = f1[i] - f2[i]
		if diff < 0:
			diff =  -diff
		diff = diff * diff * diff * 1000 
		dist = dist + diff
	return dist
def calcDistance5(f1, f2):
	dist = 0
	for i in range(len(f1)):
		diff = min(f1[i], f2[i])
		dist = dist + diff
	return dist
	

def calcDistanceCD(f1, f2, k):
	dist = 0
	l = len(f1)
	for i in range(l):
		diff = f1[(i+k)%l] * f2[i]
		dist = dist + diff
	return 1-dist
#JSD
def calcDistanceJSD(f1, f2, k):
	dist = 0
	l = len(f1)
	for i in range(l):
		dsum = f1[(i+k)%l] + f2[i]
		diff1 = 0
		diff2 = 0
		if(f1[(i+k)%l] != 0):
		 diff1 = f1[(i+k)%l] * math.log( 2*f1[(i+k)%l] / dsum ) 
		if(f2[i] != 0):
		 diff2 = f2[i] * math.log( 2*f2[i] / dsum ) 
		dist = dist + diff1 + diff2
	return dist

def calcAllFeature(rAllFeature, wAllFeature, binNum):
	allFeature = []
	for i in range(len(rAllFeature)):
		allFeature.append(calcFeature(rAllFeature[i], wAllFeature[i], binNum))
	return np.array(allFeature)

	

def find(allFeature, queryFeature):
	target = -1
	targetDist = 2
	for i in range(len(allFeature)):
		dist = calcDistance(allFeature[i], queryFeature)
		if dist < targetDist:
			target = i
			targetDist = dist
	return target

def parse(path):
	queryData = preprocess(svg2paths(queryPath)[0])
	queryRadiansFeature = calcRadiansFeature(queryData) 
	queryWeightFeature = calcWeightFeature(queryData) 
	queryFeature = calcFeature(queryRadiansFeature, queryWeightFeature)
	return queryFeature
def testPath(path, smooth=False, binNum=60):
	path = preprocess([path], smooth)
	outputData = []
	p = path
	r = calcRadiansFeature(p)
	w = calcWeightFeature(p)
        m = calcMidFeature(p)
	f = calcFeature(r, w, binNum)
#	f2 = calcFeature2(r, w, m, binNum)
#	f3 = calcFeature3(f2, binNum)
#	outputData.append([r, w, f, f2, f3, gaussian_filter1d(f, 0.2, truncate=3), 
	f4 = [ mm-1/binNum for mm in gaussian_filter1d(f, 0.6, truncate=3)]
	outputData = [r, w, f, f4, gaussian_filter1d(f, 0.2, truncate=3), 
	                           gaussian_filter1d(f, 0.4, truncate=3), 
	                           gaussian_filter1d(f, 0.5, truncate=3), 
	                           gaussian_filter1d(f, 0.6, truncate=3), 
	                           gaussian_filter1d(f, 0.8, truncate=3), p, path]
	return outputData
	
 
def readSVG(svgPath, smooth=False, binNum=60):
	start = time.time()
	start2 = time.time()
	paths = svg2paths(svgPath)[0]
	end2 = time.time()
	print " read file cost {} s".format(end2 -start2)
	start2 = time.time()
	path = preprocess(paths, smooth, 180.0/binNum) 
	end2 = time.time()
	print " preprocess cost {} s".format(end2 -start2)
	outputData = []
	p = path
	start2 = time.time()
	r = calcRadiansFeature(p)
	w = calcWeightFeature(p)
        m = calcMidFeature(p)
	f = calcFeature(r, w, binNum)
#	f2 = calcFeature2(r, w, m, binNum)
#	f3 = calcFeature3(f2, binNum)
#	outputData.append([r, w, f, f2, f3, gaussian_filter1d(f, 0.2, truncate=3), 
	f4 = [ mm-1/binNum for mm in gaussian_filter1d(f, 0.6, truncate=3)]
	outputData = [r, w, f, f4, gaussian_filter1d(f, 0.2, truncate=3), 
	                           gaussian_filter1d(f, 0.4, truncate=3), 
	                           gaussian_filter1d(f, 0.5, truncate=3), 
	                           gaussian_filter1d(f, 0.6, truncate=3), 
	                           gaussian_filter1d(f, 0.8, truncate=3), p, svgPath]
	end2 = time.time()
	print " calculate feature cost {} s".format(end2 -start2)
	end = time.time()
	print "cost {} s".format(end -start)
	return outputData, paths, path

def readDataSet0(path, smooth=False, binNum=60):
	start = time.time()
	files = [f for f in listdir(dataPath) if isfile(join(dataPath, f))]
	if ".DS_Store" in files :
		files.remove(".DS_Store") 

	#read svg from file
	svgData0 = [] 
	for p in files:
		path = svg2paths(join(dataPath, p))[0]
		svgData0.append(path)
	svgData = []
	for i in range(len(svgData0)):
		start2 = time.time()
		print i
		svgData.append(preprocess(svgData0[i], smooth, 180.0/binNum)) 
		end2 = time.time()
		print "cost {} s".format(end2 -start2)
	outputData = []
	for i in range(len(svgData)):
		p = svgData[i]
		r = calcRadiansFeature(p)
		w = calcWeightFeature(p)
		f = calcFeature(r, w, binNum)
                m = calcMidFeature(p)
		f2 = calcFeature2(r, w, m, binNum)
		v = calcVariance(f2)
#		f3 = calcFeature3(f2, binNum)
#		outputData.append([r, w, f, f2, f3, gaussian_filter1d(f, 0.2, truncate=3), 
		outputData.append([f, v, files[i]])
	end = time.time()
	print "cost {} s".format(end -start)
	return outputData
def changeToData(data):
	outputData = []
	for i in range(len(data)):
		f = data[i][0]
		fileName = data[i][-1]
		m = data[i][1]
		f2 = data[i][2]
		v = calcVariance(f2)
		outputData.append([f2, v, f, gaussian_filter1d(f, 0.2, truncate=3), 
		                                    gaussian_filter1d(f, 0.4, truncate=3), 
		                                    gaussian_filter1d(f, 0.5, truncate=3), 
		                                    gaussian_filter1d(f, 0.6, truncate=3), 
		                                    gaussian_filter1d(f, 0.7, truncate=3), 
		                                    gaussian_filter1d(f, 0.8, truncate=3), 
		                                    gaussian_filter1d(f, 0.9, truncate=3), 
		                                    gaussian_filter1d(f, 0.95, truncate=3), fileName])
	return outputData
		
	
def readDataSet(path, smooth=False, binNum=60):
	start = time.time()
	files = [f for f in listdir(dataPath) if isfile(join(dataPath, f))]
	if ".DS_Store" in files :
		files.remove(".DS_Store") 

	#read svg from file
	svgData0 = [] 
	for p in files:
		path = svg2paths(join(dataPath, p))[0]
		svgData0.append(path)
	svgData = []
	for i in range(len(svgData0)):
		start2 = time.time()
		print i
		svgData.append(preprocess(svgData0[i], smooth, 180.0/binNum)) 
		end2 = time.time()
		print "cost {} s".format(end2 -start2)
	outputData = []
	for i in range(len(svgData)):
		p = svgData[i]
		r = calcRadiansFeature(p)
		w = calcWeightFeature(p)
                m = calcMidFeature(p)
		f = calcFeature(r, w, binNum)
#		f2 = calcFeature2(r, w, m, binNum)
#		f3 = calcFeature3(f2, binNum)
#		outputData.append([r, w, f, f2, f3, gaussian_filter1d(f, 0.2, truncate=3), 
		f4 = [ mm-1/binNum for mm in gaussian_filter1d(f, 0.6, truncate=3)]
		outputData.append([r, w, f, f4, gaussian_filter1d(f, 0.2, truncate=3), 
		                                    gaussian_filter1d(f, 0.4, truncate=3), 
		                                    gaussian_filter1d(f, 0.5, truncate=3), 
		                                    gaussian_filter1d(f, 0.6, truncate=3), 
		                                    gaussian_filter1d(f, 0.7, truncate=3), 
		                                    gaussian_filter1d(f, 0.8, truncate=3), 
		                                    gaussian_filter1d(f, 0.9, truncate=3), 
		                                    gaussian_filter1d(f, 0.95, truncate=3), svgData0[i], p, files[i]])
	end = time.time()
	print "cost {} s".format(end -start)
	return outputData

def distF(f1, f2, k=0):
	target = 0 
	targetDist = 10000000 
	for i in range(binNum):
		dist = calcDistance(f1, f2, i, k)
		if dist < targetDist:
			target = i
			targetDist = dist
	return [targetDist, target]

def distF2(f1, f2, k):
	f1 = f1[k]
	f2 = f2[k]
	target = 0 
	targetDist = 10000000 
	for i in range(binNum):
		dist = calcDistance(f1, f2, i)
		if dist < targetDist:
			target = i
			targetDist = dist
	return [targetDist, target]


def find(allF, queryF, k=2, num=10):
	start = time.time()
	resultData = []
	num = min([num, len(allF)])
	for i in range(len(allF)):
		distnow = distF(allF[i], queryF, k)
		distnow.append(allF[i][-1])
		distnow.append(i)
		resultData.append(distnow)
	result = np.argsort(np.array(resultData), axis = 0)
	answer = []
	for i in range(num):
		answer.append(resultData[result[i][0]])
	end = time.time()
	print "cost {} s".format(end -start)
	return answer
def graph(g,binNum=60):
	n, bins, patches = plt.hist(g[0], binNum, range=(0,180), weights=g[1])
	plt.show()
def tagOfPath(path):
	return path[:path.index("-")]
def evaluateQ(allF, queryF, k=2, num=20):
	answer = find(allF, queryF, k, num)
		
	target = tagOfPath(queryF[-1])
	
	right = 0
	for i in range(len(answer)):
		if(target == tagOfPath(answer[i][-2])):
			right += 1
#		print "round {}, recall = {}, precision = {}".format(i, 1.0*right/20, 1.0*right/(i+1))
	print "round {}, recall = {}, precision = {}".format(i, 1.0*right/20, 1.0*right/(i+1))
        return [1.0*right/20, 1.0*right/(i+1)]
def evaluate(allF, k = 2, num = 20):
	start = time.time()
	rList = [];
	pList = [];
	for i in range(len(allF)):	
		print "test {} {}".format(i, allF[i][-1])
		r, p = evaluateQ(allF, allF[i], k, num)
		rList.append(r)
		pList.append(p)	
		end = time.time()
		print "cost {} s".format(end -start)
		print "At {}, recall = {}, precision = {}".format(i, sum(rList)/(i+1), sum(pList)/(i+1))
	end = time.time()
	print "cost {} s".format(end -start)
	print "At {}, recall = {}, precision = {}".format(num, sum(rList)/len(allF), sum(pList)/len(allF))
	return [rList,pList]
def saveTo(data, path):
	pickle.dump(data, open(path, "wb"), True)
def loadTo(path):
	return pickle.load(open(path, "rb"))
def seeRecall(path, path2):
	r = pickle.load(open(path, "rb"))
	data = pickle.load(open(path2, "rb"))
	cdf = [0]*len(r)
	for i in range(len(r)):
		if i == 0:
			cdf[i] = r[i]
		else:
			cdf[i] = (cdf[i-1] + r[i])
	for i in range(len(r)):
		cdf[i] /= i+1
	d = {}
	for i in range(len(data)):
		d[tagOfPath(data[i][-1])] = 0
	for i in range(len(data)):
		d[tagOfPath(data[i][-1])] += r[i]/20.0
	
	return [d, cdf]	
def barChart0():
	n_groups = 7 

	recall_ED = [0.3215, 0.3156, 0.299, 0.2710, 0.2856, 0.272, 0.279]
	recall_JD = [0.35928, 0.381, 0.313, 0.2821, 0.2912, 0.303, 0.282]

	fig, ax = plt.subplots()

	index = np.arange(n_groups)
	bar_width = 0.35

	opacity = 0.4
	error_config = {'ecolor': '0.3'}

	rects1 = plt.bar(index, recall_ED, bar_width,
                 alpha=opacity,
                 color='b',
                 error_kw=error_config,
                 label='ED')

	rects2 = plt.bar(index + bar_width, recall_JD, bar_width,
                 alpha=opacity,
                 color='r',
                 error_kw=error_config,
                 label='JD')

	plt.xlabel('Number of bins')
	plt.ylabel('Recall rate')
	plt.title('Recall rate by Number of bins and Distance measurement')
	plt.xticks(index + bar_width, ('10', '30', '60', '90', '120', '150', '180'))
	plt.legend()

	plt.tight_layout()
	plt.show()

def barChart():
	n_groups = 4 

	recall_NO = [0.35928, 0.381, 0.2821, 0.282]
	recall_02 = [0.3635, 0.3823, 0.2746, 0.3321]
	recall_06 = [0.3513, 0.3833, 0.3081, 0.2953]
	recall_08 = [0.3421, 0.3874, 0.3707, 0.3785]

	fig, ax = plt.subplots()

	index = np.arange(n_groups)
	bar_width = 0.35

	opacity = 0.4
	error_config = {'ecolor': '0.3'}

	rects1 = plt.bar(index*2, recall_NO, bar_width,
                 alpha=opacity,
                 color='b',
                 error_kw=error_config,
                 label='No filter')

	rects2 = plt.bar(index*2 + bar_width, recall_02, bar_width,
                 alpha=opacity,
                 color='r',
                 error_kw=error_config,
                 label='Filter, standard deviation  = 0.2')
	rects3 = plt.bar(index*2 + bar_width*2, recall_06, bar_width,
                 alpha=opacity,
                 color='g',
                 error_kw=error_config,
                 label='Filter, standard deviation  = 0.6')
	rects4 = plt.bar(index*2 + bar_width*3, recall_08, bar_width,
                 alpha=opacity,
                 color='c',
                 error_kw=error_config,
                 label='Filter, standard deviation  = 0.8')

	plt.xlabel('Number of bins')
	plt.ylabel('Recall rate')
	plt.ylim([0, 0.6])
	plt.title('Recall rate by Number of bins and filter setting')
	plt.xticks(index*2 + bar_width*2, ('10', '30', '90', '180'))
	plt.legend()

	plt.tight_layout()
	plt.show()




data = readDataSet(dataPath, False, binNum)
