from os import listdir
from os.path import isfile, join
from svgpathtools import Path, Line, QuadraticBezier, CubicBezier, Arc, svg2paths, wsvg, disvg
import math 
import cmath 
import numpy as np
from emd import emd
from fitCurves import *
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d
from scipy.ndimage.filters import gaussian_filter
import sys
 
#dataPath = "./simpleData"
dataPath = "./EData"
queryPath = "querySVG.svg"
sampleLen = 10 
maxError = 5 
binNum = 60
#Cut the line into several line segments
# input : a path 
# output : a path
def lineSegment(path):
        newPath = Path()
        for segment in path:
                if type(segment) is  Line :
                        newPath.append(segment)
                if type(segment) is  CubicBezier :
                        n = int(math.floor(segment.length()/sampleLen)) + 1
                        for i in range(n) :
                                newPath.append(Line(segment.point(1.0*i/n), segment.point(1.0*(i+1)/n)))
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
def preprocess(pathList, smooth = False):
        newPath = Path()
        for path in pathList:
		path = pathFitCurve(path, maxError, smooth)
                for segment in path:
                        newPath.append(segment)
        newPath = lineSegment(newPath)
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

def calcDistance(f1, f2):
	dist = 0
	for i in range(len(f1)):
		diff = f1[i] - f2[i]
		dist += diff * diff
	return dist

def calcDistance4(f1, f2):
	dist = 0
	for i in range(len(f1)):
		diff = f1[i] - f2[i]
		if diff < 0:
			diff =  -diff
		w = int(math.floor(diff*1000)) + 1
		diff = diff * w * w
		dist = dist + diff
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

	
def readDataSet(path, smooth=False, binNum=60):
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
		svgData.append(preprocess(svgData0[i], smooth)) 
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
		                                    gaussian_filter1d(f, 0.8, truncate=3), p, files[i]])
	return outputData

def distF(f1, f2):
	target = 0 
	targetDist = 10000000 
	for i in range(binNum):
		dist = calcDistance(np.roll(f1, i), f2)
		if dist < targetDist:
			target = i
			targetDist = dist
	return [targetDist, target]
def graph(g,binNum=60):
	n, bins, patches = plt.hist(g[0], binNum, weights=g[1])
	plt.show()

path = svg2paths(join(dataPath, sys.argv[1]))[0]
sys.argv[1]
