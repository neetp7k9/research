from os import listdir
from os.path import isfile, join
from svgpathtools import Path, Line, QuadraticBezier, CubicBezier, Arc, svg2paths, wsvg, disvg
import math 
import cmath 

dataPath = "./IconData"
queryPath = "querySVG.svg"
sampleNum = 10

def lineSegment(path):
        newPath = Path()
        for segment in path:
                if type(segment) is  Line :
                        newPath.append(segment)
                if type(segment) is  CubicBezier :
                        for i in range(sampleNum) :
                                newPath.append(Line(segment.point(i/10.0), segment.point((i+1)/10.0)))
        return newPath  

def preprocess(pathList):
        newPath = Path()
        for path in pathList:
                for segment in path:
                        newPath.append(segment)
        newPath = lineSegment(newPath)
	return newPath

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
	return radiansFeature

def calcWeightFeature(path):
	weightFeature = [ segment.length() for segment in path]
	return weightFeature
#read file according to path
files = [f for f in listdir(dataPath) if isfile(join(dataPath, f))]
if ".DS_Store" in files :
	files.remove(".DS_Store") 

#read svg from file
svgData = [ preprocess(svg2paths(join(dataPath, p))[0]) for p in files]
queryData = preprocess(svg2paths(queryPath)[0])


radiansFeature = [ calcRadiansFeature(p) for p in svgData]
queryRadiansFeature = calcRadiansFeature(queryData) 

