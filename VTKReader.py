import sys
import os
import numpy as np
from vtk import *
from vtk.util.numpy_support import vtk_to_numpy


def CompleteTikzFile(tikzFileName, Axis, CoordOverline, DataName, Data):
	
	ABC = ['a','b','c','d','f','g']
	
	#tikzFileName = BasePath + DataName
	
	f = open(tikzFileName,'w')
	f.write('\\RequirePackage{luatex85}\n')
	f.write('\\documentclass[tikz,pgf,border=1mm,14pt]{standalone}\n')
	f.write('\\usepackage{amssymb}\n')
	f.write('\\usepackage{amsmath}\n')
	f.write('\\usepackage{unicode-math}\n')
	f.write('\\usepackage{fontspec}\n')
	f.write('\n')
	f.write('\\setmainfont{Linux Libertine O}\n')
	f.write('\\defaultfontfeatures{Ligatures=TeX}\n')
	f.write('\\setsansfont[\n')
	f.write('BoldFont=LinBiolinumOB,\n')
	f.write('ItalicFont = LinBiolinumOI,\n')
	f.write('BoldItalicFont = LinLibertineOBI,\n')
	f.write(']{Linux Biolinum O}\n')
	f.write('\\setmathfont{TeXGyrePagellaMath}\n')
	f.write('\\setmathfont[range=\mathup]{Linux Libertine O}\n')
	f.write('\\setmathfont[range=\mathit/{latin,Latin,num}]{Linux Libertine O Italic}\n')
	f.write('\n')
	f.write('\\usetikzlibrary{datavisualization.formats.functions}\n')
	f.write('\\usepackage{gnuplottex}\n')
	f.write('\\usepackage[russian]{babel}\n')
	f.write('\n')
	f.write('\\begin{document}\n')
	f.write('\\begin{tikzpicture}\n')
	f.write('\\datavisualization [\n')
	f.write('scientific axes,\n')
	
	if Axis == '1':
		f.write('x axis={length=7.0cm, label={$x_2$, мм}},\n')
		tmp = '(\\bar{x}_1,x_2)'

	if Axis == '2':
		f.write('x axis={length=7.0cm, label={$x_1$, мм}},\n')
		tmp = '(x_1,\\bar{x}_2)'
	
	if DataName == 'U_X':
		f.write('y axis={length=5.0cm, label={$u_{1}'+tmp+'$, мм}},\n')
		
		
	if DataName == 'U_Y':
		f.write('y axis={length=5.0cm, label={$u_{2}'+tmp+'$, мм}},\n')
		
	if DataName == 'EpsiXX':
		f.write('y axis={length=5.0cm, label={$\\upvarepsilon_{11}'+tmp+'$}},\n')
		
	if DataName == 'EpsiXY':
		f.write('y axis={length=5.0cm, label={$\\upvarepsilon_{12}'+tmp+'$}},\n')
		
	if DataName == 'EpsiYY':
		f.write('y axis={length=5.0cm, label={$\\upvarepsilon_{22}'+tmp+'$}},\n')
		
	if DataName == 'SigmaXX':
		f.write('y axis={length=5.0cm, label={$\\upsigma_{11}'+tmp+'$, МПа}},\n')
		
	if DataName == 'SigmaXY':
		f.write('y axis={length=5.0cm, label={$\\upsigma_{12}'+tmp+'$, МПа}},\n')
		
	if DataName == 'SigmaYY':
		f.write('y axis={length=5.0cm, label={$\\upsigma_{22}'+tmp+'$, МПа}},\n')
		
	f.write('all axes={grid},\n')	
	f.write('legend ={north inside},\n')	
	
	tmpstr = '{'
	for i in range(len(Data)):
		tmpstr = tmpstr +  ABC[i] + ','
	tmpstr =  tmpstr[:-1] + '}'
	f.write('visualize as line/.list=' + tmpstr + ',\n')
	
	for i in range(len(Data)):
		f.write(ABC[i]+'={label in legend={text={$p_1=1.0$}, label in legend line coordinates={(-0.5em,0),(0,0)}}},\n')
	
	f.write('style sheet= vary dashing,\n')
	f.write(']\n')
	
	for i in range(len(Data)):
		
		tmpstr = 'data[' +'set=' + ABC[i] + '] {'
		f.write(tmpstr)
		f.write('\n')
		f.write(' x, y\n')
		
		tmp = Data[i]
		
		for j in range(np.size(CoordOverline)):
			f.write('%3.1f, %E\n' % (CoordOverline[j], tmp[j]))			
			
		f.write('\n')
		f.write('\n}\n')
		
	f.write(';\n')	
	f.write('\\end{tikzpicture}\n')
	f.write('\\end{document}')
		
	f.close()


def vtkArrraysOverLine(path, vtkFileNames, Axis, Coord, vtkArrayName):
	# load avtk file as input
	
	vtkData = []
	
	for i in range(len(vtkFileNames)):
		reader = vtk.vtkUnstructuredGridReader()
		reader.SetFileName(vtkFileNames[i])
		reader.ReadAllScalarsOn()
		reader.Update()  # Needed because of GetScalarRange
		vtkData.append(reader.GetOutput())


	vtkPoints = vtkData[0].GetPoints()


	X = np.zeros(vtkPoints.GetNumberOfPoints())
	Y = np.zeros(vtkPoints.GetNumberOfPoints())
	Z = np.zeros(vtkPoints.GetNumberOfPoints())

	for i in range(vtkPoints.GetNumberOfPoints()):
		X[i],Y[i],Z[i] = vtkPoints.GetPoint(i)

	vtkArray = []
	DataArr = []
	
	for i in range(len(vtkFileNames)):
		vtkArray.append(vtkData[i].GetPointData().GetArray(vtkArrayName))
		DataArr.append( vtk_to_numpy(vtkArray[i]))

	
	#Get the coordinates of the nodes and the scalar values

	Overline = []
	CoordOverline = []
	
	if Axis == '1':
		Overline = np.nonzero(X == float(Coord))
		CoordOverline = Y[Overline]	
	if Axis == '2':
		Overline = np.nonzero(Y == float(Coord))
		CoordOverline = X[Overline]	
		
	sortIdx = np.argsort(CoordOverline)
	
	DataArrOverLine = []
	
	
	for i in range(len(vtkFileNames)):
		tmp = DataArr[i]
		tmp = tmp[Overline]
		DataArrOverLine.append(tmp)
		tmp = tmp[sortIdx]
		DataArrOverLine[i] = tmp

	CoordOverline = CoordOverline[sortIdx]
	
	print(len(DataArrOverLine))
	
	if not os.path.exists(path +"TiKZ/"):
		os.makedirs(path +"TiKZ/")
	
	tikzFileName = path +"TiKZ/" +vtkArrayName +"_"+ Axis + "_" + Coord +".tex";
	
	print(tikzFileName)

	CompleteTikzFile(tikzFileName, Axis, CoordOverline, vtkArrayName, DataArrOverLine)

# ============================================================
# ============================================================
# ============================================================


def FindFile(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

#vtkArrrayOverLine(sys.argv[1], sys.argv[2],  sys.argv[3], sys.argv[4])

path = os.path.dirname(os.path.realpath('__file__'))

#path += "/bin/VTK/tens_beam_4_5000/"

#pArr = ['1.0', '0.9', '0.8']

#vtkFileNames = [path + "tens_beam_4_5000.vtk", path + pArr[1] + "/" +  "tens_beam_4_5000_0.9_100.0.vtk", path + pArr[2] + "/" +  "tens_beam_4_5000_0.8_100.0.vtk"]


pArr = ['1.0', '0.9', '0.8']

#path += "/bin/VTK/tens_beam_4_5000/" 

#NewPath  = path + "/" + pArr[1] + "/"

#vtkFileNames = [path + "tens_beam_4_5000.vtk", path + pArr[1] + "/" +  "tens_beam_4_5000_0.9_100.0.vtk", path + pArr[1] + "/" +  "tens_beam_4_5000_0.9_50.00.vtk"]


SearchPath = "/home/anshu/Документы/NonLocalStaticAnalysis/bin/VTK/"
vtkFileNames = ["tens_beam_8_5000.vtk", "tens_beam_8_5000_0.9_50.00.vtk", "tens_beam_8_5000_0.9_50.00_fmm.vtk"]

vtkPathToFiles = []

for i in range(len(vtkFileNames)):
	vtkPathToFiles.append(FindFile(vtkFileNames[i], SearchPath))
	
tmp = vtkPathToFiles[0]

NewPath = tmp[:-len(vtkFileNames[0])]

if len(sys.argv) > 2:
	vtkArrraysOverLine(NewPath, vtkPathToFiles, sys.argv[1], sys.argv[2], sys.argv[3])
	
	
	

	
# ============================================================
# ============================================================
# ============================================================
	
	
	
	
	
	
	
