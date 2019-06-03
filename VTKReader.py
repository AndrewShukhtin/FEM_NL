import sys
import os
import numpy as np
from vtk import *
from vtk.util.numpy_support import vtk_to_numpy


class VTKtoTIKZ:
		
	def __init__(self, vtkFileNames, SearchPath):
		self.__vtkFileNames = vtkFileNames
		self.__SearchPath = SearchPath
		self.__vtkPathToFiles = []
	
	def FindFiles(self):
		for i in range(len(self.__vtkFileNames)):
			for root, dirs, files in os.walk(self.__SearchPath):
				if self.__vtkFileNames[i] in files:
					self.__vtkPathToFiles.append( os.path.join(root, self.__vtkFileNames[i]) )
					
					
	def GenerateTikzFile(self, Axis, Coord, vtkArrayName):
		
		tmp = self.__vtkPathToFiles[0]
		self.__path = tmp[:-len(self.__vtkFileNames[0])]
		
		self.__Axis = Axis
		self.__Coord = Coord
		self.__vtkArrayName = vtkArrayName
		
		print(self.__path)
		
		CoordArrOverline = []
		DataArrOverLine = []
		
		for i in range(len(self.__vtkFileNames)):
			
			reader = vtk.vtkUnstructuredGridReader()
			reader.SetFileName(self.__vtkPathToFiles[i])
			reader.ReadAllScalarsOn()
			reader.Update()
			vtkData = reader.GetOutput()

			
			vtkPoints = vtkData.GetPoints()


			X = np.zeros(vtkPoints.GetNumberOfPoints())
			Y = np.zeros(vtkPoints.GetNumberOfPoints())
			Z = np.zeros(vtkPoints.GetNumberOfPoints())


			for j in range(vtkPoints.GetNumberOfPoints()):
				X[j],Y[j],Z[j] = vtkPoints.GetPoint(j)		
		
		
			vtkArray = vtkData.GetPointData().GetArray(vtkArrayName)
			DataArr = (vtk_to_numpy(vtkArray))
		
		
			if self.__Axis == '1':
				Overline = np.nonzero(X == float(Coord))
				CoordOverline = Y[Overline]	
			if self.__Axis == '2':
				Overline = np.nonzero(Y == float(Coord))
				CoordOverline = X[Overline]	
		
			sortIdx = np.argsort(CoordOverline)
			
			DataOverLine = DataArr[Overline]
			DataOverLine = DataOverLine[sortIdx]
			DataArrOverLine.append(DataOverLine)

			CoordOverline = CoordOverline[sortIdx]
			CoordArrOverline.append(CoordOverline)	
		
		
		if not os.path.exists(self.__path +"/TiKZ/"):
			os.makedirs(self.__path +"/TiKZ/")
		
		
		self.__tikzFileName = self.__path  +"TiKZ/" +vtkArrayName +"_"+ Axis + "_" + Coord +".tex";
		
		print(self.__tikzFileName)

		self.__WriteToTikzFile(CoordArrOverline, vtkArrayName, DataArrOverLine)
		
				
	
	def __WriteToTikzFile(self, CoordArrOverline, DataName, DataArr):
		
		ABC = ['a','b','c','d','f','g']
		
		f = open(self.__tikzFileName,'w')
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
		
		if self.__Axis == '1':
			f.write('x axis={length=7.0cm, label={$x_2$, мм}},\n')
			tmp = '(\\bar{x}_1,x_2)'

		if self.__Axis == '2':
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
		for i in range(len(DataArr)):
			tmpstr = tmpstr +  ABC[i] + ','
		tmpstr =  tmpstr[:-1] + '}'
		f.write('visualize as line/.list=' + tmpstr + ',\n')
		
		for i in range(len(DataArr)):
			f.write(ABC[i]+'={label in legend={text={$p_1=1.0$}, label in legend line coordinates={(-0.5em,0),(0,0)}}},\n')
		
		f.write('style sheet= vary dashing,\n')
		f.write(']\n')
		
		for i in range(len(DataArr)):
			
			tmpstr = 'data[' +'set=' + ABC[i] + '] {'
			f.write(tmpstr)
			f.write('\n')
			f.write(' x, y\n')
			
			Data = DataArr[i]
			Coords = CoordArrOverline[i]
			
			for j in range(np.size(Coords)):
				f.write('%3.1f, %E\n' % (Coords[j], Data[j]))
				
			f.write('\n')
			f.write('\n}\n')
			
		f.write(';\n')	
		f.write('\\end{tikzpicture}\n')
		f.write('\\end{document}')
			
		f.close()


def main():

	path = os.path.dirname(os.path.realpath('__file__'))

	SearchPath = "/home/anshu/Документы/NonLocalStaticAnalysis/bin/VTK/"
	vtkFileNames = ["tens_beam_4_5000.vtk", "tens_beam_4_5000_0.9_50.00.vtk", "tens_beam_4_5000_0.9_50.00_fmm.vtk"]

	tikz = VTKtoTIKZ(vtkFileNames, SearchPath)
	tikz.FindFiles()

	if len(sys.argv) > 2:
		tikz.GenerateTikzFile(sys.argv[1], sys.argv[2], sys.argv[3])
		

	
if __name__ == '__main__':
	main()
	

	
# ============================================================
# ============================================================
# ============================================================
	
	
	
	
	
	
	
