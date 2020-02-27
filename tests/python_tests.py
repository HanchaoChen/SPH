import vtk
import os
import numpy as np

def test_position_outboundary():
    reader = vtk.vtkXMLPolyDataReader()
    path = "data/"
    #path = "../data/"
    files = os.listdir(path)
    for file in files:
            if not os.path.isdir(file):
                 reader.SetFileName(path + file)
                 reader.Update()
                 pdata = reader.GetOutput()
                 points = pdata.GetNumberOfPoints()
                 x = np.zeros(points)
                 y = np.zeros(points)
                 z = np.zeros(points)
                 for i in range(0, points):
                         x[i],y[i],z[i] = pdata.GetPoint(i)
                         assert ((x[i] <= 20.4) and (x[i] >= -0.4))
                         assert ((y[i] <= 10.4) and (y[i] >= -0.4))
                         assert (z[i] == 0.0)

def test_position_overlap():
    reader = vtk.vtkXMLPolyDataReader()
    path = "data/"
    #path = "../data/"
    files = os.listdir(path)
    for file in files:
            if not os.path.isdir(file):
                 reader.SetFileName(path + file)
                 reader.Update()
                 pdata = reader.GetOutput()
                 points = pdata.GetNumberOfPoints()
                 for i in range(0, points):
                        assert (pdata.GetPointData().GetArray('If_topped').GetValue(i) == 0.0)

def test_density():
    reader = vtk.vtkXMLPolyDataReader()
    path = "data/"
    #path = "../data/"
    files = os.listdir(path)
    for file in files:
            if not os.path.isdir(file): 
                 reader.SetFileName(path + file)
                 reader.Update()
                 pdata = reader.GetOutput()
                 points = pdata.GetNumberOfPoints()
                 for i in range(0, points):
                          assert (pdata.GetPointData().GetArray('Density').GetValue(i) > 0.0)



def test_file_writer_output():
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName('tests/test_1.vtp')
    reader.Update()

    pdata = reader.GetOutput()
    assert pdata.GetNumberOfCells() == 10
    assert pdata.GetNumberOfPoints() == 10
    assert pdata.GetPoint(5) == (5.0, 5.0, 0.0)
    assert (pdata.GetPointData().GetArray('Velocity').GetTuple(5)
            == (5.0, -5.0, 0.0))
    assert (pdata.GetPointData().GetArray('Pressure').GetValue(5)
            == 5.0)
    
#test_position()
#test_density()
