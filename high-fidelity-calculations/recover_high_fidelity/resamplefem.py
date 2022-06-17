import os
import configparser
import vtk
import numpy
import subprocess
import math
import argparse

c3d = "/Applications/Convert3DGUI.app/Contents/bin/c3d"

## Usage
## conda activate hpmri
## python resamplefem.py

def resample(data_path, pixelsize = 16., num_scans = 30):
  #if (options.resample != None and  options.outputid != None ):
  # resample back to image
  print("resampling...")

  fileList = [data_path + 'pyruvate%06d.pvtu', data_path + 'lactate%06d.pvtu', data_path + 'pyruvatevasc%06d.pvtu']

  for idspecies in fileList:
    for idfile in range(num_scans):
      femoutputfile = idspecies % idfile 
      print('file to be processed = {}'.format(femoutputfile))
      vtkFemReader = vtk.vtkXMLPUnstructuredGridReader() 
      vtkFemReader.SetFileName(femoutputfile )
      vtkFemReader.Update() 
      
      vtkOutline  = vtk.vtkOutlineFilter()
      vtkOutline.SetInputData (vtkFemReader.GetOutput() ) 
      vtkOutline.Update()
      vtkbb = vtkOutline.GetOutput()
      femoutline = vtkbb.GetBounds()

      imageorigin = (femoutline[0],femoutline[2],femoutline[4])
      imagespacing= ((femoutline[1]-femoutline[0])/(pixelsize-1),(femoutline[3]-femoutline[2])/(pixelsize-1),(femoutline[5]-femoutline[4])/(pixelsize-1))
      
      generateimageCMD = c3d + " -create %dx%dx%d %fx%fx%fmm -origin %fx%fx%fmm -o %s/femgrid.vtk" % (pixelsize,pixelsize,pixelsize ,imagespacing[0],imagespacing[1],imagespacing[2],imageorigin[0], imageorigin[1], imageorigin[2],data_path)
      print(generateimageCMD )
      if (not os.path.isfile(data_path + '/femgrid.vtk')):
        os.system( generateimageCMD )
      
      vtkImageReader = vtk.vtkDataSetReader() 
      vtkImageReader.SetFileName(data_path + '/femgrid.vtk')
      vtkImageReader.Update() 
      
      vtkResample = vtk.vtkCompositeDataProbeFilter()
      vtkResample.SetSourceData( vtkFemReader.GetOutput() )
      vtkResample.SetInputData( vtkImageReader.GetOutput() ) 
      vtkResample.Update()
      #
      # write to disk
      # outputbase = femoutputfile.split('.').pop(0)
      outputbase = femoutputfile.replace('.pvtu', '')
      print('outputbase = {}'.format(outputbase))
      # exit()
      vtkFEMImageWriter = vtk.vtkDataSetWriter() 
      vtkFEMImageWriter.SetFileTypeToBinary() 
      vtkFEMImageWriter.SetInputData( vtkResample.GetOutput() )
      vtkFEMImageWriter.SetFileName( '%s.vtk' % outputbase  )
      vtkFEMImageWriter.Update() 
      
      # compress
      compressCMD = c3d + " %s.vtk -o %s.nii; " % (outputbase  ,outputbase  )
      print(compressCMD )
      os.system( compressCMD )

  if (not os.path.isfile(data_path + '/mask.nii')):
    maskCMD = c3d + ' -verbose ' + data_path + '/pyruvate00000?.nii -foreach -binarize -endfor -accum -add -endaccum  -binarize -o ' + data_path + '/mask.nii'
    print(maskCMD )
    os.system( maskCMD ) 
  if (not os.path.isfile(data_path + '/pyrlacmip.nii')):
    mipCMD = c3d + ' -verbose ' + data_path + '/pyruvate0000??.nii ' + data_path + '/lactate0000??.nii -accum -add -endaccum -threshold 95% inf 1 0 -o ' + data_path + '/pyrlacmip.nii'
    print(mipCMD )
    os.system( mipCMD ) 

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Run MCMC using high-fidelity QoIs as data')
  parser.add_argument('--data_path',
                      default='./0/',
                      type=str,
                      help="Path to pvtu data files")

  parser.add_argument('--pixelsize',
                      default=16.,
                      type=float,
                      help="Pixelsize")

  parser.add_argument('--num_scans',
                      default=30,
                      type=float,
                      help="Number of scans")

  args = parser.parse_args()
  
  resample(args.data_path, args.pixelsize, args.num_scans)
