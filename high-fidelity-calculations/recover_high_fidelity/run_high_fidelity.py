import sys
import os
import random
import math
import copy
import time
import numpy as np
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import dolfin as dl
import hippylib as hl
import scipy.stats as sstats
from scipy.io import loadmat
import argparse

# import hp model 
sys.path.insert(0, '../../../recover_high_fidelity/')
from hp_model_class import *

## to run
## 1. conda activate confen
## 2. run script using
## mpirun -n 6 python run_high_fidelity.py

# get comm and rank
comm = dl.MPI.comm_world
rank = comm.rank

def print_msg(msg, fls = True):
  if rank == 0:
    print(msg, flush = fls)

if __name__ == "__main__":

  ## create model (use partial hippy lib framework)

  # load mesh and subdomain
  print_msg('loading mesh and subdomains', True)
  fmesh = '../../../mesh/mesh4'
  mesh = dl.Mesh(fmesh + '.xml')
  subdomain = dl.MeshFunction("size_t", mesh, 3)
  dl.File(fmesh + '_subdomain.xml.gz') >> subdomain

  # x = mesh.coordinates()
  # x[:,:] *= VOXELSIZE
  #print('min and max coords: ', np.min(x), np.max(x))

  # FE space
  FE_polynomial = 2 # if different from 2, need to compute vasc and tiss dofs again
  Vu = dl.FunctionSpace(mesh, "Lagrange", FE_polynomial)

  # Define elements: P1 and real number
  P1  = dl.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
  R   = dl.FiniteElement("R", mesh.ufl_cell(), 0)

  # Define the state space (solution space) together with parameter space
  Vh = [Vu, None, None]

  # HP MRI parameters from the OED results
  print_msg('loading model and design parameters', True)
  def_params = loadmat('../../../def_params.mat')
  def_params = def_params['str']
  def_params = def_params['params'][0][0][0]

  # collect design params
  TRi = def_params['TRi'][0][0]
  print_msg('TRi = {}'.format(TRi))
  FaList = def_params['FaList'][0]

  # set OED design parameters flip angle pyr = 35, lac = 28
  FaList[0, :] = 35*np.pi/180.
  FaList[1, :] = 28*np.pi/180.
  print_msg('FaList = {}'.format(FaList))
  FaP = FaList[0, :]
  FaL = FaList[1, :]

  Nacq = len(TRi) + 1
  tsteps_per_acq = 20
  print_msg('Nacq = {}, tsteps_per_acq = {}'.format(Nacq, tsteps_per_acq))

  # collect model params
  DP = 20.
  DL = 20.
  LP = 0.2 #0.02
  LL = 0.2 #0.02
  kpl = def_params['ExchangeTerms'][0][0][1]
  klp = 1.e-16 #oed_data_params['ExchangeTerms'][0][1][0]
  T1P = def_params['T1s'][0][0][0]
  T1L = def_params['T1s'][0][1][0]
  kve = def_params['kve'][0][0][0]
  t0 = def_params['t0'][0][0][0]
  ve = def_params['ve'][0][0][0]
  gammaPdfA = def_params['gammaPdfA'][0][0][0]
  gammaPdfB = def_params['gammaPdfB'][0][0][0]
  scaleFactor = 10. #def_params['scaleFactor'][0][0][0]
  kve_ve_ref = kve / ve;

  print_msg('TRi = {}'.format(TRi))
  print_msg('FaP = {}'.format(FaP))
  print_msg('FaL = {}'.format(FaL))
  print_msg('DP = {}, DL = {}, \n'
            'LP = {}, LL = {}, \n'
            'kpl = {}, klp = {}, \n'
            'T1P = {}, T1L = {}, \n'
            't0 = {}, kve = {}, ve = {}, \n'
            'gammaPdfA = {}, gammaPdfB = {}, \n'
            'scaleFactor = {}'.format(DP, DL, LP, LL, kpl, klp, \
              T1P, T1L, t0, kve, ve, gammaPdfA, gammaPdfB, scaleFactor))

  # create pde problem
  print_msg('creating forward problem', True)
  fwd_result = './'
  pde = DiffReactHPModel(Vh, TRi, FaP, FaL, \
                   subdomain, \
                   tsteps_per_acq = tsteps_per_acq, nspecies = 2, \
                   model_out_path = fwd_result, save_model_sol = True, \
                   sim_out_freq = 1, sim_max_out = 20,\
                   verbosity = 2, sim_print_freq = 1, comm = comm)

  ## Run model with OED design parameters

  # set model parameters (DP, DL, kPL, kLP, T1P, T1L, LP, LL) in log-space
  print_msg('setting model parameters', True)
  params = np.log(np.array([DP, DL, kpl, klp, T1P, T1L, LP, LL]))
  pde.set_parameters_general(params)
  pde._VIF_t0 = t0
  pde.sigmaP = scaleFactor
  pde.alphaP = gammaPdfA
  pde.betaP = gammaPdfB
  pde.set_VIF()

  # run model
  print_msg('running forward problem', True)
  qoi_oed = pde.generate_state()
  pde.solve(qoi_oed)
  if rank == 0:
    np.save(fwd_result + 'qoi.npy', qoi_oed)
