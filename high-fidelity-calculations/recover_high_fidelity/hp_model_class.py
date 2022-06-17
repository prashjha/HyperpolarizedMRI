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
from pathlib import Path

STATE = 0
PARAMETER = 1

PYR = 0
LAC = 1
PYRV = 2
LACV = 3

TISSDOM = 0
VASCDOM = 1

VOXELSIZE = 1.

PENALTY = 10000.

def SubmeshDofs(submesh, V):
    'Return dofs of V that belong to cells in submesh.'
    mesh = V.mesh()
    tdim = mesh.topology().dim()
    dofmap = V.dofmap()

    submesh_dofs = set()
    parent_cell_indices = submesh.data().array('parent_cell_indices', tdim)
    for index in range(submesh.num_cells()):
        cell = parent_cell_indices[index]
        [submesh_dofs.add(dof) for dof in dofmap.cell_dofs(cell)]

    return np.array(list(submesh_dofs))

## Helper class
class GenericReprBase:
  
  def get_info(self,delim = '\n '):
    name = type(self).__name__
    
    vars_list = []
    for key, value in vars(self).items():
        if isinstance(value, np.ndarray):
            value = value.flatten()
        vars_list.append(f'{key}={value!r}')
        
    vars_str = delim.join(vars_list)
    return f'{name}({vars_str})'

  def __repr__(self,delim = '\n '):
    return self.get_info(delim) 


dl.parameters["form_compiler"]["optimize"]     = True
dl.parameters["form_compiler"]["cpp_optimize"] = True
dl.parameters["form_compiler"]["representation"] = "uflacs"
dl.parameters["form_compiler"]["cpp_optimize_flags"] = "-O3 -ffast-math -march=native"

#https://fenicsproject.org/qa/13062/mesh-partition-parameters-using-mpi/
dl.parameters["mesh_partitioner"] = "ParMETIS"

class VIFSource(dl.UserExpression):
  def __init__(self,t0,s,a,b,fa, **kwargs):
    super().__init__(**kwargs)
    self.t = 0.
    self.flip_P = fa
    self.t0 = t0
    self.s = s
    self.a = a
    self.b = b
  
  def get_val(self):
    #return np.cos(self.flip_P) * self.s * (1./self.b) * sstats.gamma.pdf((self.t - self.t0)/self.b, self.a)
    return self.s * (1./self.b) * sstats.gamma.pdf((self.t - self.t0)/self.b, self.a)

  def eval(self, values, x):
    values[0] = self.get_val()

class DiffReactHPModel(GenericReprBase):
  """
  Vh: List of function vector space and parameter space
  simulation_time: Final simulation time
  TR: Vector of time interval between successive scans. Size = N - 1 where N is the number of scans
  FaP: Flip angles at each scane for pyruvate. Size = N where N is the number of scans
  FaL: Flip angles at each scane for lactate. Size = N where N is the number of scans
  subdomain: Vascular and tissue subdomain data
  tsteps_per_acq: Number of time steps for simulation between two scans
  nspecies: Number of species (only 2 is supported presently)
  model_out_path: Path to directory where pde solution should be saved
  save_model_sol: True means we save the pde solution
  save_model_sol_step: Output model state after every 'save_model_sol_step' time steps
  sim_out_freq: Output pde state for sim when 'sim_out_freq' % 'sim_count' = 0. 'sim_count' is the
                sim number of current sim. We increment the 'sim_count' everytime the forward solve 
                is called with new parameter.
  sim_max_out: Total number of sims to be saved
  tol: Tolerance for picard iteration
  max_iter: Maximum iteration in picard iteration
  verbosity: Level of information output
  sim_print_freq: Print information for sim when 'print_freq' % 'sim_count' = 0
  comm: MPI comm in case this is run in parallel
  """
  def __init__(self, Vh, TR, FaP, FaL, \
               subdomain, \
               tsteps_per_acq = 10, nspecies = 2, \
               model_out_path = './fwd_result/', save_model_sol = False, save_model_sol_step = 2, \
               sim_out_freq = 10, sim_max_out = 20,\
               tol = 1.0e-12, max_iter = 100, \
               verbosity = 1, sim_print_freq = 1, comm = None):
      
    ## comm
    self.comm = comm
    self.rank = 0
    self.size = 0
    if self.comm is not None:
      self.rank = self.comm.rank
      self.size = self.comm.size
    
    ## Vh groups parameter space and function space
    self.Vh = Vh 
    self.mcmc_use = False
    if Vh[PARAMETER] is not None:
      self.mcmc_use = True
    
    ## control parameters
    self.Nacq = len(TR) + 1 # number of scans is number of interval plus 1
    self.TR = TR
    self.FaP = FaP
    self.FaL = FaL
    self.tsteps_per_acq = tsteps_per_acq
    
    if len(self.FaP) != len(self.FaL) or len(FaP) != self.Nacq:
      raise ValueError("Flip angles data is invalid.")
    
    self.time = 0. # current time of simulation
    self.dt = 0.0001
    self.dt_expr = dl.Constant(self.dt) # dolfin constant for dt as dt changes during the simulation
    
    ## number of species
    self.nspecies = nspecies #NSPECIES
    
    ## output related parameters
    self.model_out_path = model_out_path
    self.sim_model_out_path = None
    self.save_model_sol = save_model_sol
    self.save_model_sol_step = save_model_sol_step
    
    self.file_model_sol = None
    self.file_model_sol_vasc = None
    self.var_file_names = ['pyruvate', 'lactate', 'pyruvate_v', 'lactate_v']
    self.var_names = ['Pyruvate', 'Lactate', 'Pyruvate V', 'Lactate V']
    
    self.sim_count = 0 # current number of sim (how many forward solves have been called so far)
    self.sim_out_freq = sim_out_freq # useful when interfacing this with hipplib mcmc
    self.sim_max_out = sim_max_out   # restrict total number of independent simulation outputs 
    self.sim_out_current = 0
    if self.mcmc_use:
      self.sim_out_freq = 1
      self.sim_max_out = 10000

    # count total output per simulation
    self.out_count = 0
    
    # verbosity 
    self.verbosity = verbosity
    self.sim_print_freq = sim_print_freq if self.mcmc_use else 1
    
    ## parameters
    if self.mcmc_use:
      self.param_dim = Vh[PARAMETER].dim()
      self.m = dl.Function(Vh[PARAMETER])
    else:
      self.param_dim = None
      self.m = None
    
    self.DP = dl.Constant(math.log(2.))
    self.DL = dl.Constant(math.log(2.))
    self.kPL = dl.Constant(math.log(0.1))
    self.kLP = dl.Constant(math.log(0.00001))
    self.T1P = dl.Constant(math.log(43.))
    self.T1L = dl.Constant(math.log(33.))
    
    # fixed parameters
    self.LP = dl.Constant(math.log(0.02))
    self.LL = dl.Constant(math.log(0.02))
    self._VIF_t0 = 0.  # offset by arbitrarily fixed TR_start value
    self.sigmaP, self.alphaP, self.betaP = 10., 2.5, 4.5
    
    # different ways to consider VIF source
    # use gamma function (must set 't' before assembly)
    self.VIF_fn_P = None
    self.set_VIF()

    ## solver properties
    if self.size > 0:
      self.solver = dl.PETScKrylovSolver('gmres', 'hypre_euclid')
    else:
      self.solver = dl.PETScKrylovSolver("gmres", "ilu")
    self.solver.parameters["nonzero_initial_guess"] = True
    self.tol = tol
    self.max_iter = max_iter
    self.explicit_time_solver = True
    
    ## state variables
    self.u_0 = [dl.Function(self.Vh[STATE]) for i in range(self.nspecies)]
    self.u = [dl.Function(self.Vh[STATE], name = self.var_names[i]) for i in range(self.nspecies)]
    self.uv = [dl.Function(self.Vh[STATE], name = self.var_names[i+2]) for i in range(self.nspecies)]
    
    ## test and trial functions
    self.help = dl.Function(Vh[STATE])
    self.p = dl.TestFunction(self.Vh[STATE])
    self.u_trial = dl.TrialFunction(self.Vh[STATE])
    
    # mass matrix
    self.M = dl.assemble(self.u_trial*self.p*dl.dx)
    
    # for integral QoI calculation
    self.z = dl.assemble(dl.Constant(1.0)*self.p*dl.dx)
    
    ## bilinear and linear forms
    self.F = [None]*(self.nspecies+1)

    self.H = [None]*(self.nspecies+1)
    self.b = [None]*(self.nspecies+1)
    
    self.a = [None]*(self.nspecies+1)
    self.L = [None]*(self.nspecies+1)
    
    ## subdomain
    self.mesh = Vh[STATE].mesh()
    self.subdomain = subdomain
    self.dx_sub = dl.Measure('dx', domain=self.mesh, subdomain_data=self.subdomain)
    
    self.volume_tiss = dl.assemble(dl.Constant(1.0)*dl.Constant(1.0)*self.dx_sub(TISSDOM)) 
    self.volume_vasc = dl.assemble(dl.Constant(1.0)*dl.Constant(1.0)*self.dx_sub(VASCDOM))
    self.volume = self.volume_tiss + self.volume_vasc
    self.ev_vol_ratio = self.volume_tiss / self.volume

      
  def print_info(self, to_screen = True, to_file = False, fname = 'model.txt'):
    if self.rank == 0:
      info_str = self.get_info()
      if to_file:
        f = open(fname, 'w')
        f.write(info_str)
        f.close()
      if to_screen:
        print(info_str)
          
  def set_VIF_t0(self, t0):
    self._VIF_t0 = t0
      
  def print_msg(self, msg, lvl):
    if self.rank == 0:
      if self.sim_count % self.sim_print_freq == 0 and lvl <= self.verbosity:
        print(msg, flush=True)
          
  def set_VIF(self):
    if self.FaP is not None and len(self.FaP) >= 1:
      self.VIF_fn_P = VIFSource(t0 = self._VIF_t0, s = self.sigmaP, \
                                a = self.alphaP, b = self.betaP,\
                                fa = self.FaP[0], degree=2)
    else:
      raise ValueError("Flip angle vector is either None or has no elements.")
          
  def TR_to_time(self, TR):
    t = []
    ct = 0.
    for a in TR:
      ct += a
      t.append(ct)
      
    return t

  def reset_sol(self):

    for i in range(self.nspecies):
      self.u[i].vector().zero()
      self.u_0[i].vector().zero()
      self.u[i].vector().apply("")
      self.u_0[i].vector().apply("")

    for i in range(self.nspecies):
      self.uv[i].vector().zero()
      self.uv[i].vector().apply("")
          
  def generate_state(self):
    """ Return a vector in the shape of the fwd model (QoI) output. """
    return np.empty((self.nspecies, self.Nacq))

  def generate_pde_state(self):
    """ Return a list of vectors that correspons to the model state vectors. """
    return [dl.Function(self.Vh[STATE]).vector() for i in range(self.nspecies)]

  def generate_parameter(self):
    """ Return a vector in the shape of the parameter. """
    if self.Vh[PARAMETER] is not None:
      return dl.Function(self.Vh[PARAMETER]).vector()
    else:
      raise ValueError("HippyLib parameter space is not supplied.")

  def init_parameter(self, m):
    """ Initialize the parameter. """
    dummy = self.generate_parameter()
    m.init( dummy.mpi_comm(), dummy.local_range() )
      
  def set_parameters(self, param):
    """ Replace parameters with new parameters. """

    if self.Vh[PARAMETER] is None:
      raise ValueError("HippyLib parameter space is not supplied. If running outside "
                         "hippylib framework set parameters "
                         "using set_parameters_general with list parameters as input. "
                         "Parameters are assumed "
                         "to be in log-space.")
    
    self.m.vector().zero()
    self.m.vector().axpy(1., param)
        
    if self.param_dim == 2:
      (self.DL, self.kPL) = dl.split(self.m)
    elif self.param_dim == 6:
      (self.DP, self.DL, self.kPL, self.kLP, self.T1P, self.T1L) = dl.split(self.m)
    else:
      raise IndexError("The parameter dimension should be 2 or 6.")
          
  def set_parameters_general(self, param):
    """ Replace parameters with new parameters. 
        If number of params
        1. is 2 then update DL and kPL
        2. is 6 then update DP, DL, kPL, kLP, T1P, T1L
        3. is 8 then update DP, DL, kPL, kLP, T1P, T1L, LP, LL
    """
    
    if len(param) == 2:
      self.DL.assign(param[0])
      self.kPL.assign(param[1])
    
    if len(param) == 6:
      self.DP.assign(param[0])
      self.DL.assign(param[1])
      self.kPL.assign(param[2])
      self.kLP.assign(param[3])
      self.T1P.assign(param[4])
      self.T1L.assign(param[5])
        
        
    if len(param) == 8:
      self.DP.assign(param[0])
      self.DL.assign(param[1])
      self.kPL.assign(param[2])
      self.kLP.assign(param[3])
      self.T1P.assign(param[4])
      self.T1L.assign(param[5])
      self.LP.assign(param[6])
      self.LL.assign(param[7])
        
          
  def assign_fn_vector(self, u1, u2, scale = None):
    """ Copy values of u2 vector into u1. """
    u1.vector().zero()
    if scale is not None:
      u1.vector().axpy(scale, u2.vector())
    else:
      u1.vector().axpy(1., u2.vector())

    u1.vector().apply("")
  
  def assign_fn_vectors(self, u1, u2, scale_vec = None):
    """ Copy values of u2 vector into u1. """
    for i in range(self.nspecies):
      if scale_vec is not None:
        self.assign_fn_vector(u1[i], u2[i], scale = scale_vec[i])
      else:
        self.assign_fn_vector(u1[i], u2[i], scale = 1.)
          
  def set_varf_forms(self):
    """ Set variational form. """
    
    if self.nspecies != 2:
      raise ValueError('Number of species must be two for this variational form')

    ## tissue
    self.F[PYR] = (self.u_trial-self.u_0[PYR])*self.p*dl.dx \
             + self.dt_expr*dl.exp(self.DP)*dl.inner(dl.grad(self.u_trial), dl.grad(self.p))*dl.dx \
             - self.dt_expr*(-1./dl.exp(self.T1P))*self.u_trial*self.p*dl.dx \
             - self.dt_expr*(-dl.exp(self.kPL))*self.u_trial*self.p*dl.dx \
             - self.dt_expr*dl.exp(self.kLP)*self.u[LAC]*self.p*dl.dx \
             - self.dt_expr*dl.exp(self.LP)*self.VIF_fn_P*self.p*self.dx_sub(VASCDOM) 

             # kpl on only tissue domain
             #- self.dt_expr*(-dl.exp(self.kPL))*self.u_trial*self.p*self.dx_sub(TISSDOM) \
             #- self.dt_expr*dl.exp(self.kLP)*self.u[LAC]*self.p*self.dx_sub(TISSDOM) \
    
    self.F[LAC] = (self.u_trial-self.u_0[LAC])*self.p*dl.dx \
             + self.dt_expr*dl.exp(self.DL)*dl.inner(dl.grad(self.u_trial), dl.grad(self.p))*dl.dx \
             - self.dt_expr*(-1./dl.exp(self.T1L))*self.u_trial*self.p*dl.dx \
             - self.dt_expr*(-dl.exp(self.kLP))*self.u_trial*self.p*dl.dx \
             - self.dt_expr*dl.exp(self.kPL)*self.u[PYR]*self.p*dl.dx

             # kpl on tissue domain only
             #- self.dt_expr*(-dl.exp(self.kLP))*self.u_trial*self.p*self.dx_sub(TISSDOM) \
             #- self.dt_expr*dl.exp(self.kPL)*self.u[PYR]*self.p*self.dx_sub(TISSDOM)

    self.F[2] = self.u_trial*self.p*dl.dx - self.VIF_fn_P*self.p*self.dx_sub(VASCDOM)
        
    for i in range(self.nspecies+1):
      self.a[i], self.L[i] = dl.lhs(self.F[i]), dl.rhs(self.F[i])
          
  def evaluate(self):
    """ Compute QoI. 
        QoI = (1/|\Omega|) * \int_\Omega \phi (x, t) dx
    """
    a = []
    for i in range(self.nspecies):
      a.append(self.u[i].vector().inner(self.z)/self.volume)
        
    return np.array(a)

  def set_save_files(self, save_fdir = None, save_vasc = True):
    """ Create output directory and files for saving state of model. """
    if save_fdir is not None:
      self.sim_model_out_path = save_fdir + '/'
    else:
      self.sim_model_out_path = self.model_out_path + str(self.sim_count) + '/'
        
    if self.rank == 0:
      Path(self.sim_model_out_path).mkdir(parents=True, exist_ok=True)
        
    self.file_model_sol = [None]*self.nspecies
    self.file_model_sol_vasc = [None]*self.nspecies
    
    for i in range(self.nspecies):
      self.file_model_sol[i] = dl.File(self.sim_model_out_path + \
                                         self.var_file_names[i] + '.pvd', "compressed")
    if save_vasc:
      self.file_model_sol_vasc[0] = dl.File(self.sim_model_out_path + \
                                         self.var_file_names[0] + 'vasc' + '.pvd', "compressed")
          
  def save_state(self, save_vasc = True):
    """ Save the state of pde system to file. """
    
    self.print_msg('  Saving model solution', 2)
    
    # write states to file
    for i in range (self.nspecies):
      self.file_model_sol[i] << (self.u[i], self.time)
    if save_vasc:
      self.file_model_sol_vasc[0] << (self.uv[0], self.time)

    self.out_count += 1
          
  def solve_vascular_fields(self):
    #self.upv_np[:] = self.VIF_fn_P.get_val()
    #self.upv.vector().vec()[self.vasc_dofs] = self.upv_np
    #self.upv.vector().vec()[self.vasc_dofs] = self.VIF_fn_P.get_val()

    # solve
    i = 2
    #self.print_msg(self.H[i], 0)
    #self.print_msg(self.b[i], 0)

    if self.H[i] is None:
      self.H[i] = dl.assemble(self.a[i])
    #else:
    #  dl.assemble(self.a[i], tensor = self.H[i])
    
    if self.b[i] is None:
      self.b[i] = dl.assemble(self.L[i])
    else:
      dl.assemble(self.L[i], tensor = self.b[i])      
    
    self.solver.solve(self.H[i], self.uv[0].vector(), self.b[i])

  def solve(self, out, save_fdir = None):
    """ 
       Takes in the specified parameter and solves the forward model and computes the QoI (model output). 
    """
    
    start = time.time()
    
    # reset the solution vector
    self.reset_sol()
    
    # set VIF function
    self.set_VIF()

    # time stepping
    self.time = 0.0
    self.out_count = 0
    saved_state_current_sim = False # did we save any of the time step of current sim

    # set variational form to get VIF at initial time (for output)
    self.VIF_fn_P.t = self.time
    self.set_varf_forms()
    self.solve_vascular_fields()

    # save initial model state
    if self.save_model_sol and self.sim_count % self.sim_out_freq == 0 \
        and self.sim_out_current < self.sim_max_out:

      # create save files (call once)
      if saved_state_current_sim == False:
        self.set_save_files(save_fdir, save_vasc = True)

      # update flag
      saved_state_current_sim = True
      
      # actually save files
      self.print_msg('  Output time: {}, Output number: {}'.format(self.time, self.out_count + 1), 2)
      self.save_state(save_vasc = True)

    # save at first scan (t=0 i.e. initial time)
    out[:, 0] = self.evaluate()
    
    # loop over TRs
    self.print_msg('loop over scans', 1)
    tacq = 0.
    for iacq in range(1, self.Nacq):

      # solution period is [t_{i-1}, t_i] 
      # tacq is t_i and time_diff is t_i - t_{i-1}
      # NOTE: TR[0] = second acq - first acq, TR[1] = third acq - second acq, etc
      time_diff = self.TR[iacq-1] 
      tacq = tacq + time_diff

      # set time step
      self.dt = time_diff/self.tsteps_per_acq
      self.dt_expr.assign(self.dt) # update dt in the assembly function

      self.print_msg('Solving at acq step: {}, tacq: {}, time diff: {}'.format(iacq, tacq, time_diff), 2)

      # flip angles ({i-1}th scan angle to modify the initial condition for ith scan)
      flip_P = self.FaP[iacq-1]
      flip_L = self.FaL[iacq-1]
      
      # reset old solution from the last solution
      self.assign_fn_vector(self.u_0[0], self.u[0], np.cos(flip_P))
      self.assign_fn_vector(self.u_0[1], self.u[1], np.cos(flip_L))
          
      # copy old solution to current solution
      self.assign_fn_vector(self.u[0], self.u_0[0], 1.)
      self.assign_fn_vector(self.u[1], self.u_0[1], 1.)
      
      # set flip angle in VIF function 
      # self.VIF_fn_P.flip_P = flip_P

      local_tstep_count = 0
      
      # loop over time steps
      for itstep in range(self.tsteps_per_acq):
          
        self.time = self.time + self.dt
        
        self.print_msg('time: {:6.4f}'.format(self.time), 2)

        # update old time step solution
        self.assign_fn_vectors(self.u_0, self.u)
        
        ## solve
        # step 1 - assemble systems
        
        # for transient source, set current time
        self.VIF_fn_P.t = self.time
        
        # since parameter is updated, we need to update the variational forms
        self.set_varf_forms()

        for i in range(self.nspecies):

          if self.H[i] is None:
            self.H[i] = dl.assemble(self.a[i])
          else:
            if local_tstep_count == 0:
              dl.assemble(self.a[i], tensor = self.H[i])

          if self.b[i] is None:
            self.b[i] = dl.assemble(self.L[i])
          else:
            dl.assemble(self.L[i], tensor = self.b[i])
        
        # step 2 - solve
        for i in range(self.nspecies):
          self.solver.solve(self.H[i], self.u[i].vector(), self.b[i])

        self.solve_vascular_fields()

        local_tstep_count = local_tstep_count + 1
        
      # check if current time is exactly equal to acquisition time
      if np.abs(self.time - tacq) > 1.e-10:
        self.print_msg('current time {} and acquisition time {} do not match'.format(self.time, tacq), 1)
      #self.time = tacq + time_diff
      
      # output qoi at current time (ie tacq + time_diff)
      out[:, iacq] = self.evaluate()
      self.print_msg('  QoI added: {}'.format(out[:, iacq]), 2)
      
      # save model state
      if saved_state_current_sim:
        self.print_msg('  Output time: {}, Output number: {}'.format(self.time, self.out_count + 1), 2)
        self.save_state(save_vasc = True)

    self.print_msg('Total output: {}'.format(self.out_count), 2)

    # increment the counter for number of times the forward model is solved
    self.sim_count += 1
    if saved_state_current_sim:
      self.sim_out_current += 1
    
    # save QoI
    if saved_state_current_sim and self.rank == 0:
      np.save(self.sim_model_out_path + 'qoi.npy', out)
        
    end = time.time()
    if self.rank == 0:
      self.print_msg('compute time {}s'.format((end-start)), 2)
          
          
  def solveFwd(self, out, x):
    """ 
       Takes in the specified parameter and solves the forward model and computes the QoI (model output). 
    """
    
    self.print_msg('Solve with parameters: {}, sim count: {}'.format(np.exp(x[PARAMETER].get_local()), self.sim_count+1), 1)
    
    # replace parameter with specified parameters
    self.set_parameters(x[PARAMETER]) 
    
    # solve
    self.solve(out)
