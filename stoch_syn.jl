
# load libraries
using PDMP,Sundials

######
# PARAMETERS
######

const N_ampa = 100 # number of AMPARs
gamma_ampa = 20.0e-12 # (S) AMPA single channel conductance
E_ampa = 0.0 # mV
N_nmda = 20 # number of NMDARs
gamma_nmda = 20.0e-12 # (S) NMDA single channel conductance
E_nmda = E_ampa

Cs = 1e-6*1e-8*1; # F; spine head capacitance
Cd = 1e-6*1e-8*500; # F; dendrite capacitance
gdend = 1/(100e6); # Siemens; dendrite leak conductance
Eleak = -70e-3; # V; Leak conductance reversal potential
R_neck = 500.0e6 # Ohm

# collect parameters
parms = Vector{Float64}([N_ampa, gamma_ampa, E_ampa, N_nmda, gamma_nmda, E_nmda, R_neck, Cs, Cd, gdend, Eleak]);

###########
# DEFINE DYNAMICS
#########

# Continuous variables
function F_ss!(xcdot::Vector{Float64},xc::Vector{Float64}, xd::Array{Int64}, t::Float64, parms::Vector)
  # vector field used for the continuous variables Vs and Vd
  Vs = xc[1];
  Vd = xc[2];

  No_ampa = xd[2];

  gamma_ampa = parms[2];
  E_ampa = parms[3];
  Rneck = parms[7];
  Cs = parms[8];
  Cd = parms[9];
  gdend = parms[10];
  Eleak = parms[11];

  xcdot[1] = ( -(Vs-Vd)/Rneck - No_ampa*gamma_ampa*(Vs-E_ampa) )/Cs;
  # xcdot[1] = ( -(Vs-Vd)/Rneck )/Cs;
  xcdot[2] = ( -(Vd-Vs)/Rneck - gdend*(Vd-Eleak) )/Cd;

end

# Discrete stochastic variables: Number of open ampa
function R_ss(xc, xd, t, parms, sum_rate::Bool)
  # rate function for each transition
  # in this case,  the transitions are xd->xd+2 or xd->xd-2
  N_ampa = parms[1];
  if sum_rate==false
    return vec([ xd[1]*100, xd[2]*100 ])
  else
    return vec([10,10])
  end
end

# matrix of jumps for the discrete variables, analogous to chemical reactions
const nu = [[-1 1];[1 -1]]

#######
# INITIAL CONDITIONS
######
xc0 = vec([-70e-3, -70e-3])
xd0 = vec([N_ampa, 0]); # Initial number of AMPArs closed, open



# parameters
tf = 500e-3

########
# RUN PROGRAM
########

# compile the program:
dummy =  PDMP.pdmp(100,xc0,xd0,F_ss!,R_ss,nu,parms,0.0,tf,false)

# compute a trajectory
result =  @time PDMP.pdmp(100,xc0,xd0,F_ss!,R_ss,nu,parms,0.0,tf,false)

using Plots, GR
gr()
Plots.plot(result.time, result.xd[2,:],line=:step,xlabel = "Time (s)",ylabel = "Number AMPARs",label="Xd")
Plots.plot(result.time, result.xc[1,:],label="Vsp")
Plots.plot!(result.time, result.xc[2,:],label="Vdend")
