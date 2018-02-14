
# load libraries
using PDMP,Sundials

# Continuous variables:
# 1) Spine voltage
# 2) Dendrite voltage
# 3) Spine calcium concentration
#
# Discrete variables
# 1) N_ampa_closed
# 2) N_ampa_open
# 3) N_nmda_closed
# 4) N_nmda_open
# 5) N_nmda_B1
# 6) N_nmda_B2
# 7) Car_m0h0
# 8) Car_m1h0
# 9) Car_m0h1
# 10) Car_open
# 11) Cat_m0h0
# 12) Cat_m1h0
# 13) Cat_m0h1
# 14) Cat_open

# NMDAR model from Jahr and Stevens J Neurosci (1990)


######
# PARAMETERS
######

# AMPARs
N_ampa = 100 # number of AMPARs
gamma_ampa = 20.0e-12 # (S) AMPA single channel conductance
E_ampa = 0.0 # mV

# NMDARs
N_nmda = 20 # number of NMDARs
gamma_nmda = 20.0e-12 # (S) NMDA single channel conductance
E_nmda = E_ampa
Mg = 1 # uM

# Calcium channels
E_cav = 10e-3 # mV
N_car = 5 # number of R-type calcium channels
gamma_car = 17e-12  # (S) R-type single channel conductance
N_cat = 5 # number of T-type calcium channels
gamma_cat = 10e-12  # (S) T-type single channel conductance

# Electrical properties of spine
Cs = 1e-6*1e-8*1; # F; spine head capacitance
Cd = 1e-6*1e-8*500; # F; dendrite capacitance
gdend = 1/(100e6); # Siemens; dendrite leak conductance
Eleak = -70e-3; # V; Leak conductance reversal potential
R_neck = 500.0e6 # Ohm

glu_base = 1e-4 # glutamate pulse
glu_width = 2.0e-4 # (ms) glutamate pulse width

eta_ca_nmdar = 0.1 # scaling constant determining magnitude of calcium influx relative to NMDAR current
eta_ca_cav = 1 # scaling constant determining magnitude of calcium influx relative to CaV current
ca_tau = 0.014 # s, spine calcium decay time constant (pumping)

# collect parameters
parms = Vector{Float64}([N_ampa, gamma_ampa, E_ampa, N_nmda, gamma_nmda, E_nmda, R_neck, Cs, Cd, gdend, Eleak, glu_base, glu_width, eta_ca_nmdar, ca_tau, Mg]);

###########
# DEFINE DYNAMICS
#########

# Continuous variables
function F_ss!(xcdot::Vector{Float64},xc::Vector{Float64}, xd::Array{Int64}, t::Float64, parms::Vector)
  # vector field used for the continuous variables Vs and Vd
  Vs = xc[1];
  Vd = xc[2];
  Ca = xc[3];

  No_ampa = xd[2];
  No_nmda = xd[4];
  No_car = xd[10];
  No_cat = xd[14];

  gamma_ampa = parms[2];
  E_ampa = parms[3];
  gamma_nmda = parms[5];
  E_nmda = parms[6];
  Rneck = parms[7];
  Cs = parms[8];
  Cd = parms[9];
  gdend = parms[10];
  Eleak = parms[11];

  eta_ca_nmdar = parms[14];
  ca_tau = parms[15];

  xcdot[1] = ( -(Vs-Vd)/Rneck
               - No_ampa*gamma_ampa*(Vs-E_ampa)
               - No_nmda*gamma_nmda*(Vs-E_nmda)
               - No_car*gamma_car*(Vs-E_cav)
               - No_cat*gamma_cat*(Vs-E_cav) )/Cs;
  # xcdot[1] = ( -(Vs-Vd)/Rneck )/Cs;
  xcdot[2] = ( -(Vd-Vs)/Rneck - gdend*(Vd-Eleak) )/Cd;
  xcdot[3] = -Ca/ca_tau
            + eta_ca_nmdar*No_nmda*gamma_nmda*(E_nmda-Vs)
            + eta_ca_cav*No_car*gamma_car*(E_cav-Vs)
            + eta_ca_cav*No_cat*gamma_cat*(E_cav-Vs);

end

# Discrete stochastic variables: Number of open ampa
function R_ss(xc, xd, t, parms, sum_rate::Bool)
  # rate function for each transition
  Vs = xc[1];

  N_ampa = parms[1];
  glu = parms[12] # Glutamate present? 0 = no, 1 = yes
  Mg = parms[16];

  beta_m_r_star = 1/(4e-4)
  minf_m_r_star = 1/(1+exp((0.003-0.01)/0.008))
  alpha_m_r_star = beta_m_r_star*minf_m_r_star/(1-minf_m_r_star)
  tau_m_r = 1/(alpha_m_r_star + beta_m_r_star)
  minf_r = 1/(1+exp((0.003-Vs)/0.008))
  alpha_m_r = minf_r/tau_m_r
  beta_m_r = (1-minf_r)/tau_m_r

  tau_h_r = 100e-3
  hinf_r = 1/(1 + exp((Vs+0.039)/0.0092))
  alpha_h_r = hinf_r/tau_h_r
  beta_h_r = (1-hinf_r)/tau_h_r

# Add in correct T-type rates
  beta_m_t_star = 1/(1e-3)
  minf_m_t_star = 1/(1+exp((-0.032+0.02)/0.007))
  alpha_m_t_star = beta_m_t_star*minf_m_t_star/(1-minf_m_t_star)
  tau_m_t = 1/(alpha_m_t_star + beta_m_t_star)
  minf_t = 1/(1+exp((-0.032-Vs)/0.007))
  alpha_m_t = minf_t/tau_m_t
  beta_m_t = (1-minf_t)/tau_m_t

  tau_h_t = 50e-3
  hinf_t = 1/(1 + exp((Vs+0.07)/0.0065))
  alpha_h_t = hinf_t/tau_h_t
  beta_h_t = (1-hinf_t)/tau_h_t


  if sum_rate==false
    return vec([ xd[1]*glu*3000,
                 xd[2]*500,
                 xd[3]*glu*3000,
                 xd[4]*1000*exp(-2.847), # A
                 xd[4]*1000*exp(-0.016*Vs - 2.91), # a1
                 xd[5]*1000*exp(0.009*Vs + 1.22), # b1
                 xd[5]*1000*exp(-0.693), # B1
                 xd[4]*Mg*1000*exp(-0.045*Vs - 6.97), # a2
                 xd[6]*1000*exp(0.017*Vs + 0.96), # b2
                 xd[6]*1000*exp(-3.101), # B2
                 xd[7]*alpha_m_r,
                 xd[8]*beta_m_r,
                 xd[7]*alpha_h_r,
                 xd[9]*beta_h_r,
                 xd[8]*alpha_h_r,
                 xd[10]*beta_h_r,
                 xd[9]*alpha_m_r,
                 xd[10]*beta_m_r,
                 xd[11]*alpha_m_t,
                 xd[12]*beta_m_t,
                 xd[11]*alpha_h_t,
                 xd[13]*beta_h_t,
                 xd[12]*alpha_h_t,
                 xd[14]*beta_h_t,
                 xd[13]*alpha_m_t,
                 xd[14]*beta_m_t])
  else
    return       xd[1]*glu*3000 +
                 xd[2]*500 +
                 xd[3]*glu*3000 +
                 xd[4]*1000*exp(-2.847) +
                 xd[4]*1000*exp(-0.016*Vs - 2.91) +
                 xd[5]*1000*exp(0.009*Vs + 1.22) +
                 xd[5]*1000*exp(-0.693) +
                 xd[4]*Mg*1000*exp(-0.045*Vs - 6.97) +
                 xd[6]*1000*exp(0.017*Vs + 0.96) +
                 xd[6]*1000*exp(-3.101) +
                 xd[7]*alpha_m_r +
                 xd[8]*beta_m_r +
                 xd[7]*alpha_h_r +
                 xd[9]*beta_h_r +
                 xd[8]*alpha_h_r +
                 xd[10]*beta_h_r +
                 xd[9]*alpha_m_r +
                 xd[10]*beta_m_r +
                 xd[11]*alpha_m_t +
                 xd[12]*beta_m_t +
                 xd[11]*alpha_h_t +
                 xd[13]*beta_h_t +
                 xd[12]*alpha_h_t +
                 xd[14]*beta_h_t +
                 xd[13]*alpha_m_t +
                 xd[14]*beta_m_t
  end
end

# matrix of jumps for the discrete variables
const nu = [[-1 1 0 0 0 0 0 0 0 0 0 0 0 0]; # AMPA C -> O
            [1 -1 0 0 0 0 0 0 0 0 0 0 0 0]; # AMPA O -> C
            [0 0 -1 1 0 0 0 0 0 0 0 0 0 0]; # NMDA C -> O
            [0 0 1 -1 0 0 0 0 0 0 0 0 0 0]; # NMDA O -> C
            [0 0 0 -1 1 0 0 0 0 0 0 0 0 0]; # NMDA O -> B1
            [0 0 0 1 -1 0 0 0 0 0 0 0 0 0]; # NMDA B1 -> O
            [0 0 1 0 -1 0 0 0 0 0 0 0 0 0]; # NMDA B1 -> C
            [0 0 0 -1 0 1 0 0 0 0 0 0 0 0]; # NMDA O -> B2
            [0 0 0 1 0 -1 0 0 0 0 0 0 0 0]; # NMDA B2 -> O
            [0 0 1 0 0 -1 0 0 0 0 0 0 0 0]; # NMDA B2 -> C
            [0 0 0 0 0 0 -1 1 0 0 0 0 0 0]; # CaR m0h0 -> m1h0
            [0 0 0 0 0 0 1 -1 0 0 0 0 0 0]; # CaR m1h0 -> m0h0
            [0 0 0 0 0 0 -1 0 1 0 0 0 0 0]; # CaR m0h0 -> m0h1
            [0 0 0 0 0 0 1 0 -1 0 0 0 0 0]; # CaR m0h1 -> m0h0
            [0 0 0 0 0 0 0 -1 0 1 0 0 0 0]; # CaR m1h0 -> O
            [0 0 0 0 0 0 0 1 0 -1 0 0 0 0]; # CaR O -> m1h0
            [0 0 0 0 0 0 0 0 -1 1 0 0 0 0]; # CaR m0h1 -> O
            [0 0 0 0 0 0 0 0 1 -1 0 0 0 0]; # CaR O -> m0h1
            [0 0 0 0 0 0 0 0 0 0 -1 1 0 0]; # CaT m0h0 -> m1h0
            [0 0 0 0 0 0 0 0 0 0 1 -1 0 0]; # CaT m1h0 -> m0h0
            [0 0 0 0 0 0 0 0 0 0 -1 0 1 0]; # CaT m0h0 -> m0h1
            [0 0 0 0 0 0 0 0 0 0 1 0 -1 0]; # CaT m0h1 -> m0h0
            [0 0 0 0 0 0 0 0 0 0 0 -1 0 1]; # CaT m1h0 -> O
            [0 0 0 0 0 0 0 0 0 0 0 1 0 -1]; # CaT O -> m1h0
            [0 0 0 0 0 0 0 0 0 0 0 0 -1 1]; # CaT m0h1 -> O
            [0 0 0 0 0 0 0 0 0 0 0 0 1 -1]] # CaT O -> m0h1

# function to run between events
function stim_events(xc0,xd0,parms,event_times, event_indices)
  XC = copy(xc0)
  XD = copy(xd0)
  ts = 0.0
  tt = copy(ts)
  glu_width = parms[13]

  for (countloop, tevent) in enumerate(event_times)
    if event_indices[countloop]==1
      parms[12] = 1e-4 # set glutamate off
      res =  PDMP.pdmp(100,xc0,xd0,F_ss!,R_ss,nu,parms,ts,tevent,false,ode=:cvode)
      XC=hcat(XC,copy(res.xc[:,2:end]))
      XD=hcat(XD,copy(res.xd[:,2:end]))
      tt=vcat(tt,vec(res.time[2:end]))
      xc0 = XC[:,end]
      xd0 = XD[:,end]

      parms[12] = 1.0 # set glutamate on
      tend = tevent + glu_width
      res =  PDMP.pdmp(100,xc0,xd0,F_ss!,R_ss,nu,parms,tevent,tend,false,ode=:cvode)
      XC=hcat(XC,copy(res.xc[:,2:end]))
      XD=hcat(XD,copy(res.xd[:,2:end]))
      tt=vcat(tt,vec(res.time[2:end]))
      xc0 = XC[:,end]
      xd0 = XD[:,end]
      ts = tend

    elseif event_indices[countloop]==0
      parms[12] = 1e-4 # set glutamate off
      res =  PDMP.pdmp(100,xc0,xd0,F_ss!,R_ss,nu,parms,ts,tevent,false,ode=:cvode)
      XC=hcat(XC,copy(res.xc[:,2:end]))
      XD=hcat(XD,copy(res.xd[:,2:end]))
      tt=vcat(tt,vec(res.time[2:end]))
      xc0 = XC[:,end]
      xd0 = XD[:,end]
      ts = tevent
    end
  end

  return tt,XC,XD

end

#######
# INITIAL CONDITIONS
######
xc0 = vec([-70e-3, -70e-3, 0])
xd0 = vec([N_ampa, 0, N_nmda, 0, 0, 0, N_car, 0, 0, 0, N_cat, 0, 0, 0]); # Initial states of AMPA/NMDArs

# parameters
tf = 100 # ms

# Event times
event_times = 1e-3*collect(5:5:tf)
event_times = [10e-3 100e-3]
event_indices = zeros(UInt8,event_times)
event_indices[1] = 1

########
# RUN PROGRAM
########

# compile the program:
dummy =  PDMP.pdmp(100,xc0,xd0,F_ss!,R_ss,nu,parms,0.0,tf*1e-3,false)

# compute a trajectory
tt,XC,XD = stim_events(xc0,xd0,parms,event_times,event_indices);

ntrials = 2
results_time = Vector(ntrials)
results_XC = Vector(ntrials)
results_XD = Vector(ntrials)
for countloop = 1:ntrials
  results_time[countloop],results_XC[countloop],results_XD[countloop] = stim_events(xc0,xd0,parms,event_times,event_indices);
end

using Plots, GR
gr(reuse=false)

# Plots.plot(tt*1000, XD[2,:],line=:step,xlabel = "Time (ms)",ylabel = "Number open",label="AMPA")
# Plots.plot!(tt*1000, XD[4,:],line=:step,xlabel = "Time (ms)",ylabel = "Number open",label="NMDA")
# Plots.plot(tt*1000, XC[1,:],label="Vsp")
# Plots.plot!(tt*1000, XC[2,:],label="Vdend")
# Plots.plot(tt*1000,XC[3,:],label="Spine [Ca2+]")

p = Plots.plot(results_time[1]*1000, results_XD[1][2,:],line=:step, linecolor=:blue)
for i = 1:ntrials
  Plots.plot!(results_time[i]*1000, results_XD[i][2,:],line=:step,linecolor=:blue)
  Plots.plot!(results_time[i]*1000, results_XD[i][4,:],line=:step,linecolor=:red)
end
p

gr(reuse=false)
q = Plots.plot(results_time[1]*1000, results_XC[1][3,:],line=:step, linecolor=:green)
for i = 2:ntrials
  Plots.plot!(results_time[i]*1000, results_XC[i][3,:],line=:step,linecolor=:green)
end
q
