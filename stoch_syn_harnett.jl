#=
Model of a single synapse connected to a dendrite with stochastic AMPARAs, NMDARs, VGCCs. Stimulate with glutamate pulse.
Emulate Harnett et al Nature 2012 experiments. Can either activate synapse via glutamate pulse, or inject a current into the dendrite. NMDArs are blocked so only source of calcium is VGCCs.
=#


# load libraries
using PDMP
using Sundials

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
N_nmda = 0 # number of NMDARs
gamma_nmda = 20.0e-12 # (S) NMDA single channel conductance
E_nmda = E_ampa
Mg = 1 # uM

# Calcium channels
E_cav = 45e-3 # mV
N_car = 3 # number of R-type calcium channels
gamma_car = 17e-12  # (S) R-type single channel conductance
N_cat = 3 # number of T-type calcium channels
gamma_cat = 10e-12  # (S) T-type single channel conductance

# Electrical properties of spine
Cs = 1e-6*1e-8*1; # F; spine head capacitance
Cd = 1e-6*1e-8*500; # F; dendrite capacitance
gdend = 1/(100e6); # Siemens; dendrite leak conductance
Eleak = -70e-3; # V; Leak conductance reversal potential
R_neck = 500.0e6 # Ohm

glu_base = 1e-4 # glutamate pulse
glu_width = 2.0e-4 # (ms) glutamate pulse width

eta_ca_nmdar = 0.1e15 # scaling constant determining magnitude of calcium influx relative to NMDAR current
eta_ca_cav = 1e15 # scaling constant determining magnitude of calcium influx relative to CaV current
ca_tau = 0.014 # s, spine calcium decay time constant (pumping)

dye_kf = 800 # For OGB1 from Bartol et al 2015, assuming a [Ca] = 1uM.
dye_kb = 160 # For OGB1 from Bartol et al 2015
dye_tau = 1./(dye_kf + dye_kb)

A_bap = 400e-12 # Amps; Amplitude of dend current mimicking a BAP
tau_rise_bap = 0.2e-3 # s; rise time constant of BAP current
tau_decay_bap = 2e-3 # s; decay time constant of BAP current
c_bap = # constant to make peak current equal 1

# saving saving in PDMP
sampling_rate = 10. # it helps when total_rate is equal to zero...

# collect parameters
parms = Vector{Float64}([N_ampa, gamma_ampa, E_ampa, N_nmda, gamma_nmda, E_nmda, R_neck, Cs, Cd, gdend, Eleak, glu_base, glu_width, eta_ca_nmdar, ca_tau, Mg, eta_ca_cav, E_cav, gamma_car, gamma_cat, dye_kf, dye_kb, A_bap, tau_rise_bap, tau_decay_bap, sampling_rate]);

###########
# DEFINE DYNAMICS
#########

# Continuous variables
function F_ss!(xcdot::Vector{Float64},xc::Vector{Float64}, xd::Array{Int64}, t::Float64, parms::Vector)
  # vector field used for the continuous variables Vs and Vd
  Vs = xc[1];
  Vd = xc[2];
  Ca = xc[3];
  B = xc[4];
  I_bap = A_bap*c_bap*(xc[6]-xc[5]);

  No_ampa = xd[2];
  No_nmda = xd[4];
  No_car  = xd[10];
  No_cat  = xd[14];

  gamma_ampa = parms[2];
  E_ampa     = parms[3];
  gamma_nmda = parms[5];
  E_nmda     = parms[6];
  Rneck      = parms[7];
  Cs         = parms[8];
  Cd         = parms[9];
  gdend      = parms[10];
  Eleak      = parms[11];

  eta_ca_nmdar = parms[14];
  ca_tau       = parms[15];
  eta_ca_cav = parms[17];
  E_cav = parms[18];
  gamma_car = parms[19];
  gamma_cat = parms[20];
  dye_kf = parms[21];
  dye_kb = parms[22];
  tau_decay_bap = parms[24];
  tau_decay_bap = parms[25];

  # Spine voltage
  xcdot[1] = ( -(Vs-Vd)/Rneck
               - No_ampa*gamma_ampa*(Vs-E_ampa)
               - No_nmda*gamma_nmda*(Vs-E_nmda)
               - No_car*gamma_car*(Vs-E_cav)
               - No_cat*gamma_cat*(Vs-E_cav) )/Cs;

  # Dendrite voltage
  xcdot[2] = ( -(Vd-Vs)/Rneck - gdend*(Vd-Eleak) +  I_bap)/Cd;

  # Spine calcium
  xcdot[3] = ( -Ca/ca_tau
            + eta_ca_nmdar*No_nmda*gamma_nmda*(E_nmda-Vs)
            + eta_ca_cav*No_car*gamma_car*(E_cav-Vs)
            + eta_ca_cav*No_cat*gamma_cat*(E_cav-Vs) );

  # Spine calicum dye/buffer
  xcdot[4] = -dye_kb*B + dye_kf*Ca;

  # Dendritic injected current mimicking BAP/EPSP waveform
  xcdot[5] = - xc[5]/tau_rise_bap
  xcdot[6] = - xc[6]/tau_decay_bap

end

# Discrete stochastic variables: Number of open ampa
function R_ss(rate, xc, xd, t, parms, sum_rate::Bool)
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
    rate[1] = xd[1]*glu*5000
    rate[2] = xd[2]*500
    rate[3] = xd[3]*glu*5000
    rate[4] = xd[4]*1000*exp(-2.847) # A
    rate[5] = xd[4]*1000*exp(-0.016*Vs - 2.91) # a1
    rate[6] = xd[5]*1000*exp(0.009*Vs + 1.22) # b1
    rate[7] = xd[5]*1000*exp(-0.693) # B1
    rate[8] = xd[4]*Mg*1000*exp(-0.045*Vs - 6.97) # a2
    rate[9] = xd[6]*1000*exp(0.017*Vs + 0.96) # b2
    rate[10] = xd[6]*1000*exp(-3.101) # B2
    rate[11] = xd[7]*alpha_m_r
    rate[12] = xd[8]*beta_m_r
    rate[13] = xd[7]*alpha_h_r
    rate[14] = xd[9]*beta_h_r
    rate[15] = xd[8]*alpha_h_r
    rate[16] = xd[10]*beta_h_r
    rate[17] = xd[9]*alpha_m_r
    rate[18] = xd[10]*beta_m_r
    rate[19] = xd[11]*alpha_m_t
    rate[20] = xd[12]*beta_m_t
    rate[21] = xd[11]*alpha_h_t
    rate[22] = xd[13]*beta_h_t
    rate[23] = xd[12]*alpha_h_t
    rate[24] = xd[14]*beta_h_t
    rate[25] = xd[13]*alpha_m_t
    rate[26] = xd[14]*beta_m_t
    rate[27] = parms[end]
    return 0.
  else
    return       xd[1]*glu*5000 +
                 xd[2]*500 +
                 xd[3]*glu*5000 +
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
                 xd[14]*beta_m_t + parms[end]
  end
end

# matrix of jumps for the discrete variables
nu =       [[-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]; # AMPA C -> O
            [1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0]; # AMPA O -> C
            [0 0 -1 1 0 0 0 0 0 0 0 0 0 0 0]; # NMDA C -> O
            [0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0]; # NMDA O -> C
            [0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0]; # NMDA O -> B1
            [0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0]; # NMDA B1 -> O
            [0 0 1 0 -1 0 0 0 0 0 0 0 0 0 0]; # NMDA B1 -> C
            [0 0 0 -1 0 1 0 0 0 0 0 0 0 0 0]; # NMDA O -> B2
            [0 0 0 1 0 -1 0 0 0 0 0 0 0 0 0]; # NMDA B2 -> O
            [0 0 1 0 0 -1 0 0 0 0 0 0 0 0 0]; # NMDA B2 -> C
            [0 0 0 0 0 0 -1 1 0 0 0 0 0 0 0]; # CaR m0h0 -> m1h0
            [0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0]; # CaR m1h0 -> m0h0
            [0 0 0 0 0 0 -1 0 1 0 0 0 0 0 0]; # CaR m0h0 -> m0h1
            [0 0 0 0 0 0 1 0 -1 0 0 0 0 0 0]; # CaR m0h1 -> m0h0
            [0 0 0 0 0 0 0 -1 0 1 0 0 0 0 0]; # CaR m1h0 -> O
            [0 0 0 0 0 0 0 1 0 -1 0 0 0 0 0]; # CaR O -> m1h0
            [0 0 0 0 0 0 0 0 -1 1 0 0 0 0 0]; # CaR m0h1 -> O
            [0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0]; # CaR O -> m0h1
            [0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0]; # CaT m0h0 -> m1h0
            [0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0]; # CaT m1h0 -> m0h0
            [0 0 0 0 0 0 0 0 0 0 -1 0 1 0 0]; # CaT m0h0 -> m0h1
            [0 0 0 0 0 0 0 0 0 0 1 0 -1 0 0]; # CaT m0h1 -> m0h0
            [0 0 0 0 0 0 0 0 0 0 0 -1 0 1 0]; # CaT m1h0 -> O
            [0 0 0 0 0 0 0 0 0 0 0 1 0 -1 0]; # CaT O -> m1h0
            [0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0]; # CaT m0h1 -> O
            [0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0];# CaT O -> m0h1
            [0 0 0 0 0 0 0 0 0 0 0 0 0  0 1]] # Poisson sampling
# function to run between events
function stim_events(xc0,xd0,parms,event_times, event_indices)
  XC = copy(xc0)
  XD = copy(xd0)
  ts = 0.0
  tt = copy(ts)
  glu_width = parms[13]

  for (countloop, tevent) in enumerate(event_times)
    if event_indices[countloop]==1 # Activate synapses
      parms[12] = 1e-4 # set glutamate off
      res =  PDMP.pdmp!(xc0,xd0,F_ss!,R_ss,nu,parms,ts,tevent,ode=:lsoda,n_jumps = 1000)
      XC = hcat(XC,copy(res.xc[:,2:end]))
      XD = hcat(XD,copy(res.xd[:,2:end]))
      tt = vcat(tt,vec(res.time[2:end]))
      xc0 = XC[:,end]
      xd0 = XD[:,end]

      parms[12] = 1.0 # set glutamate on
      tend = tevent + glu_width
      res =  PDMP.pdmp!(xc0,xd0,F_ss!,R_ss,nu,parms,tevent,tend,ode=:lsoda,n_jumps = 1000)
      XC = hcat(XC,copy(res.xc[:,2:end]))
      XD = hcat(XD,copy(res.xd[:,2:end]))
      tt = vcat(tt,vec(res.time[2:end]))
      xc0 = XC[:,end]
      xd0 = XD[:,end]
      ts = tend

    elseif event_indices[countloop]==2 # BAP
      parms[12] = 1e-4 # set glutamate off
      xc0[5] = 1.0;
      xc0[6] = 1.0;
      res =  PDMP.pdmp!(xc0,xd0,F_ss!,R_ss,nu,parms,ts,tevent,ode=:lsoda,n_jumps = 1000)
      XC = hcat(XC,copy(res.xc[:,2:end]))
      XD = hcat(XD,copy(res.xd[:,2:end]))
      tt = vcat(tt,vec(res.time[2:end]))
      xc0 = XC[:,end]
      xd0 = XD[:,end]
      ts = tevent

    elseif event_indices[countloop]==0 #
      parms[12] = 1e-4 # set glutamate off
      res =  PDMP.pdmp!(xc0,xd0,F_ss!,R_ss,nu,parms,ts,tevent,ode=:lsoda,n_jumps = 1000)
      XC = hcat(XC,copy(res.xc[:,2:end]))
      XD = hcat(XD,copy(res.xd[:,2:end]))
      tt = vcat(tt,vec(res.time[2:end]))
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
xc0 = vec([-70e-3, -70e-3, 0, 0, 0]) # Initial states of continuous variables
xd0 = vec([N_ampa, 0, N_nmda, 0, 0, 0, 0, 0, N_car, 0, 0, 0, N_cat, 0, 0]); # Initial states of discrete variables

# parameters
tf = 50 # ms

# Event times
event_times = 1e-3*collect(5:5:tf)
#event_times = [10e-3 100e-3]
event_indices = zeros(UInt8,event_times)
event_indices[1] = 2 # 0 = passively record states / 1 = EPSP / 2 = BAP

########
# RUN PROGRAM
########

# compile the program:
dummy =  PDMP.pdmp!(xc0,xd0,F_ss!,R_ss,nu,parms,0.0,tf*1e-3,n_jumps = 100)

# compute a trajectory
tt,XC,XD = stim_events(xc0,xd0,parms,event_times,event_indices);

ntrials = 50
results_time = Vector(ntrials)
results_XC = Vector(ntrials)
results_XD = Vector(ntrials)
for countloop = 1:ntrials
  results_time[countloop],results_XC[countloop],results_XD[countloop] = stim_events(xc0,xd0,parms,event_times,event_indices);
end

using Plots, GR
gr(reuse=false)
# plotly()

##########
# POST ANALYSIS
#########

peak_Ca = zeros(ntrials,1)
peak_B = zeros(ntrials,1)
peak_Vspine = zeros(ntrials,1)
peak_Vdend = zeros(ntrials,1)
for i = 1:ntrials
  peak_Vspine[i] = maximum(results_XC[i][1,:])
  peak_Vdend[i] = maximum(results_XC[i][2,:])
  peak_Ca[i] = maximum(results_XC[i][3,:])
  peak_B[i] = maximum(results_XC[i][4,:])
end

mean_peak_Vspine = mean(peak_Vspine - Eleak);
mean_peak_Vdend = mean(peak_Vdend - Eleak);
mean_peak_B = mean(peak_B);

cv_peak_Vspine = std(peak_Vspine)/mean_peak_Vspine
cv_peak_Vdend = std(peak_Vdend)/(mean(peak_Vdend - Eleak))
cv_peak_Ca = std(peak_Ca)/mean(peak_Ca)
cv_peak_B = std(peak_B)/mean(peak_B)



#########
# PLOTTING
#########

# Plots.plot(tt*1000, XD[2,:],line=:step,xlabel = "Time (ms)",ylabel = "Number open",label="AMPA")
# Plots.plot!(tt*1000, XD[4,:],line=:step,xlabel = "Time (ms)",ylabel = "Number open",label="NMDA")
# Plots.plot(tt*1000, XC[1,:],label="Vsp")
# Plots.plot!(tt*1000, XC[2,:],label="Vdend")
# Plots.plot(tt*1000,XC[3,:],label="Spine [Ca2+]")

p = Plots.plot(results_time[1]*1000, results_XD[1][2,:],line=:step,label="AMPAR",xlabel="time (ms)",ylabel="Number of channels", linecolor=:dodgerblue)
Plots.display(Plots.plot!(results_time[1]*1000, results_XD[1][4,:],line=:step,label="NMDAR",linecolor=:red))
Plots.display(Plots.plot!(results_time[1]*1000, results_XD[1][10,:],line=:step,label="R-type",linecolor=:magenta))
Plots.display(Plots.plot!(results_time[1]*1000, results_XD[1][14,:],line=:step,label="T-type",linecolor=:green))
for i = 2:ntrials
  Plots.plot!(results_time[i]*1000, results_XD[i][2,:],line=:step,label="",linecolor=:blue)
  Plots.display(Plots.plot!(results_time[i]*1000, results_XD[i][4,:],line=:step,label="",linecolor=:red))
  Plots.display(Plots.plot!(results_time[i]*1000, results_XD[i][10,:],line=:step,label="",linecolor=:magenta))
  Plots.display(Plots.plot!(results_time[i]*1000, results_XD[i][14,:],line=:step,label="",linecolor=:green))
end
p


# Plot Calcium
q = Plots.plot(results_time[1]*1000, results_XC[1][3,:], label="Ca2+",xlabel="time (ms)",ylabel="a.u.", linecolor=:green,reuse=false)
for i = 2:ntrials
  Plots.plot!(results_time[i]*1000, results_XC[i][3,:],label="",linecolor=:green)
end
q

# Plot Buffer/dye
q2 = Plots.plot(results_time[1]*1000, results_XC[1][4,:], label="Dye",xlabel="time (ms)",ylabel="a.u.", linecolor=:green,reuse=false)
for i = 2:ntrials
  Plots.plot!(results_time[i]*1000, results_XC[i][4,:],label="",linecolor=:green)
end
q2

# Plot voltage
r = Plots.plot(results_time[1]*1000, 1e3*results_XC[1][1,:], label="Vspine",xlabel="time (ms)",ylabel="Voltage (mV)", linecolor=:dodgerblue)
Plots.plot!(results_time[1]*1000, 1e3*results_XC[1][2,:], label="Vdend", linecolor=:red)
for i = 2:ntrials
  Plots.plot!(results_time[i]*1000, 1e3*results_XC[i][1,:], label="", linecolor=:dodgerblue)
  Plots.plot!(results_time[i]*1000, 1e3*results_XC[i][2,:], label="", linecolor=:red)
end
r


Plots.scatter(1e3*peak_Vspine,1e3*peak_Vdend,xlim=(-70,-40),ylim=(-70,-65),xlabel="Peak Vspine (mV)",ylabel="Peak Vdend (mV)")
Plots.scatter(1e3*peak_Vspine,peak_Ca,xlim=(-70,-40),xlabel=" Peak Vspine (mV)",ylabel="Peak Calcium")
Plots.scatter(peak_Ca,peak_B,xlabel=" Peak Ca",ylabel="Peak Dye")


Plots.plot(results_time[1]*1000, results_XD[1][10,:],line=:step,label="",linecolor=:magenta)
Plots.display(Plots.plot!(results_time[1]*1000, results_XD[1][14,:],line=:step,label="",linecolor=:green))
for i = 2:ntrials
  Plots.display(Plots.plot!(results_time[i]*1000, results_XD[i][10,:],line=:step,label="",linecolor=:magenta))
  Plots.display(Plots.plot!(results_time[i]*1000, results_XD[i][14,:],line=:step,label="",linecolor=:green))
end
