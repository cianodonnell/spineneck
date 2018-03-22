###########
# FUNCTIONS DESCRIBING SYNAPTIC DYNAMICS
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
    rate[1] = xd[1]*glu*2500
    rate[2] = xd[2]*1000
    rate[3] = xd[3]*glu*2500
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
    return       xd[1]*glu*2500 +
                 xd[2]*1000 +
                 xd[3]*glu*2500 +
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
      parms[12] = 0 # set glutamate off
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
      parms[12] = 0 # set glutamate off
      res =  PDMP.pdmp!(xc0,xd0,F_ss!,R_ss,nu,parms,ts,tevent,ode=:lsoda,n_jumps = 1000)
      XC = hcat(XC,copy(res.xc[:,2:end]))
      XD = hcat(XD,copy(res.xd[:,2:end]))
      tt = vcat(tt,vec(res.time[2:end]))
      xc0 = XC[:,end]
      xd0 = XD[:,end]
      xc0[5] = 1.0;
      xc0[6] = 1.0;
      ts = tevent

    elseif event_indices[countloop]==0 #
      parms[12] = 0 # set glutamate off
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
