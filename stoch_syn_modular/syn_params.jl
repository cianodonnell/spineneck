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

# Electrical properties of spine and dendrite
Cs = 1e-6*1e-8*1; # F; spine head capacitance
Cd = 1e-6*1e-8*1000; # F; dendrite capacitance
gdend = 1/(20e6); # Siemens; dendrite leak conductance
Eleak = -70e-3; # V; Leak conductance reversal potential
R_neck = 500.0e6 # Ohm

glu_base = 0 # glutamate pulse
glu_width = 2.0e-4 # (ms) glutamate pulse width

eta_ca_nmdar = 0.1e15 # scaling constant determining magnitude of calcium influx relative to NMDAR current
eta_ca_cav = 1e15 # scaling constant determining magnitude of calcium influx relative to CaV current
ca_tau = 0.014 # s, spine calcium decay time constant (pumping)

dye_kf = 800 # For OGB1 from Bartol et al 2015, assuming a [Ca] = 1uM.
dye_kb = 160 # For OGB1 from Bartol et al 2015
dye_tau = 1./(dye_kf + dye_kb)

A_bap = 400e-12 # Amps; Amplitude of dend current mimicking a BAP
tau_rise_bap = 0.1e-3 # s; rise time constant of BAP current
tau_decay_bap = 1e-3 # s; decay time constant of BAP current
tpeak = (tau_decay_bap*tau_rise_bap/(tau_decay_bap-tau_rise_bap))*log(tau_decay_bap/tau_rise_bap); # peak time of BAP current
c_bap = 1./(exp(-tpeak/tau_decay_bap) - exp(-tpeak/tau_rise_bap))# constant to make peak current equal 1

# saving saving in PDMP
sampling_rate = 10000. # it helps when total_rate is equal to zero...

# collect parameters
parms = Vector{Float64}([N_ampa, gamma_ampa, E_ampa, N_nmda, gamma_nmda, E_nmda, R_neck, Cs, Cd, gdend, Eleak, glu_base, glu_width, eta_ca_nmdar, ca_tau, Mg, eta_ca_cav, E_cav, gamma_car, gamma_cat, dye_kf, dye_kb, A_bap, tau_rise_bap, tau_decay_bap, sampling_rate])
