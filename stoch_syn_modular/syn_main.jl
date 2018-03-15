# Main file that runs stochastic model of synapse.

using PDMP, Sundials, Plots, GR, Interpolations

include("syn_params.jl") # load model parameters
include("syn_dynamics.jl") # load functions describing model dynamics

#######
# SIMULATION PARAMETERS
######

# Initial conditions
xc0 = vec([-70e-3, -70e-3, 0, 0, 0, 0]) # Initial states of continuous variables
xd0 = vec([N_ampa, 0, N_nmda, 0, 0, 0, 0, 0, N_car, 0, 0, 0, N_cat, 0, 0]); # Initial states of discrete variables

########
# RUN PROGRAM
########
# compile the program
dummy =  PDMP.pdmp!(xc0,xd0,F_ss!,R_ss,nu,parms,0.0,50e-3,n_jumps = 100)


# Simulate EPSP experiment
include("harnett_exps.jl") # Load experiment simulation functions
tstim = 5e-3;
tend = 50e-3;
ntrials = 10;
results_time,results_XC,results_XD = simEPSPs(tstim,tend,ntrials);

# Interpolate sim data onto regular time array for averaging
include("interp1.jl")
dt = 0.1e-3
tvec = collect(0:dt:tend); # New regular time array
results_XC_grid = Array{Any}(ntrials);
results_XD_grid = Array{Any}(ntrials);
for i = 1:ntrials
  results_XC_grid[i] = zeros(Float64,6,length(tvec));
  for j = 1:length(xc0)
    results_XC_grid[i][j,:] = interp1(results_time[i], results_XC[i][j,:], tvec)
  end
end

results_XC_reordered = convertSimResults(results_XC_grid);

# Get mean and std of dynamics
nC = length(xc0);
mean_results_XC = Array{Any}(nC);
std_results_XC = Array{Any}(nC);
for i = 1:nC
  mean_results_XC[i] = mean(results_XC_reordered[i],2);
  std_results_XC[i] = std(results_XC_reordered[i],2);
end



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

mean_peak_Vspine = mean(peak_Vspine - Eleak)
mean_peak_Vdend = mean(peak_Vdend - Eleak)
mean_peak_B = mean(peak_B)

cv_peak_Vspine = std(peak_Vspine)/mean_peak_Vspine
cv_peak_Vdend = std(peak_Vdend)/(mean(peak_Vdend - Eleak))
cv_peak_Ca = std(peak_Ca)/mean(peak_Ca)
cv_peak_B = std(peak_B)/mean(peak_B)

#########
# PLOTTING
#########
