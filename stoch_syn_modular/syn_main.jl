# Main file that runs stochastic model of synapse.

using PDMP, Sundials, Plots, GR, Interpolations, StatsBase

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
ntrials = 100;
results_time,results_XC,results_XD = simEPSPs(tstim,tend,ntrials);

# Interpolate sim data onto regular time array for averaging
include("interp1.jl")
dt = 0.1e-3
tvec = collect(0:dt:tend); # New regular time array
nt = length(tvec)
results_XC_grid = Array{Any}(ntrials);
results_XD_grid = Array{Any}(ntrials);
for i = 1:ntrials
  results_XC_grid[i] = zeros(Float64,6,length(tvec));
  for j = 1:length(xc0)
    results_XC_grid[i][j,:] = interp1(results_time[i], results_XC[i][j,:], tvec)
  end
end

# reorder so that outer list is by variable rather than by trial (default)
results_XC_reordered_EPSP = convertSimResults(results_XC_grid);

# Get mean, std, quantiles of continuous variable dynamics
nC = length(xc0);
mean_results_XC_EPSP = Array{Any}(nC);
std_results_XC_EPSP = Array{Any}(nC);
quantiles_results_XC_EPSP = Array{Any}(nC);
for i = 1:nC
  mean_results_XC_EPSP[i] = mean(results_XC_reordered_EPSP[i],2);
  std_results_XC_EPSP[i] = std(results_XC_reordered_EPSP[i],2);

  quantiles_results_XC_EPSP[i] = zeros(nt,3);
  for t = 1:nt
    quantiles_results_XC_EPSP[i][t,:] = Base.quantile(results_XC_reordered_EPSP[i][t,:],[0.25 0.5 0.75]);
  end
end


##########
# POST ANALYSIS
#########

peak_Ca_EPSP = zeros(ntrials,1)
peak_B_EPSP = zeros(ntrials,1)
peak_Vspine_EPSP = zeros(ntrials,1)
peak_Vdend_EPSP = zeros(ntrials,1)
for i = 1:ntrials
  peak_Vspine_EPSP[i] = maximum(results_XC_reordered_EPSP[1][:,i])
  peak_Vdend_EPSP[i] = maximum(results_XC_reordered_EPSP[2][:,i])
  peak_Ca_EPSP[i] = maximum(results_XC_reordered_EPSP[3][:,i])
  peak_B_EPSP[i] = maximum(results_XC_reordered_EPSP[4][:,i])
end

mean_peak_Vspine_EPSP = mean(peak_Vspine_EPSP - Eleak)
mean_peak_Vdend_EPSP = mean(peak_Vdend_EPSP - Eleak)
mean_peak_Ca_EPSP = mean(peak_Ca_EPSP)
mean_peak_B_EPSP = mean(peak_B_EPSP)

cv_peak_Vspine_EPSP = std(peak_Vspine_EPSP)/mean_peak_Vspine_EPSP
cv_peak_Vdend_EPSP = std(peak_Vdend_EPSP)/mean_peak_Vdend_EPSP
cv_peak_Ca_EPSP = std(peak_Ca_EPSP)/mean_peak_Ca_EPSP
cv_peak_B_EPSP = std(peak_B_EPSP)/mean_peak_B_EPSP

#########
# DENDRITIC CURRENT INJECTION EXPERIMENT
#########
tstim = 5e-3;
tend = 50e-3;
ntrials = 100;
A_bap = 40e-12;
tau_rise_bap = 0.1e-3;
tau_decay_bap = 1e-3;
results_time,results_XC,results_XD = simIdend(tstim,tend,ntrials,A_bap,tau_rise_bap,tau_decay_bap);

results_XC_grid = Array{Any}(ntrials);
results_XD_grid = Array{Any}(ntrials);
for i = 1:ntrials
  results_XC_grid[i] = zeros(Float64,6,length(tvec));
  for j = 1:length(xc0)
    results_XC_grid[i][j,:] = interp1(results_time[i], results_XC[i][j,:], tvec)
  end
end

# reorder so that outer list is by variable rather than by trial (default)
results_XC_reordered_Idend = convertSimResults(results_XC_grid);

# Get mean, std, quantiles of continuous variable dynamics
nC = length(xc0);
mean_results_XC_Idend = Array{Any}(nC);
std_results_XC_Idend = Array{Any}(nC);
quantiles_results_XC_Idend = Array{Any}(nC);
for i = 1:nC
  mean_results_XC_Idend[i] = mean(results_XC_reordered_Idend[i],2);
  std_results_XC_Idend[i] = std(results_XC_reordered_Idend[i],2);

  quantiles_results_XC_Idend[i] = zeros(nt,3);
  for t = 1:nt
    quantiles_results_XC_Idend[i][t,:] = Base.quantile(results_XC_reordered_Idend[i][t,:],[0.25 0.5 0.75]);
  end
end

# POST ANALYSIS

peak_Ca_Idend = zeros(ntrials,1)
peak_B_Idend = zeros(ntrials,1)
peak_Vspine_Idend = zeros(ntrials,1)
peak_Vdend_Idend = zeros(ntrials,1)
for i = 1:ntrials
  peak_Vspine_Idend[i] = maximum(results_XC_reordered_Idend[1][:,i])
  peak_Vdend_Idend[i] = maximum(results_XC_reordered_Idend[2][:,i])
  peak_Ca_Idend[i] = maximum(results_XC_reordered_Idend[3][:,i])
  peak_B_Idend[i] = maximum(results_XC_reordered_Idend[4][:,i])
end

mean_peak_Vspine_Idend = mean(peak_Vspine_Idend - Eleak)
mean_peak_Vdend_Idend = mean(peak_Vdend_Idend - Eleak)
mean_peak_Ca_Idend = mean(peak_Ca_Idend)
mean_peak_B_Idend = mean(peak_B_Idend)

cv_peak_Vspine_Idend = std(peak_Vspine_Idend)/mean_peak_Vspine_Idend
cv_peak_Vdend_Idend = std(peak_Vdend_Idend)/mean_peak_Vdend_Idend
cv_peak_Ca_Idend = std(peak_Ca_Idend)/mean_peak_Ca_Idend
cv_peak_B_Idend = std(peak_B_Idend)/mean_peak_B_Idend


#########
# PLOTTING
#########
gr(reuse=true)
# Plots.plot(tvec*1e3,mean_results_XC[4],label="mean")
Plots.plot(tvec*1e3,quantiles_results_XC_EPSP[4],label=["75%" "median" "25%"],color=:dodgerblue,line=[:dot :solid :dot],xlabel="Time (ms)",ylabel="Ca2+ dye")
Plots.plot!(tvec*1e3,quantiles_results_XC_Idend[4],label=["75%" "median" "25%"],color=:red,line=[:dot :solid :dot])

q = Plots.plot(tvec*1e3,results_XC_reordered_EPSP[4][:,1],color=:dodgerblue,label="EPSP",xlabel="Time (ms)",ylabel="Ca2+ dye")
Plots.plot!(tvec*1e3,results_XC_reordered_Idend[4][:,1],color=:red,label="Idend")
for i = 2:20
  Plots.plot!(tvec*1e3,results_XC_reordered_EPSP[4][:,i],color=:dodgerblue,label="")
  Plots.plot!(tvec*1e3,results_XC_reordered_Idend[4][:,i],color=:red,label="")
end
q

# Plot dendritic voltage
r = Plots.plot(tvec*1e3,1e3*results_XC_reordered_EPSP[2][:,1],color=:dodgerblue,label="EPSP",xlabel="Time (ms)",ylabel="Vdend (mV)")
Plots.plot!(tvec*1e3,1e3*results_XC_reordered_Idend[2][:,1],color=:red,label="Idend")
for i = 2:20
  Plots.plot!(tvec*1e3,1e3*results_XC_reordered_EPSP[2][:,i],color=:dodgerblue,label="")
  Plots.plot!(tvec*1e3,1e3*results_XC_reordered_Idend[2][:,i],color=:red,label="")
end
r

# Plot spine voltage
s = Plots.plot(tvec*1e3,1e3*results_XC_reordered_EPSP[1][:,1],color=:dodgerblue,label="EPSP",xlabel="Time (ms)",ylabel="Vspine (mV)")
Plots.plot!(tvec*1e3,1e3*results_XC_reordered_Idend[1][:,1],color=:red,label="Idend")
for i = 2:20
  Plots.plot!(tvec*1e3,1e3*results_XC_reordered_EPSP[1][:,i],color=:dodgerblue,label="")
  Plots.plot!(tvec*1e3,1e3*results_XC_reordered_Idend[1][:,i],color=:red,label="")
end
s

# Plots.plot(tvec*1e3,c_bap*(results_XC_reordered_Idend[6][:,1]-results_XC_reordered_Idend[5][:,1]),color=:dodgerblue,label="",xlabel="Time (ms)",ylabel="Idend")
