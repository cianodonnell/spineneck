########
# Functions for simulating Harnett experiments
#######

function simEPSPs(tstim,tend,ntrials)
  # Event times
  event_times = [0, tstim, tend]
  event_indices = zeros(event_times)
  event_indices[2] = 1 # 0 = passively record states / 1 = EPSP / 2 = BAP

  # Run multiple realisations, save results
  results_time = Vector(ntrials)
  results_XC = Vector(ntrials)
  results_XD = Vector(ntrials)
  for countloop = 1:ntrials
    results_time[countloop],results_XC[countloop],results_XD[countloop] = stim_events(xc0,xd0,parms,event_times,event_indices);
  end

  return results_time, results_XC, results_XD

end

function convertSimResults(results) # restructure sim results so outer list is by variable rather than by trial
  ntrials = length(results);
  nvariables = length(results[1][:,1]);
  nt = length(results[1][1,:]);
  results_reordered = Array{Any}(nvariables);
  for i = 1:nvariables
      results_reordered[i] = zeros(nt,ntrials);
      for j = 1:ntrials
        results_reordered[i][:,j] = results[j][:,i];
      end
  end

  return results_reordered
end
