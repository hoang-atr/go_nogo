function event_timing = parse_event_timing(all_event_timings,Ntrial,trial_duration)
event_timing = cell(Ntrial,1);
for n = 1:Ntrial
    t0 = (n-1)*trial_duration;
    t1 = n*trial_duration;
    
    ix = (all_event_timings>=t0 & all_event_timings<t1);
    event_timing{n} = all_event_timings(ix)-t0;
end