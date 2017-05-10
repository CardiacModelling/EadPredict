function PlotEADClassificationExamples()
%% Look at the test traces used for checking afterdepolarisation detection
% Input files: test/TestTraces/*
% Input MATLAB: none
% Output files: ~/Documents/
% Output MATLAB: graphs

prefix = '/home/scratch/workspace/BethM/test/TestTraces/*.dat';

single_with = 'ESW_tentusscher_model_2006_endomembrane_L_type_calcium_current_conductance17.dat';
multiple_without = 'EMO_shannon_wang_puglisi_weber_bers_2004_model_updatedmembrane_fast_sodium_current_shift_inactivation50.dat';

time_start = 1;
time_end = 20000;

subplot(1,3,1)
filename = ['/home/scratch/workspace/BethM/test/TestTraces/',single_with];
data = importdata(filename);
times = data.data(time_start:time_end,1);
voltages = data.data(time_start:time_end,2);
plot(times,voltages,'k','Linewidth',2)
ylabel('Voltage (mV)')
xlabel('Time (ms)')
xlim([0 2000])
title('Single with repolarisation')
tight(gca)


subplot(1,3,3)
filename = ['/home/scratch/workspace/BethM/test/TestTraces/',multiple_without];
data = importdata(filename);
times = data.data(time_start:time_end,1);
voltages = data.data(time_start:time_end,2);
plot(times,voltages,'k','Linewidth',2)
ylabel('Voltage (mV)')
xlabel('Time (ms)')
xlim([0 2000])
title('Multiple without repolarisation')
tight(gca)

%% extra one
time_start = 60001;
time_end = 80000;
subplot(1,3,2)
filename = '/home/scratch/testoutput/SingleRuns/ohara_rudy_2011membrane_L_type_calcium_current_conductance28/ohara_rudy_2011membrane_L_type_calcium_current_conductance28.dat';
data = importdata(filename);
times = data.data(time_start:time_end,1)-time_start*0.1;
voltages = data.data(time_start:time_end,2);
plot(times,voltages,'k','Linewidth',2)
ylabel('Voltage (mV)')
xlabel('Time (ms)')
graphtitle = 'Multiple with repolarisation';
title(graphtitle)
xlim([0 2000])
tight(gca)


tidyprint(20,6,'Tox_Res_Paper/Graphs/Examples/Trio')
end

function tight(ax)

outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1) - 0.01;
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4) - 0.1;
ax.Position = [left bottom ax_width ax_height];

end
