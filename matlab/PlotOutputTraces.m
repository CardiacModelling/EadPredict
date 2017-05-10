function PlotOutputTraces(drug_name)

start = getenv('CHASTE_TEST_OUTPUT');
conc = 1;

% undrugged trace
control_data = importdata(strcat(start,'drug_effect_traces_control'));
control_voltage_trace = control_data.data(:,1);

% times
times = [1:length(control_voltage_trace)].*0.1;

% drugged trace
data = importdata(strcat(start,'drug_effect_traces_', drug_name, num2str(conc)));
voltage_trace = data.data(:,1);

% plot
plot(times,control_voltage_trace)
hold on
plot(times,voltage_trace);
legend({'control',drug_name})
xlabel('Time (ms)')
ylabel('Membrane voltage (mV)')
title(drug_name)

tidyprint(10,10,drug_name)

%    late_sodium_trace = data.data(:,2);
%    calcium_trace = data.data(:,3);
%    herg_trace = solution.GetAnyVariable("membrane_rapid_delayed_rectifier_potassium_current");
%    ito_trace = solution.GetAnyVariable("membrane_transient_outward_current");
%    ik1_trace = solution.GetAnyVariable("membrane_inward_rectifier_potassium_current");
%    fast_sodium_trace = solution.GetAnyVariable("membrane_fast_sodium_current");
%    iks_trace = solution.GetAnyVariable("membrane_slow_delayed_rectifier_potassium_current_conductance");


end