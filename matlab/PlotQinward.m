function PlotQinward()

start = getenv('CHASTE_TEST_OUTPUT');
data = importdata(strcat(start,'Tox_Res_Paper/collated_data.tsv'),'\t');
%redferns = data.data(:,1);
drugnames = data.textdata;
colormap('summer')
% remove trailing whitespace
drugnames = strtrim(drugnames);
concs = {'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25'};
colors = hsv(2);
for i=0:1
    data = importdata(strcat(start,'qinward_traces_ajmaline', num2str(i)));
    voltage_trace = data.data(:,1);
    late_sodium_trace = data.data(:,2);
    calcium_trace = data.data(:,3);

    plot(voltage_trace, 'color',colors(i+1,:));
    hold on
    xlabel('Time (ms)')
    ylabel('Voltage (mV)')
    title('Effect of ajmaline on APD')
    legend({'control','with ajmaline'})

    %    subplot(1,2,1)
%    hold on
%    plot(late_sodium_trace,'color',colors(i+1,:))
%    subplot(1,2,2)
%    hold on
%    plot(calcium_trace,'color',colors(i+1,:))
end
% legend(concs)
% title('calcium current (pA/uF)')
% xlabel('time')
% 
% subplot(1,2,1)
% title('late sodium current (pA/uF)')
% xlabel('time')
% 
tidyprint(40,20,'Tox_Res_Paper/Graphs/ajmaline_trace')

end