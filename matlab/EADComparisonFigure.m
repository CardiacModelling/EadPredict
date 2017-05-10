function EADComparisonFigure()
%% Make graphs of membrane voltage to illustrate differences in EAD thresholds for different drugs.
% Input files created by TestEADs.hpp
% Input files: testoutput/ADCheck/*

close all

controls = {'1_drug_effect_traces_cisapride_0_1','1_drug_effect_traces_cisapride_0_9','1_drug_effect_traces_cisapride_0_27','1_drug_effect_traces_cisapride_0_31'};
cisapride = {'1_drug_effect_traces_cisapride_1_1','1_drug_effect_traces_cisapride_1_9','1_drug_effect_traces_cisapride_1_27','1_drug_effect_traces_cisapride_1_31'};
nitrendipine = {'1_drug_effect_traces_nitrendipine_1_1','1_drug_effect_traces_nitrendipine_1_9','1_drug_effect_traces_nitrendipine_1_27','1_drug_effect_traces_nitrendipine_1_31'};
filenames = {nitrendipine;controls;cisapride};
drugnames = {'nitrendipine','controls','cisapride'};
change = {'g_{CaL} \times 1.00','g_{CaL} \times 9','g_{CaL} \times 25','g_{CaL} \times 28'};
start = getenv('CHASTE_TEST_OUTPUT');

for i=1:4
    subplot(4,4,4*(i-1)+1)
    axis off
    axis([0 2 0 2])
    text(1,1,change{i})
    for j=1:3
        if i + j == 5
            linecolor = 'r';
        else
            linecolor = 'k';
        end
        subplot(4,4,4*(i-1)+j+1)
        openfile = strcat(start,'Tox_Res_Paper/',filenames{j}{i});
        a = importdata(openfile);
        %% append the voltages to a big matrix
        half_size = floor(length(a.data(:,1))/2);
        offset = 15000;
        %quarter_size = floor(length(a.data(:,1))/4);
        voltages = a.data(1:half_size,1);
        times = (1:half_size)*0.1;
        %% make graph
        plot(times,voltages,linecolor,'LineWidth',2)
        %ylim([-100 50])
        %xlim([0 1000])
        set(gca,'box','off')
        axis off
        if i==1
            title(drugnames{j})
        end
    end
end


%% add scale bars
subplot(4,4,1)
hold on 
xoffset = 1.5;
text(xoffset,-0.4,'1 s')
text(xoffset-1,0.1,'50 mV')
plot([xoffset xoffset],[-0.1 0.4],'k','LineWidth',1,'Clipping','off')
plot([xoffset xoffset+1],[-0.1 -0.1],'k','LineWidth',1,'Clipping','off')

%% print to file
tidyprint(20,15,'Tox_Res_Paper/Graphs/EADComparisonFigure')
