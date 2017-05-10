function Plot_All_Thresholds()

% take in data from file
start = getenv('CHASTE_TEST_OUTPUT');
data = importdata(strcat(start,'Tox_Res_Paper/all_results.tsv'),'\t');
classifications = importdata('classification_results.tsv','\t');
boundaries = importdata('classification_boundaries.tsv','\t');
protocol_names = {'OHara APD90', 'I_CaL EAD threshold OHara', 'I_CaL 2ko EAD threshold OHara','I_CaL three_quarters ko EAD threshold OHara','I_Kr EAD threshold OHara', ...
    'I_Kr 2ko EAD threshold OHara', 'I_Kr three_quarters ko EAD threshold OHara', 'Grandi APD90', 'hERG IC50 divided by EFTPC_max'};
protocol_ylabel = {'APD90 (ms)','I_{CaL} threshold (\times)','I_{CaL} threshold (\times)','I_{CaL} threshold (\times)','I_{Kr} threshold (\times)','I_{Kr} threshold (\times)', ...
    'I_{Kr} threshold (\times)','APD90 (ms)', 'hERG IC_{50} / EFTPC_{max}'};
protocol_titles = {'OHara APD90', 'I_{CaL} EAD threshold OHara', 'I_{CaL} 2ko EAD threshold OHara','I_{CaL} 3/4 ko EAD threshold OHara','I_{Kr} EAD threshold OHara', ...
    'I_{Kr} 2ko EAD threshold OHara', 'I_{Kr} 3/4 ko EAD threshold OHara', 'Grandi APD90', 'hERG IC_{50} / EFTPC_{max}'};

% get redfern categories and merge categories 1 and 2
redferns = data.data(:,10);
redferns(redferns==1)=2;
allnums = ones(length(redferns),1);

for protocol=1:9
    figure;
    plot(data.data(redferns==2,protocol),allnums(redferns==2),'x','Color','r');
    hold on;
    plot(data.data(redferns==3,protocol),allnums(redferns==3),'x','Color',[1,0.5,0.1]);
    plot(data.data(redferns==4,protocol),allnums(redferns==4),'yx');
    plot(data.data(redferns==5,protocol),allnums(redferns==5),'gx');
    legend({'Category 2 or 1','Category 3','Category 4','Category 5'})
    line([boundaries(1,protocol),boundaries(1,protocol)],[0 2],'Color',[1,0.5,0.1])
    line([boundaries(2,protocol),boundaries(2,protocol)],[0 2],'Color','y')
    line([boundaries(3,protocol),boundaries(3,protocol)],[0 2],'Color','g')
    
    xlabel(protocol_ylabel{protocol});
    title(protocol_titles{protocol});
    % print
    tidyprint(20,15,strcat('Tox_Res_Paper/Graphs/Boundaries_Classified/',protocol_names{protocol}))
    close all
end