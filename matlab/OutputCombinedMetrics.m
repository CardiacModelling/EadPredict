function OutputCombinedMetrics()
start = getenv('CHASTE_TEST_OUTPUT');
data = importdata(strcat(start,'Tox_Res_Paper/collated_data.tsv'),'\t');
column_headers = importdata(strcat(start,'Tox_Res_Paper/collated_data_key.dat'),'\t');
column_headers = strsplit(column_headers{1},'\t');
redferns = data.data(:,1);

%% scale data to combine
controls = [286.466, 230.688, 8.82E-05, 208.083, 8.61E-05, 17.4121, 24.1268, 0.153015];

% make sure everything's where you expect it to be
expected_headers = {'Drug name','Redferns', 'APD90','Grandi L-S APD50',...
    'Grandi L-S Cai',    'OHara L-S APD50',    'OHara L-S Cai',    'INa shift',...
    'ICaL increase',    'IKr block',    'Herg block / Cmax'};
scaled_metrics = zeros(length(redferns),length(controls));
if sum(strcmp(column_headers,expected_headers)) == length(column_headers)
    
    for i=1:length(controls)
        scaled_metrics(:,i) = data.data(:,i+1) - controls(i);
        if abs(max(scaled_metrics(:,i))) > abs(min(scaled_metrics(:,i)))
            scaled_metrics(:,i) = scaled_metrics(:,i) ./ abs(max(scaled_metrics(:,i)));
        else
            scaled_metrics(:,i) = scaled_metrics(:,i) ./ abs(min(scaled_metrics(:,i)));
        end
    end
end
dlmwrite(strcat(start,'Tox_Res_Paper/','scaled_metrics.tsv'),scaled_metrics,'\t')
f = fopen(strcat(start,'Tox_Res_Paper/','scaled_metrics_key.tsv'),'w');
for i=1:length(controls)
   fprintf(f, '%s\t',column_headers{i+2});
end
fclose(f);
end
