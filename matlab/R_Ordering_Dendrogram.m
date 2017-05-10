function R_Ordering_Dendrogram()
% use following collate_data.py, OutputCombinedMetrics.m, and
% OrderDendrogram.R

% take in data from file
start = getenv('CHASTE_TEST_OUTPUT');
data = importdata(strcat(start,'Tox_Res_Paper/collated_data.tsv'),'\t');
redferns = data.data(:,1);
column_headers = importdata(strcat(start,'Tox_Res_Paper/collated_data_key.dat'),'\t');
column_headers = strsplit(column_headers{1},'\t');
scaled_metrics = importdata(strcat(start,'Tox_Res_Paper/scaled_metrics.tsv'),'\t');
ordering_folder = '/home/scratch/workspace/BethM/test/MatlabAnalysis/R/';
ordering_filenames = {'order_APD90', 'order_INa_EADs', 'order_ICaL_EADs', 'order_IKr_EADs', 'order_APD90_EADs', 'order_EADs', 'order_APD90_ICaL', 'order_Grandi_LS','order_OHara_LS','order_hERG_cmax'};
graph_titles = {'APD90', 'INa EADs', 'ICaL EADs', 'IKr EADs', 'APD90 & EADs', 'EADs', 'APD90 & ICaL', 'Grandi LS','OHara LS','hERG/Cmax'};
combinations = {1, 6, 7, 8, [1 6 7 8], [6 7 8], [1 6], [2 3], [4 5]};
leaf_names = cell(length(data.textdata),1);
for drug = 1:length(data.textdata)
    leaf_names{drug} = [data.textdata{drug},' ',num2str(redferns(drug))];
end
optimal_order = sort(redferns);

for i=1:length(combinations)
    leaf_order = importdata(strcat(ordering_folder, ordering_filenames{i}));
    PlotDendrogram(scaled_metrics(:,combinations{i}),redferns,leaf_names,leaf_order,graph_titles{i});
    goodness = sum((optimal_order - redferns(leaf_order)).^2);
    sprintf('%s: %d',graph_titles{i},goodness)
end

%% hERG block / Cmax
leaf_order = importdata(strcat(ordering_folder,ordering_filenames{end}));
PlotDendrogram(data.data(:,end),redferns,leaf_names,leaf_order,graph_titles{end});
goodness = sum((optimal_order - redferns(leaf_order)).^2);
sprintf('%s: %d',graph_titles{end},goodness)

end

function PlotDendrogram(metric,redferns,leaf_names,leaf_order, graph_title)

tree = (pdist(metric));
links = linkage(tree);
figure;
[h,~,~] = log_dendrogram(links,length(redferns),'Orientation','left','Labels',leaf_names,'reorder',leaf_order);
title(graph_title)
set(h,'LineWidth',2)
set(h,'Color','k')
set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5])
set(gca, 'XMinorTick', 'on')
xlabel('Log Euclidean distance')
end