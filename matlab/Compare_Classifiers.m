function [svm_results, lda_results, output_headers] = Compare_Classifiers()

%% take in data from file
start = getenv('CHASTE_TEST_OUTPUT');
data = importdata(strcat(start,'Tox_Res_Paper/collated_data.tsv'),'\t');
column_headers = importdata(strcat(start,'Tox_Res_Paper/collated_data_key.tsv'),'\t');
column_headers = strsplit(column_headers{1},'\t');
output_headers = {};

%% check that the data haven't changed
expected_headers = {'Drug name','Redferns', 'APD90','Grandi L-S APD50',...
    'Grandi L-S Cai',    'OHara L-S APD50',    'OHara L-S Cai',    'INa shift',...
    'ICaL increase',    'IKr block',    'Herg block / Cmax'};
if sum(strcmp(column_headers,expected_headers)) ~= length(column_headers)
    return
end

%% pool redferns
redferns = data.data(:,1);
redferns(redferns==1)=2;

%% intialise results matrix
lda_results = zeros(size(data.data,1),size(data.data,2));
svm_results = zeros(size(data.data,1),size(data.data,2));
results_counter = 1;

%% take L-S measures and combine into one matrix each
% data.data(:,3:4)
%
trainingdata = data.data(:,3:4);
% loop through each data point
for data_point = 1:length(trainingdata)
    test_point = trainingdata(data_point,:);
    trainingmatrix = [trainingdata(1:data_point-1,:); trainingdata(data_point+1:end,:)];
    groupsmatrix = [redferns(1:data_point-1); redferns(data_point+1:length(redferns))];
    % classify using LDA
    class = classify(test_point, trainingmatrix, groupsmatrix,'linear');
    lda_results(data_point,results_counter) = class;
    % classify using SVM
    svm = fitcecoc(trainingmatrix,groupsmatrix);
    svm_class = predict(svm,test_point);
    svm_results(data_point,results_counter) = svm_class;
end
results_counter = results_counter + 1;
% add header
output_headers = {output_headers{:},'Grandi Lancaster-Sobie'};

trainingdata = data.data(:,5:6);
% loop through each data point
for data_point = 1:length(trainingdata)
    test_point = trainingdata(data_point,:);
    trainingmatrix = [trainingdata(1:data_point-1,:); trainingdata(data_point+1:end,:)];
    groupsmatrix = [redferns(1:data_point-1); redferns(data_point+1:length(redferns))];
    % classify using LDA
    class = classify(test_point, trainingmatrix, groupsmatrix,'linear');
    lda_results(data_point,results_counter) = class;
    % classify using SVM
    svm = fitcecoc(trainingmatrix,groupsmatrix);
    svm_class = predict(svm,test_point);
    svm_results(data_point,results_counter) = svm_class;
end
results_counter = results_counter + 1;
% add header
output_headers = {output_headers{:},'OHara Lancaster-Sobie'};

%% do the others
for intervention = 2:(length(column_headers)-1)
    % ignore L-S measures
    if intervention >2 && intervention < 7
        continue
    end
    trainingdata = data.data(:,intervention);
    % loop through each data point
    for data_point = 1:length(trainingdata)
        test_point = trainingdata(data_point);
        trainingmatrix = [trainingdata(1:data_point-1); trainingdata(data_point+1:end)];
        groupsmatrix = [redferns(1:data_point-1); redferns(data_point+1:length(redferns))];
        % classify using LDA
        class = classify(test_point, trainingmatrix, groupsmatrix,'linear');
        lda_results(data_point,results_counter) = class;
        % classify using SVM
        svm = fitcecoc(trainingmatrix,groupsmatrix);
        svm_class = predict(svm,test_point);
        svm_results(data_point,results_counter) = svm_class;
        
    end
    output_headers = {output_headers{:},column_headers{intervention+1}};
    results_counter = results_counter + 1;
end

%% combinations

combos = {[7 8 9],[2 7 8 9],[2 8]};
for this_combo = 1:length(combos)
    trainingdata = data.data(:,combos{this_combo});
    % loop through each data point
    for data_point = 1:length(trainingdata)
        test_point = trainingdata(data_point,:);
        trainingmatrix = [trainingdata(1:data_point-1,:); trainingdata(data_point+1:end,:)];
        groupsmatrix = [redferns(1:data_point-1); redferns(data_point+1:length(redferns))];
        % classify using LDA
        class = classify(test_point, trainingmatrix, groupsmatrix,'linear');
        lda_results(data_point,results_counter) = class;
        % classify using SVM
        svm = fitcecoc(trainingmatrix,groupsmatrix);
        svm_class = predict(svm,test_point);
        svm_results(data_point,results_counter) = svm_class;
        
    end
    this_header = '';
    for column = 1:length(combos{this_combo})
        this_header = strcat(this_header,' ',column_headers{combos{this_combo}(column)+1},' + ');
    end
    % remove trailing '+'
    this_header = this_header(1:end-1);
    output_headers = {output_headers{:},this_header};
    results_counter = results_counter + 1;
end




%% create output files
lda_filename = strcat(start,'lda_results.tsv');
lda_file = fopen(lda_filename,'w+');
svm_filename = strcat(start,'svm_results.tsv');
svm_file = fopen(svm_filename,'w+');

%% output data
for drug=1:size(data.data,1)
    % drug name
    fprintf(lda_file, '%s', data.textdata{drug});
    fprintf(svm_file, '%s', data.textdata{drug});
    for metric=1:(length(output_headers))
        % results
        fprintf(lda_file, '\t%d', lda_results(drug,metric));
        fprintf(svm_file, '\t%d', svm_results(drug,metric));
    end
    fprintf(lda_file, '\n');
    fprintf(svm_file, '\n');
end
fclose(lda_file);
fclose(svm_file);

%% output column headers
key_filename = strcat(start,'lda_svm_key.tsv');
key_file = fopen(key_filename,'w+');
fprintf(key_file,'Drug name');
for header=1:length(output_headers)
   fprintf(key_file, '\t%s', output_headers{header});
end
fclose(key_file);
end