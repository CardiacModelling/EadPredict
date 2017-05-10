function results = LancasterSobieClassify()

% take in data from file
start = getenv('CHASTE_TEST_OUTPUT');
data = importdata(strcat(start,'Tox_Res_Paper/all_results.tsv'),'\t');
% get redfern categories and merge categories 1 and 2
redferns = data.data(:,10);
redferns(redferns==1)=2;
grandi = importdata(strcat(start,'Tox_Res_Paper/grandi_l_s.tsv'),'\t');
ohara = importdata(strcat(start,'Tox_Res_Paper/ohara_l_s.tsv'),'\t');

% intialise results matrix
grandi_results = zeros(size(grandi.data,1),1);
trainingdata = grandi.data;
% loop through each data point
for data_point = 1:size(trainingdata,1)
    test_point = trainingdata(data_point,:);
    trainingmatrix = [trainingdata(1:data_point-1,:); trainingdata(data_point+1:end,:)];
    groupsmatrix = [redferns(1:data_point-1); redferns(data_point+1:length(redferns))];
    % classify
    class = classify(test_point, trainingmatrix, groupsmatrix,'linear');
    grandi_results(data_point) = class;
end
%% find boundaries
[~, ~, ~, ~, coeffs] = classify([1 1], trainingdata,redferns, 'linear');
grandi_boundary23 = -coeffs(1,2).const/coeffs(1,2).linear;
grandi_boundary34 = -coeffs(2,3).const/coeffs(2,3).linear;
grandi_boundary45 = -coeffs(3,4).const/coeffs(3,4).linear;
grandi_boundaries= [grandi_boundary23; grandi_boundary34; grandi_boundary45];

%% output to file
dlmwrite(strcat(start,'grandi_results.tsv'),grandi_results,'\t')
dlmwrite(strcat(start,'grandi_boundaries.tsv'),grandi_boundaries,'\t')


% intialise results matrix
ohara_results = zeros(size(ohara.data,1),1);
trainingdata = ohara.data;
% loop through each data point
for data_point = 1:size(trainingdata,1)
    test_point = trainingdata(data_point,:);
    trainingmatrix = [trainingdata(1:data_point-1,:); trainingdata(data_point+1:end,:)];
    groupsmatrix = [redferns(1:data_point-1); redferns(data_point+1:length(redferns))];
    % classify
    class = classify(test_point, trainingmatrix, groupsmatrix,'linear');
    ohara_results(data_point) = class;
end
%% find boundaries
[~, ~, ~, ~, coeffs] = classify([1 1], trainingdata,redferns, 'linear');
ohara_boundary23 = -coeffs(1,2).const/coeffs(1,2).linear;
ohara_boundary34 = -coeffs(2,3).const/coeffs(2,3).linear;
ohara_boundary45 = -coeffs(3,4).const/coeffs(3,4).linear;
ohara_boundaries= [ohara_boundary23; ohara_boundary34; ohara_boundary45];

%% output to file
dlmwrite(strcat(start,'ohara_results.tsv'),ohara_results,'\t')
dlmwrite(strcat(start,'ohara_boundaries.tsv'),ohara_boundaries,'\t')