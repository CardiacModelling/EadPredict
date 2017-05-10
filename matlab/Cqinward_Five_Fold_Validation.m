function Cqinward_Five_Fold_Validation()

% get redferns
start = getenv('CHASTE_TEST_OUTPUT');
data = importdata('redferns.dat');
redferns = data.data(:,1);
drugnames = data.textdata;

% remove trailing whitespace
drugnames = strtrim(drugnames);

% initalise vectors
lda_results = zeros(length(drugnames),26);
svm_results = zeros(length(drugnames),26);

% create groups
randomised_drugs = randperm(length(drugnames));
group_1 = randomised_drugs(1:7);
group_2 = randomised_drugs(8:14);
group_3 = randomised_drugs(15:23);
group_4 = randomised_drugs(24:32);
group_5 = randomised_drugs(33:41);
groups = {group_1,group_2,group_3,group_4,group_5};

% loop through concentrations
for j = 1:25
    %% get cqinward results
    trainingdata = zeros(length(drugnames),1);
    for i=1:length(drugnames)
        cqinwarddata = importdata(strcat(start,'cqinward_results_',drugnames{i},'.tsv'));
        cqinwards = cqinwarddata.data(end-j+1,4);
        trainingdata(i) = cqinwards;
    end
    
    for test_group = 1:length(groups)
        %% classify
        training_group = setdiff(1:length(trainingdata), groups{test_group});
        trainingmatrix = trainingdata(training_group);
        test_point = trainingdata(groups{test_group});
        groupsmatrix = redferns(training_group);
        % classify using LDA
        class = classify(test_point, trainingmatrix, groupsmatrix,'linear');
        for i=1:length(class)
            lda_results(groups{test_group}(i),j) = class(i);
        end
        % classify using SVM
        svm = fitcecoc(trainingmatrix,groupsmatrix);
        svm_class = predict(svm,test_point);
        for i=1:length(svm_class)
            svm_results(groups{test_group}(i),j) = svm_class(i);
        end
    end
end

j = 26;
%% combine em all together
trainingdata = zeros(length(drugnames),25);
%% get cqinward results
for i=1:length(drugnames)
    cqinwarddata = importdata(strcat(start,'cqinward_results_',drugnames{i},'.tsv'));
    cqinwards = cqinwarddata.data(2:end,4);
    trainingdata(i,:) = cqinwards;
end

%% classify
for test_group = 1:length(groups)
    %% classify
    training_group = setdiff([1:length(trainingdata)], groups{test_group});
    trainingmatrix = trainingdata(training_group,:);
    test_point = trainingdata(groups{test_group},:);
    groupsmatrix = redferns(training_group);
    % classify using LDA
    class = classify(test_point, trainingmatrix, groupsmatrix,'linear');
    for i=1:length(class)
        lda_results(groups{test_group}(i),j) = class(i);
    end
    % classify using SVM
    svm = fitcecoc(trainingmatrix,groupsmatrix);
    svm_class = predict(svm,test_point);
    for i=1:length(svm_class)
        svm_results(groups{test_group}(i),j) = svm_class(i);
    end
end

dlmwrite(strcat(start,'5_cqinward_classification_lda_results.tsv'),lda_results,'\t')
dlmwrite(strcat(start,'5_cqinward_classification_svm_results.tsv'),svm_results,'\t')

end