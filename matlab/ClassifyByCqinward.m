function ClassifyByCqinward()

start = getenv('CHASTE_TEST_OUTPUT');
data = importdata('redferns.dat');
redferns = data.data(:,1);
drugnames = data.textdata;
% remove trailing whitespace
drugnames = strtrim(drugnames);
%colors = {[1 0 0],[1 0.5 0],[1 1 0],[0 1 0]};
lda_results = zeros(length(drugnames),26);
svm_results = zeros(length(drugnames),26);


for j = 1:25
    trainingdata = zeros(length(drugnames),1);
    %% get cqinward results
    for i=1:length(drugnames)
        cqinwarddata = importdata(strcat(start,'cqinward_results_',drugnames{i},'.tsv'));
        cqinwards = cqinwarddata.data(end-j+1,4);
        trainingdata(i) = cqinwards;
    end
    
    %% classify
    for data_point = 1:length(redferns)
        test_point = trainingdata(data_point,:);
        trainingmatrix = [trainingdata(1:data_point-1,:); trainingdata(data_point+1:end,:)];
        groupsmatrix = [redferns(1:data_point-1); redferns(data_point+1:length(redferns))];
        % classify using LDA
        class = classify(test_point, trainingmatrix, groupsmatrix,'linear');
        lda_results(data_point,j) = class;
        % classify using SVM
        svm = fitcecoc(trainingmatrix,groupsmatrix);
        svm_class = predict(svm,test_point);
        svm_results(data_point,j) = svm_class;
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
for data_point = 1:length(redferns)
    test_point = trainingdata(data_point,:);
    trainingmatrix = [trainingdata(1:data_point-1,:); trainingdata(data_point+1:end,:)];
    groupsmatrix = [redferns(1:data_point-1); redferns(data_point+1:length(redferns))];
    % classify using LDA
    class = classify(test_point, trainingmatrix, groupsmatrix,'linear');
    lda_results(data_point,j) = class;
    % classify using SVM
    svm = fitcecoc(trainingmatrix,groupsmatrix);
    svm_class = predict(svm,test_point);
    svm_results(data_point,j) = svm_class;
end



dlmwrite(strcat(start,'cqinward_classification_lda_results.tsv'),lda_results,'\t')
dlmwrite(strcat(start,'cqinward_classification_svm_results.tsv'),svm_results,'\t')

end