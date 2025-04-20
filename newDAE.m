[~, XTrain, Ytraining_regression, XValidation, YValidation] = Link5_data_Generation(0.95, [0,30], 5000);
Layers = [
    imageInputLayer([128 1 1], "Name", "imageinput", "Normalization", "none")
% 第一层：减少神经元数量 + 添加 Dropout
    fullyConnectedLayer(64, "Name", "fc_1")        
    batchNormalizationLayer("Name", "bn_1")
    reluLayer("Name", "relu_1")
    dropoutLayer(0.4, "Name", "drop_1")            % 添加 40% 的 Dropout
    
    % 第二层：进一步缩减维度 + 分组全连接（Grouped FC）
    fullyConnectedLayer(32, "Name", "fc_2")        
    batchNormalizationLayer("Name", "bn_2")
    reluLayer("Name", "relu_2")
    dropoutLayer(0.4, "Name", "drop_2")            % 添加 40% 的Dropout
    
    % 第三层：轻量化瓶颈层
    fullyConnectedLayer(16, "Name", "fc_3")         
    batchNormalizationLayer("Name", "bn_3")
    reluLayer("Name", "relu_3")  

    % 第四层：  分组全连接（Grouped FC）
    fullyConnectedLayer(32, "Name", "fc_4")        
    batchNormalizationLayer("Name", "bn_4")
    reluLayer("Name", "relu_4")
    dropoutLayer(0.2, "Name", "drop_2") 

    % 第五层：  分组全连接（Grouped FC）
    fullyConnectedLayer(64, "Name", "fc_5")        
    batchNormalizationLayer("Name", "bn_5")
    reluLayer("Name", "relu_5")
    dropoutLayer(0.2, "Name", "drop_2") 

% 第四层：保持原始输出维度
    fullyConnectedLayer(128, "Name", "fc_6")
    regressionLayer("Name", "regressionoutput")
];

Options = trainingOptions('adam', ...  
    'MaxEpochs', 100, ...
    'MiniBatchSize', 64, ...
    'InitialLearnRate', 0.005, ...
    'LearnRateSchedule',  'piecewise', ...
    'LearnRateDropFactor', 0.7, ...
    'LearnRateDropPeriod', 8, ...
    'L2Regularization', 0.05, ...
    'ValidationData', {XValidation, YValidation}, ...
    'ValidationFrequency', 50, ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress');

[DNN_Trained, info] = trainNetwork(XTrain, Ytraining_regression, Layers, Options);
save('DAE_Link5.mat','DNN_Trained');