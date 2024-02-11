
data=X;
% data_corelation=[data ;Newntl((startday-1)*24+1:24*d,:)'];
% corrplot(data_corelation(:,1:end-72)')
[row,col] =size(data);


        numTimeStepsTrain = floor(0.9*numel(data(1,:)))-1;
%        dataTrain = data(1:numTimeStepsTrain+1);

%        data=[data;K];
        dataTrain = data(:,1:numTimeStepsTrain+1); %TRAIN
        
        dataTest = data(:,numTimeStepsTrain+1:end); %TEST
%        dataTest_realwobtm = data_realwobtm(:,numTimeStepsTrain+1:end); %TEST
        
        numFeatures = row; %
        numResponses = row; %        
        mu=zeros(numFeatures,1);sig=zeros(numFeatures,1);
        for i=1:numFeatures
            mu(i) = mean(dataTrain(i,:));
            sig(i) = std(dataTrain(i,:));
        end
        dataTrainStandardized = (dataTrain - mu) ./ sig;
        XTrain = dataTrainStandardized(:,1:end-1);%remain1
        YTrain = dataTrainStandardized(:,2:end);% emit first one

        numHiddenUnits = 200;

        layers = [ ...
            sequenceInputLayer(numFeatures)
            lstmLayer(numHiddenUnits)
            fullyConnectedLayer(numResponses)
            regressionLayer];
        options = trainingOptions('adam', ...
            'MaxEpochs',150, ...
            'GradientThreshold',1, ...
            'InitialLearnRate',0.005, ...
            'LearnRateSchedule','piecewise', ...
            'LearnRateDropPeriod',125, ...
            'LearnRateDropFactor',0.2, ...
            'Verbose',0, ...
            'Plots','training-progress');
        net(j) = trainNetwork(XTrain,YTrain,layers,options);
%        load lstmnet0213
        dataTestStandardized = (dataTest - mu) ./ sig;
        XTest = dataTestStandardized(:,1:end-1);
        net(j) = predictAndUpdateState(net(j),XTrain);%update last point
        [net(j),YPred] = predictAndUpdateState(net(j),YTrain(:,end));

        numTimeStepsTest = numel(XTest(1,:));
        for i = 2:numTimeStepsTest
            [net(j),YPred(:,i)] = predictAndUpdateState(net(j),YPred(:,i-1),'ExecutionEnvironment','cpu');
        end
        YPred = sig.*YPred + mu;
        YPred_DA=YPred;
        YTest = dataTest(:,2:end);
%        YTest1 = dataTest_realwobtm(:,2:end);

        rmse1 = sqrt(mean((YPred(1,:)-YTest(1,:)).^2));
        nrmse1=rmse1/(max(YTest(1,:))-min(YTest(1,:)));
        
%         rmse11 = sqrt(mean((YPred(1,:)+data_pv(end-70:end)-data_realntl(end-70:end)).^2));
%         nrmse11=rmse11/(max(data_realntl(end-70:end))-min(data_realntl(end-70:end)));
        
        figure
        plot(dataTrain(1,1:end-1))
        hold on
        idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
        plot(idx,[data(1,numTimeStepsTrain) YPred(1,:)],'.-')
        hold off
        xlabel("time")
        ylabel("power kw")
        title("Forecast")
        legend(["Observed" "Forecast"])
        figure
        subplot(2,1,1)
        plot(YTest(1,:))
        hold on
        plot(YPred(1,:),'.-')
        hold off
        legend(["Observed" "Forecast"])
        ylabel("Cases")
        title("Forecast")
        hold on

        %resetnet
        net(j) = resetState(net(j));
        net(j) = predictAndUpdateState(net(j),XTrain);
        YPred = [];
        numTimeStepsTest = numel(XTest(1,:));
        XTest(1,:)=XTest(1,:).*(0.8+0.2*rand(numTimeStepsTest,1))';
        for i = 1:numTimeStepsTest
            [net(j),YPred(:,i)] = predictAndUpdateState(net(j),XTest(:,i),'ExecutionEnvironment','cpu');
        end
        YPred = sig.*YPred + mu;
        rmse2 = sqrt(mean((YPred(1,:)-YTest(1,:)).^2));
        nrmse2=rmse2/(max(YTest(1,:))-min(YTest(1,:)));
        
%         
%         rmse3 = sqrt(mean((YPred(1,:)-YTest1(1,:)).^2));
%         nrmse3=rmse2/(max(YTest1(1,:))-min(YTest1(1,:)));
        
%         rmse22 = sqrt(mean((YPred(1,:)+data_pv(end-70:end)-data_realntl(end-70:end)).^2));
%         nrmse22=rmse22/(max(data_realntl(end-70:end))-min(data_realntl(end-70:end)));
%%        

        subplot(2,1,2)
        plot(YTest(1,:))
        hold on
        plot(YPred(1,:),'.-')
        hold off
        legend(["Observed" "Predicted"])
        ylabel("Cases")
        title("Forecast with Updates")
        Results{j}={net(j),YPred,YPred_DA,YTest,nrmse1,nrmse2};   
        RMSE_all(j,:)=[nrmse1 nrmse2];



