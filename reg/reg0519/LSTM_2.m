clear
  load('ACE4s.mat');
  load('RegUD.mat');
  load('solar5min0919.mat');
  for i=1:8640
      for j=1:5
        solar1min(5*(i-1)+j,:)=solar5min(i,:);
      end
  end
  ACE4s=T_reg_4sec{1:end-1,2};a=length(ACE4s);k_5_15=15;
  a1=reshape(ACE4s,k_5_15,a/k_5_15);ACE1min=sum(a1)'/k_5_15;   

  ACE4sday=reshape(ACE4s,length(ACE4s)/30,30);
  RegUp=RegPro{:,6};
  RegDown=RegPro{:,5};
  RegPro1=[RegUp RegDown];
%   for i=1:2880
%       for j=1:225
%         RegPro2(225*(i-1)+j,:)=RegPro1(i,:);
% %         RegPro2(3*(i-1)+2,:)=RegPro1(i,:);
% %         RegPro2(3*(i-1)+3,:)=RegPro1(i,:);
%       end
%     %RegPro2(4*(i-1)+4,:)=RegPro1(i,:);
%   end
  
%     k_5_15=75;%75 or 225
%     a=length(ACE4s);a1=reshape(ACE4s,k_5_15,a/k_5_15);ACE_STARmin=sum(a1)'/k_5_15;   
% %     a2=reshape(REUP4sec,k_5_15,a/k_5_15);REUPmin=sum(a2)'/k_5_15;   
% %     a3=reshape(REUD4sec,k_5_15,a/k_5_15);REUDmin=sum(a3)'/k_5_15;   
%     a4=reshape(ACE4s,k_5_15,a/k_5_15);ACEmin=sum(a4)'/k_5_15;  
% %solar=readtable('C:\Users\lxh180005\Desktop\solar ramp\1-min Resolution Quantiles\Sept,2019\IBM_raw_CA_Topaz_reordered.csv');
% %load('Reg1920.mat')
% %T_reg_4sec1920.ACEStar=T_reg_4sec1920{:,7}-T_reg_4sec1920{:,8}+T_reg_4sec1920{:,9}; %plot(T_reg_4sec1920{1:43200,10})
% %ACE_1min_1920=T_reg_4sec1920{:,7}-T_reg_4sec1920{:,7}+T_reg_4sec1920{:,7};

 data=ACE1min(1:21600)';
numTimeStepsTrain = floor(0.9*numel(data));

dataTrain = data(1:numTimeStepsTrain+1);
dataTest = data(numTimeStepsTrain+1:end);
mu = mean(dataTrain);
sig = std(dataTrain);

dataTrainStandardized = (dataTrain - mu) / sig;
XTrain = dataTrainStandardized(1:end-1);
YTrain = dataTrainStandardized(2:end);
numFeatures = 1;
numResponses = 1;
numHiddenUnits = 200;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
    regressionLayer];
options = trainingOptions('adam', ...
    'MaxEpochs',250, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.005, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','training-progress');
net = trainNetwork(XTrain,YTrain,layers,options);
dataTestStandardized = (dataTest - mu) / sig;
XTest = dataTestStandardized(1:end-1);
net = predictAndUpdateState(net,XTrain);
[net,YPred] = predictAndUpdateState(net,YTrain(end));

numTimeStepsTest = numel(XTest);
for i = 2:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');
end
YPred = sig*YPred + mu;
YTest = dataTest(2:end);
rmse = sqrt(mean((YPred-YTest).^2));
figure
plot(dataTrain(1:end-1))
hold on
idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
plot(idx,[data(numTimeStepsTrain) YPred],'.-')
hold off
xlabel("Month")
ylabel("Cases")
title("Forecast")
legend(["Observed" "Forecast"])
figure
subplot(2,1,1)
plot(YTest)
hold on
plot(YPred,'.-')
hold off
legend(["Observed" "Forecast"])
ylabel("Cases")
title("Forecast")

subplot(2,1,2)
stem(YPred - YTest)
xlabel("Month")
ylabel("Error")
title("RMSE = " + rmse)

net = resetState(net);
net = predictAndUpdateState(net,XTrain);
YPred = [];
numTimeStepsTest = numel(XTest);
for i = 1:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,XTest(:,i),'ExecutionEnvironment','cpu');
end
YPred = sig*YPred + mu;
rmse = sqrt(mean((YPred-YTest).^2));
figure
subplot(2,1,1)
plot(YTest)
hold on
plot(YPred,'.-')
hold off
legend(["Observed" "Predicted"])
ylabel("Cases")
title("Forecast with Updates")

subplot(2,1,2)
stem(YPred - YTest)
xlabel("Month")
ylabel("Error")
title("RMSE = " + rmse)
%end


