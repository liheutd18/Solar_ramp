% %new data
 %T_reg_4sec1920 = readtable('C:\Users\lxh180005\Desktop\solar ramp\ACE_REG_2019_2020_received 2020_10_31_TRANSFORMED FOR JOSEPHINE.xlsx');
 load('Reg1920.mat');
 solardata = readtable('C:\Users\lxh180005\Desktop\solar ramp\OASIS Data\OASIS Data\SLD_REN_FCST_Feb_May_2020_Solar.csv', 'Delimiter',',');
 winddata = readtable('C:\Users\lxh180005\Desktop\solar ramp\OASIS Data\OASIS Data\SLD_REN_FCST_Feb_May_2020_Wind.csv', 'Delimiter',',');
 loaddata = readtable('C:\Users\lxh180005\Desktop\solar ramp\OASIS Data\OASIS Data\SLD_FCST_Feb_May_2020.csv', 'Delimiter',',');
 regupdata = readtable('C:\Users\lxh180005\Desktop\solar ramp\OASIS Data\20200501_20200531_AS_REQ_DAM_RegUp+Dn_20210529_03_10_20_v1_JosephineTransformed.xlsx','Sheet','RegUp');
 regup=table2array(regupdata(:,14));regup=regup';%may regulation procurement
 regup_11=repelem(regup,12,1);regup_11=regup_11';

% solardata.datetime = datetime(solardata.INTERVALSTARTTIME_GMT, 'InputFormat', 'yyyy-MMM-dd-HH:mm:ss.SSS');
% winddata.datetime = datetime(winddata.INTERVALSTARTTIME_GMT, 'InputFormat', 'yyyy-MMM-dd-HH:mm:ss.SSS');
% loaddata.datetime = datetime(loaddata.INTERVALSTARTTIME_GMT, 'InputFormat', 'yyyy-MMM-dd-HH:mm:ss.SSS');
% 
% idx_solar_May = month(solardata.datetime) == 5;
% solar = solardata(idx_solar_May,:);

solar=solardata(77839-18*3:104568,12);
solar_1=reshape(table2array(solar),3,8928);
%solar_1=repelem(solar_1,1,5);
solar_11=solar_1(1,:);
solar_11_ramp=[0,diff(solar_11)];

wind=winddata(51893-18*2:69712,12);
wind_1=reshape(table2array(wind),2,8928);
%wind_1=repelem(wind_1,1,5);
wind_11=wind_1(1,:);
wind_11_ramp=[0,diff(wind_11)];

load=loaddata(25957-9:34875,12);
load_1=reshape(table2array(load),1,8928);
%load_11=repelem(load_1,1,5);
load_11_ramp=[0,diff(load_1)];

Ace_star=T_reg_4sec1920{1:end,7}+T_reg_4sec1920{1:end,8}-T_reg_4sec1920{1:end,9};
Ace_star=Ace_star(306721:351360);%May 1min ACE*
a2=reshape(Ace_star,5,8928);Ace_star_5min=sum(a2)'/5; Ace_star_5min=Ace_star_5min';

%Ace_star'
%X=[Ace_star_5min;solar_11;solar_11_ramp;wind_11;wind_11_ramp;load_1;load_11_ramp]; %plot'
X=[solar_11;solar_11_ramp;wind_11;wind_11_ramp;load_1;load_11_ramp];
X_max=max(X,[],2);
X_norm=X./X_max;

% Ace_star_allmin=[];
% for j=1:24
%     for i=1:31
%             Ace_star_allmin(j,:)=[Ace_star_allmin; Ace_star(1+(j-1)*60+(i-1)*24*60:60+j*60+i*24*60)];
%     end
% end
% a2=reshape(Ace_star,15,length(Ace_star)/15);ACEmin=sum(a2)'/15;
% Ace_star=ACEmin;

k1=60;
a1=length(Ace_star);

Ace_star_24=reshape(Ace_star,k1,a1/k1);
Ace_star_24_min=min(Ace_star_24);
Ace_star_24_max=max(Ace_star_24);

for i=1:744
[f,x] = ecdf(Ace_star_24(:,i));
Ace_star_24_max_qua(i) = interp1(f, x, 0.985);
Ace_star_24_min_qua(i) = interp1(f, x, 0.015);
end

X1=[solar_11_ramp;solar_11];

%X=[Ace_star_5min;solar_11;solar_11_ramp;wind_11;wind_11_ramp;load_1;load_11_ramp]; %plot'

output_real=cell2mat(output)+cell2mat(error);
MAPE=(1/length(output_real))*sum(abs(cell2mat(error))./abs(output_real));
nMAPE=(1/length(output_real))*sum(abs(cell2mat(error))./(abs(output_real)*(max(output_real)-min(output_real))));
%%

clearvars Pareto
k1=12;
output_fore=cell2mat(output);
a1=length(output_fore(13:end));

Ace_fore_24=reshape(output_fore(13:end),k1,a1/k1);
Ace_fore_24_min=min(Ace_fore_24);

Ace_fore_24_max=max(Ace_fore_24);

%%%% change location
Ace_fore_24_max=[200 200 Ace_fore_24_max];
ACE_tran=reshape(Ace_fore_24_max,24,31);
ACE_tran16=ACE_tran(1:6,:); a16 = sort(ACE_tran16(:));
ACE_tran2124=ACE_tran(21:24,:); a2124 = sort(ACE_tran2124(:));

nightcap = 491; % 468 491
ACE_tran16=nightcap*ones(6,31);ACE_tran2124=nightcap*ones(4,31);


%Ace_fore_24_max(Ace_fore_24_max<491)=491;
%Ace_fore_24_max(Ace_fore_24_max<0)=0;

%%%
beta = 1;
extra_base = 0;
 
for i= 1:2500
    %nightcap = 350 + i;
    for d=1:31
        for t=1:24
            Ace_fore_24_max_new(1+(d-1)*24:6+(d-1)*24) = nightcap;
            Ace_fore_24_max_new(7+(d-1)*24:20+(d-1)*24) = (1+ 0.001* (i-501)) * Ace_fore_24_max(7+(d-1)*24:20+(d-1)*24) + 75;
        end
    end
%      Ace_fore_24_max_new = (1+ 0.001* (i-501)) * Ace_fore_24_max ;
    
%%%%
%Ace_fore_24_max=[0 0 Ace_fore_24_max];

Estimation=repelem(Ace_fore_24_max_new,1,12);
exceed1=(output_real>Estimation(13:end));
exceed1_total=sum(exceed1)/length(output_real);

procurement_est=sum(Ace_fore_24_max_new);
procurement_real=sum(regup);
Estimation_real=repelem(regup,1,12);
exceed2=(output_real>Estimation_real(13:end));
exceed2_total=sum(exceed2)/length(output_real);
procurement_real=sum(regup(577:end));
procurement_best=sum(Ace_star_24_max(577:end));
procurement_fore=sum(Ace_fore_24_max_new(577:end));

Pro_real=repelem(regup,1,12);
exceed1=(output_real(end-2015:end)>Pro_real(end-2015:end));
exceed1_total=sum(exceed1)/length(output_real(end-2015:end));

Pro_best=repelem(Ace_star_24_max,1,12);
exceed2=(output_real(end-2015:end)>Pro_best(end-2015:end));
exceed2_total=sum(exceed2)/length(output_real(end-2015:end));

Pro_fore=repelem(Ace_fore_24_max_new,1,12);
exceed3=(output_real(end-2015:end)>Pro_fore(end-2015:end));
exceed3_total=sum(exceed3)/length(output_real(end-2015:end));

PRO=[procurement_real,procurement_best,procurement_fore];
EXCEED=[exceed1_total,exceed2_total,exceed3_total];
Pareto(i,:)=[PRO EXCEED];
clearvars PRO EXCEED;
end
Pareto1 = Pareto;
Pareto (:,1:3) = Pareto (:,1:3)/1000;
%figure
scatter (Pareto(1:i,6),Pareto(1:i,3),5 ,'filled'); 
hold on;
scatter (Pareto(1,4),Pareto(1,1),200 ,'filled','d');
hold on;
scatter (Pareto(1,5),Pareto(1,2),300 ,'filled','p');

% scatter (Pareto(501,6),Pareto(501,3),50 ,'filled'); 
% scatter (Pareto(1,6),Pareto(1,3),50 ,'filled'); 
% scatter (Pareto(2500,6),Pareto(2500,3),50 ,'filled'); 

%clearvars Pareto;
% plot(cell2mat(output))
% hold on
% plot(output_real)

% k2=24;
% a2=length(Ace_star_24_min);
% ACE24hour_min=reshape(Ace_star_24_min,k2,a2/k2);
% ACE24hour_max=reshape(Ace_star_24_max,k2,a2/k2);
% for i=1:24
%     subplot(4,6,i)
%     histogram(ACE24hour_max(i,:))
% end

% ACE_STAR4sec=T_reg_4sec{2:end,2};REUP4sec=T_reg_4sec{2:end,4};REUD4sec=T_reg_4sec{2:end,5};ACE4sec=T_reg_4sec{2:end,3};
% k_5_15=225;%75 or 225
% a=length(ACE_STAR4sec);a1=reshape(ACE_STAR4sec,k_5_15,a/k_5_15);ACE_STARmin=sum(a1)'/k_5_15;   
% a2=reshape(REUP4sec,k_5_15,a/k_5_15);REUPmin=sum(a2)'/k_5_15;   
% a3=reshape(REUD4sec,k_5_15,a/k_5_15);REUDmin=sum(a3)'/k_5_15;   
% a4=reshape(ACE4sec,k_5_15,a/k_5_15);ACEmin=sum(a4)'/k_5_15;
% RegPro = readtable('C:\Users\lxh180005\Desktop\solar ramp\OASIS Data\20200201_20200301_AS_REQ_DAM_RegUp+Dn_20210529_03_10_20_v1_Josephine Transformed.xlsx');
% 
% for i=1:720
%     RegPro2(4*(i-1)+1,:)=RegPro1(i,:);
%     RegPro2(4*(i-1)+2,:)=RegPro1(i,:);
%     RegPro2(4*(i-1)+3,:)=RegPro1(i,:);
%     RegPro2(4*(i-1)+4,:)=RegPro1(i,:);
% end
% 
% RegPro3=RegPro2(end-296:end,:);
% exceedbaseline=YTest'>RegPro3(:,2);
% c1=sum(exceedbaseline);c2=sum(RegPro3(:,2));
% a=4;
% YPred1=YPred*a;
% 
% tmp_error_nl = reshape(YPred1(1:296), 4, length(YPred1(1:296))/4);
% error_max = max(tmp_error_nl, [], 1);
% error_min = min(tmp_error_nl, [], 1);error_min=-abs(error_min);
% Error=[error_max' error_min'];
% for i=1:74
%     Error1(4*(i-1)+1,:)=Error(i,:);
%     Error1(4*(i-1)+2,:)=Error(i,:);
%     Error1(4*(i-1)+3,:)=Error(i,:);
%     Error1(4*(i-1)+4,:)=Error(i,:);
% end
% Error1(297,:)=Error1(296,:);
% exceednew=YTest'>Error1(:,1);
% b1=sum(exceednew);b2=sum(Error1(:,1));
% plot(Error1)