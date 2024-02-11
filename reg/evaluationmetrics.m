%function f=evaluationmetrics()
%last 7 days
procurement_real=sum(regup(528:end));
procurement_best=sum(Ace_star_24_max(528:end));
procurement_fore=sum(Ace_fore_24_max(528:end));

Pro_real=repelem(regup,1,12);
exceed1=(output_real(end-2015:end)>Pro_real(end-2015:end));
exceed1_total=sum(exceed1)/length(output_real(end-2015:end));

Pro_best=repelem(Ace_star_24_max,1,12);
exceed2=(output_real(end-2015:end)>Pro_best(end-2015:end));
exceed2_total=sum(exceed2)/length(output_real(end-2015:end));

Pro_fore=repelem(Ace_fore_24_max,1,12);
exceed3=(output_real(end-2015:end)>Pro_fore(end-2015:end));
exceed3_total=sum(exceed3)/length(output_real(end-2015:end));

PRO=[procurement_real,procurement_best,procurement_fore];
EXCEED=[exceed1_total,exceed2_total,exceed3_total];
Pareto=[PRO;EXCEED];

%scatter(PRO,EXCEED)

% k1=12;
% output_fore=cell2mat(output);
% a1=length(output_fore(13:end));
% 
% Ace_fore_24=reshape(output_fore(13:end),k1,a1/k1);
% Ace_fore_24_min=min(Ace_fore_24);
% 
% Ace_fore_24_max=max(Ace_fore_24);
% Ace_fore_24_max=[0 0 Ace_fore_24_max];
% 
% Estimation=repelem(Ace_fore_24_max,1,12);
% exceed1=(output_real>Estimation(13:end));
% exceed1_total=sum(exceed1)/length(output_real);
% 
% procurement_est=sum(Ace_fore_24_max);
% procurement_real=sum(regup);
% Estimation_real=repelem(regup,1,12);
% exceed2=(output_real>Estimation_real(13:end));
% exceed2_total=sum(exceed2)/length(output_real);

%end