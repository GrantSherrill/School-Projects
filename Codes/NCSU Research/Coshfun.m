%% This top part has been commendeered for the error calculation on the 042722 experiment
%Numbers being divided in denominator may seem random, but they are the
%analytical gradient contained in the mu verification 05/02/22 file 
er1=[0.01, 0.014 ];
er2=[0.0122,	0.0122];
er3=[0.0209	,0.0131,	0.0153	,0.0142];
er4=[0.0167,	0.0282,	0.0195];
er5=[0.0203	,0.0235,	0.0299];
er6=[0.0408,	0.0234];
er7=[0.034	,0.033];


er1=2*sqrt(sum(er1.^2))/(length(er1)*0.375);
er2=2*sqrt(sum(er2.^2))/(length(er2)*1.505);
er3=2*sqrt(sum(er3.^2))/(length(er3)*2.316);
er4=2*sqrt(sum(er4.^2))/(length(er4)*3.192);
er5=2*sqrt(sum(er5.^2))/(length(er5)*3.799);
er6=2*sqrt(sum(er6.^2))/(length(er6)*5.407);
er7=2*sqrt(sum(er7.^2))/(length(er7)*7.466);

er=[er1,er2,er3,er4,er5,er6,er7];
%% 




mo=[0.17,0.047,0.042,0.0645,0.107,0.086,.105]';
B=[8,25,35,45,52,70,95]';

% B2=[10.6,7.6,4.86,2.16,7.73,6.38,5,1.34,2.7,4,4.7]';
% mo2=[ 56.26769627, 28.88888889 ,17.36943907,12.56568779,34.34108527,13.29896907,10.81395349,9.580573951,20.33505155,25.85139319,8.865979381]';
% 
% 



b=linspace(-1000,1000,1000);
b0=500;
sat=1;
Mag=sat.*(coth(b./b0)-1./(b./b0));


figure()
hold on
plot(B,mo,'.','MarkerSize',30)
% plot(B2,mo2,'.','MarkerSize',24)
%plot(b,Mag)

legend('\mu Values')
    title('\mu Values as a Function of Field ')
    xlabel('Field(mT)')
    ylabel('Moment pJ/T')
    h = gca;
    set(h,'FontSize',[40]);
hold off

%% 

linguess=[1,1];
modelfun = @(lcon) linfun(lcon , mo,B);
[lcons] = fminsearch(modelfun,linguess);
linapprox = lcons(1)+( b.* lcons(2));

modelfun2 = @(lcon) linfun(lcon , mo2,B2);
[lcons2] = fminsearch(modelfun2,linguess);
linapprox2 = lcons2(1)+( b.* lcons2(2));

figure()
hold on
plot(B,mo,'.','MarkerSize',30)
%plot(B2,mo2,'.','MarkerSize',30)
plot(b,linapprox,'-','LineWidth',1.8,'Color','k')
%plot(B2,linapprox2,'-','MarkerSize',1.8)

legend('\mu Values', 'Linear Fit' )
    title('\mu Values as a Function of Field ')
    xlabel('Field(mT)')
    ylabel('Moment fJ/T')
    h = gca;
    set(h,'FontSize',[40]);
hold off

%% 
coshguess=[73];


modelfun = @(ccon) Coshfun(ccon , mo,B);
[ccons] = fminsearch(modelfun,coshguess);
coshapprox =sat.*(coth(B./(ccons(1)))-1./(B./ccons(1)));

% modelfun2 = @(ccon) Coshfun(ccon , mo2,B2);
% [ccons2] = fminsearch(modelfun2,coshguess);
% coshapprox2 = sat.*(coth(B2./(ccons2(1)))-1./(B2./ccons2(1)));

coshex=sat.*(coth(b./(ccons(1)))-1./(b./ccons(1)));
% coshex2=sat.*(coth(b./(ccons2(1)))-1./(b./ccons2(1)));
hold on
plot(B,mo,'.','MarkerSize',24)
plot(B,coshapprox,'-','LineWidth',1.8,'Color','k')
plot(b,coshex,'-','LineWidth',1.8,'Color','k')
errorbar(B,mo,er,'Linestyle','none')
legend('\mu Values', 'Coth Fit')
    title('\mu Values as a Function of Field ' )
    xlabel('Field(mT)')
    ylabel('Moment pJ/T')
    xlim([0 300])
    ylim([0 1])
    h = gca;
    set(h,'FontSize',[40]);



hold off
%% 
coshguess=[500,1000];


modelfun = @(ccon) Coshfun2(ccon , mo,B);
[ccons] = fminsearch(modelfun,coshguess);
coshapprox =ccons(2).*(coth(B./(ccons(1)))-1./(B./ccons(1)));



coshex=ccons(2).*(coth(b./(ccons(1)))-1./(b./ccons(1)));

hold on
plot(B,mo,'.','MarkerSize',24)

plot(B,coshapprox,'-','LineWidth',1.8,'Color','k')
plot(b,coshex,'-','LineWidth',1.8,'Color','b')

hold off
%% 


quadguess=[1,1,1];
modelfun = @(qcon) quadfun(qcon , mo,B);
[qcons] = fminsearch(modelfun,quadguess);
quadapprox = qcons(1)+( B.* qcons(2))+((B.^2).*qcons(3));

modelfun2 = @(qcon) quadfun(qcon , mo2,B2);
[qcons2] = fminsearch(modelfun2,quadguess);
quadapprox2 = qcons2(1)+( B2.* qcons2(2))+((B2.^2) .*qcons2(3));
figure()
hold on
plot(B,mo,'.','MarkerSize',24)
plot(B2,mo2,'.','MarkerSize',24)
plot(B,quadapprox,'-','LineWidth',1.8,'Color','k')
plot(B2,quadapprox2,'-','MarkerSize',1.8)
hold off


%% 

function J = Coshfun(par,obs,b)

sat=1;
y=sat.*(coth(b./par(1))-1./(b./par(1)));
J = (y-obs)'*((y-obs));

end 

function J = Coshfun2(par,obs,b)

y=par(2).*(coth(b./par(1))-1./(b./par(1)));
J = (y-obs)'*((y-obs));

end 

function J = linfun(um,obs,rho)

y=um(1)+(rho.*um(2));
J = (y-obs)'*((y-obs));

end 


function J = quadfun(um,obs,rho)

y=um(1)+(rho.*um(2))+(rho.^2 .*um(3)) ;

J = (y-obs)'*((y-obs));

end