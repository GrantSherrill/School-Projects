clear;clc;close all;
%Input these with each file run
% Reminder to change graph title with each run
addpath('matlab/210708')
file1 = '11_46_51_frames.tifResults.csv';
file2 = '210708 70glycerol1 11_46_51 stage pos Z Step   0.csv';
file3= '11_46_51_times.txt';  %keep in mind this is how to get from frames to ms
name=['210708 70glycerol11'  '.svg']; 
name2=['210708 70glycerol11' num2str(2) '.svg' ];
%saveName='210608 11_38_03.mat';
frequency=2;%in hz
omega=2*pi*frequency;
fr=120;
time=6;  %in sec, note thisis the total time, the time plus the wait time
waittime=1; %in sec
amplitude= .1; %in um
r=3.4/2; %radius of the bead in um
k=196.56; %Find this by running the calibration program first, this is in pN/um
pixeltoum=.0985; % pixel conversion to um 
beaddata = importdata(strcat(file1)); %Stage position according to the tweezer in x
beadx=beaddata.data(:,4);
beadx=beadx*pixeltoum; %convert pixels to microns

Waitfrac= waittime/time;


Stagedata= importdata(file2);
Stagex=Stagedata.data(:,2);
Stagetime=Stagedata.data(:,1);
Stagemove=Stagetime(round(Waitfrac*length(Stagetime)):length(Stagetime))*1/1000; %1/1000 can be applied here or multiply by 1000 in the frequency. The data is output in ms, and frequency is s^-1

Timedata=importdata(strcat(file3));
movetime=Timedata(round(Waitfrac*(length(Timedata)-Timedata(1,2))):length(Timedata)-(length(Timedata)-fr*time+Timedata(1,2)),1)*1/fr; %this gives a column of time in ms when the stage is moving 

%l=length(Stagex);
%% Split data sets 
Stagexp=Stagex(1:round(Waitfrac*length(Stagex))); %p referes to passive, a to active
Stagexa=Stagex(round(Waitfrac*length(Stagex)):length(Stagex));

Stageavg=mean(Stagexa);
Stagexn=Stagexa-Stageavg;%stage x normalized

beadxp=beadx(1:round(Waitfrac*length(beadx)));
beadxa=beadx(round(Waitfrac*(length(Timedata)-Timedata(1,2))):fr*time-Timedata(1,2));
%% find trap position
trapx=mean(beadxp);%in um
delx=-(beadxa-trapx); %VERY IMPORTANT: this is our delta x, our effective movement away from the trap position, which is what our spring force is given by
%% This is where the code differs from the calibration data. As this point, we have k, and so are wondering about the phase shift

f=k*delx; %Multiply by 1*10^-6 to get from um to m (delx is in um, k is in N/m), and then by 1*10^12 to get from N to pN
fdis=max(f);
xdis=max(Stagexn);
fig = figure('WindowState','maximized');
pause(1);
fig.Position
subplot(2,2,1);
hold on
yyaxis left
title('Stage data and tracked data')
xlabel('Time(s)')
ylabel('Stage Displacement (um)')
plot(Stagemove,Stagexn)

yyaxis right
plot(movetime, delx)
ylabel('Bead Movement (um)')
legend('Stage Position', 'Bead Position');
hold off
myFittype = fittype('a*sin(b*x-p)+c','dependent',{'y'},'independent',{'x'},... %Fit both A and p and c
    'coefficients',{'a','p','c'},'problem',{'b'});
[fitstage,stagestat]=fit(Stagemove, Stagexa, myFittype,'problem',{omega},'StartPoint', [xdis,0,0]);
%plot(fitstage,Stagemove,Stagexa)
stagecoeff= coeffvalues(fitstage);
stageci=confint(fitstage);
% 
% fitbead=fit(movetime, delx, myFittype, 'problem', {omega},'StartPoint', [1,0,trapx]);
% %plot(fitbead,movetime,delx)
% beadcoeff=coeffvalues(fitbead);
% 

[fitbead2,beadstat]=fit(movetime, f, myFittype, 'problem', {omega},'StartPoint', [fdis,0,trapx]); 
%plot(fitbead2,movetime,f)
beadcoeff2=coeffvalues(fitbead2);
beadci=confint(fitbead2);

phasediff=abs(stagecoeff(1,2)-beadcoeff2(1,2));%If this returns a value larger than pi,2 this can be investigated

if phasediff>pi
    phasediff=phasediff-pi;
end
if phasediff>pi/2
    phasediff=abs(pi-phasediff);
    
end
%% Now thaat we have the phase shift, we want to construct our producables- G', G'', and tan(delta)
%remember that G'=fmax*cos(delphi)/6*pi*r*xmax

fmax=abs(beadcoeff2(1,1));
xmax=abs(stagecoeff(1,1)); 

beadup=beadci(1,1)*sin(omega*movetime-beadci(1,2))+beadci(1,3);
beaddown=beadci(2,1)*sin(omega*movetime-beadci(2,2))+beadci(2,3);

gprime=fmax*cos(phasediff)/(6*pi*r*xmax); %Note here that we are in pN/um^2, which to get to N/m^2 we would multiply by 1*10^-12*1*10^6*1*10^6 meaning this value is correct for SI.
gdprime=fmax*sin(phasediff)/(6*pi*r*xmax);
tandel=sin(phasediff)/cos(phasediff);

if frequency*(time-waittime)<=15
    subplot(2,2,2)

    plot(fitstage,Stagemove,Stagexa)
    title('Stage Fit')
    xlabel('Time (s)')
    ylabel('Displacement (um)')

    subplot(2,2,4)


    plot(fitbead2,movetime,f)
    ylabel('Force (pN)')
    title('Force Fit')
    xlabel('Time(s)')

    subplot(2,2,3)
    hold on

    yyaxis left
    xlim([0 time]);
    plot(fitstage,'b')
    title('Stage data and tracked data')
    xlabel('Time(s)')
    ylabel('Stage Displacement (um)')



    yyaxis right
    xlim([0 time]);
    plot(fitbead2)
    ylabel('Force (pN)')
   
    text(1,fmax/2, sprintf('Tandel=%f', tandel))
    
    plot(movetime, beadup, '--') 
    plot(movetime, beaddown, '--')
    legend('Stage Displacement', 'Force','','');
    hold off
elseif frequency*(time-waittime)>15
  subplot(2,2,2)

    plot(fitstage,Stagemove,Stagexn)
    xlim([waittime 10/frequency+waittime]);
    title('Stage Fit')
    xlabel('Time (s)')
    ylabel('Displacement (um)')

    subplot(2,2,4)


    plot(fitbead2,movetime,f)
    xlim([waittime 10/frequency+waittime]);
    ylabel('Force (pN)')
    title('Force Fit')
    xlabel('Time(s)')

    subplot(2,2,3)
    hold on

    yyaxis left
    xlim([0 10/frequency]);
    plot(fitstage,'b')
    title('Stage data and tracked data')
    xlabel('Time(s)')
    ylabel('Stage Displacement (um)')



    yyaxis right
    xlim([0 10/frequency]);
    plot(fitbead2)
    ylabel('Force (pN)')
    
    text(1,fmax/2, sprintf('Tandel=%f', tandel))
    plot(movetime, beadup, '--') 
    plot(movetime, beaddown, '--')
    legend('Stage Displacement', 'Force','','');
    hold off  
end
% saveas(gcf, name);



if frequency*(time-waittime)<=15
    hold on
    fig = figure('WindowState','maximized');
pause(1);
fig.Position

    yyaxis left
    xlim([0 time]);
    plot(fitstage,'b')
    title('Stage data and tracked data')
    xlabel('Time(s)')
    ylabel('Stage Displacement (um)')



    yyaxis right
    xlim([0 time]);
    plot(fitbead2)
    ylabel('Force (pN)')
    legend('Stage Position', 'Force');
    
    dim=[.15 .2 .2 .2];
    str={sprintf('Tandel=%f', tandel), sprintf('Gprime=%f', gprime),sprintf('Gdprime=%f', gdprime),sprintf('rsqr=%f', beadstat.rsquare)};
    annotation('textbox',dim,'BackgroundColor','w','String', str, 'FitBoxToText','on');
    hold off
elseif frequency*(time-waittime)>15
    hold on
    fig = figure('WindowState','maximized');
pause(1);
fig.Position

    yyaxis left
    xlim([0 10/frequency]);
    plot(fitstage,'b')
    title('Stage data and tracked data')
    xlabel('Time(s)')
    ylabel('Stage Displacement (um)')



    yyaxis right
    xlim([0 10/frequency]);
    plot(fitbead2)
    ylabel('Force (pN)')
    
    legend('Stage Position', 'Force');
   
    dim=[.15 .2 .2 .2];
    str={sprintf('Tandel=%f', tandel), sprintf('Gprime=%f', gprime),sprintf('Gdprime=%f', gdprime),sprintf('rsqr=%f', beadstat.rsquare)};
    annotation('textbox',dim,'BackgroundColor','w','String', str, 'FitBoxToText','on');
    hold off
end 
% saveas(gcf, name2);













