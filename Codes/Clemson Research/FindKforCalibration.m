
clear;clc;close all;
%Input these with each file run
addpath('matlab/210708')
file1 = '12_28_47_frames.tifResults.csv';
file2 = '210708 70glycerol10 12_28_47 stage pos Z Step   0.csv';
file3= '12_28_47_times.txt';    %keep in mind this is how to get from frames to ms
frequency=2;%in hz
omega=2*pi*frequency;
fr=120;
time=6;  %in sec
waittime=1; %in sec
amplitude= .25; %in um
pixeltoum=.0985; % pixel conversion to um 
beaddata = importdata(strcat(file1)); %Stage position according to the tweezer in x
beadx=beaddata.data(:,4);

beadx=beadx*pixeltoum; %convert pixels to microns

Waitfrac= waittime/time;


Stagedata= importdata(file2);
Stagex=Stagedata.data(:,2);
Stagetime=Stagedata.data(:,1);
Stagemove=Stagetime(round(Waitfrac*length(Stagetime)):length(Stagetime))*1/1000; %1/1000 can be applied here or multiply by 1000 in the frequency. The data is put in ms, and frequency is s^-1

Timedata=importdata(strcat(file3));
movetime=Timedata((Waitfrac*length(Timedata)):length(Timedata)-(length(Timedata)-fr*time+Timedata(1,2)),1)*1/fr; %this gives a column of time in ms when the stage is moving 

%l=length(Stagex);
%% Split data sets 
Stagexp=Stagex(1:Waitfrac*length(Stagex)); %p referes to passive, a to active
Stagexa=Stagex(Waitfrac*length(Stagex):length(Stagex));

Stageavg=mean(Stagexa);
Stagexn=Stageavg-Stagexa;

beadxp=beadx(1:Waitfrac*length(beadx));
beadxa=beadx(Waitfrac*length(beadx):length(Timedata)-(length(Timedata)-fr*time+Timedata(1,2)));
%% find trap position
trapx=mean(beadxp);
delx=-(beadxa-trapx); %VERY IMPORTANT: this is our delta x, our effective movement away from the trap position, which is what our spring force is given by
maxx=max(delx);
%% Note this is the first program for calibrating the data, so instead of finding the force using f=kx, k needs to be defined


myFittype = fittype('a*sin(b*x-p)+c','coefficients',{'a','p','c'},'problem',{'b'} );
testfitstage=fit(Stagemove, Stagexn, myFittype,'problem',{omega},'StartPoint', [amplitude,0,trapx]);
%plot(testfitstage,Stagemove,Stagexn)
stagecoeff= coeffvalues(testfitstage);
figure()
testfitbead=fit(movetime, delx, myFittype, 'problem', {omega},'StartPoint', [maxx,0,trapx]);
plot(testfitbead,movetime,delx)
beadcoeff=coeffvalues(testfitbead);

phasecheck=abs((stagecoeff(1,2)-beadcoeff(1,2))-pi/2); %Double check this is close to 0, this is needed for the approximation that water is viscous and should be 90° out of phase

% figure()
% hold on
% yyaxis left
% plot(Stagemove,Stagexn)
% ylabel('Displacement (um)')
% 
% yyaxis right
% plot(movetime, delx,'.')
% ylabel('Displacement (um)')
% xlabel('Time (s)')
% 
% if phasecheck>pi
%     phasecheck=abs(phasecheck-pi);
% end
% if phasecheck>pi/2
%     phasecheck=abs(pi-phasecheck);
%     
% end
%% Now that we have our amplitudes, and we have double checked that our two data sets are approximately out of phase we can find our k value. 
%we use the equation k=6*pi*r*mu*v/x where v and x are the fit functions.
%They are completely out of phase and have the same wavelngth so the
%sinusoids are canceled (due to d/dx sin(x)=sin(x+pi/2)
%note here both amplitudes are in terms of um so the units cancel
r = 3.4e-06/2; %in m
mu= .0225; %in N*s/m^2

k=abs(6*pi*omega*r*mu*stagecoeff(1,1)/beadcoeff(1,1))*10^6; %convert to pN/um



