
clear;clc;close all;
%Input these with each file run
file1 = 'Traj_11_38_03_frames.tifResults.csv';
file2 = '210608 0_ glycerol    11_38_03 stage pos Z Step   0.csv';
file3= '11_38_03_times.txt';  %keep in mind this is how to get from frames to ms
saveName='210608 11_38_03.mat';
frequency=1;%in hz
wavelength=2*pi*frequency;
time=10;  %in sec
waittime=5; %in sec
amplitude= 5; %in um
pixeltoum=.0985; % pixel conversion to um 
beaddata = importdata(strcat(file1)); %Stage position according to the tweezer in x
beadx=beaddata.data(:,4);
beadx=beadx*pixeltoum; %convert pixels to microns

Waitfrac= waittime/time;


Stagedata= importdata(file2);
Stagex=Stagedata.data(:,2);
Stagetime=Stagedata.data(:,1);
Stagemove=Stagetime(Waitfrac*length(Stagetime):length(Stagetime))*1/1000; %1/1000 can be applied here or multiply by 1000 in the frequency. The data is put in ms, and frequency is s^-1

Timedata=importdata(strcat(file3));
movetime=Timedata(Waitfrac*length(Timedata):length(Timedata),2)*1/1000; %this gives a column of time in ms when the stage is moving 

l=length(Stagex);
%% Split data sets 
Stagexp=Stagex(1:Waitfrac*length(Stagex)); %p referes to passive, a to active
Stagexa=Stagex(Waitfrac*length(Stagex):length(Stagex));

Stageavg=mean(Stagexa);
Stagexn=Stageavg-Stagexa;

beadxp=beadx(1:Waitfrac*length(beadx));
beadxa=beadx(Waitfrac*length(beadx):length(beadx));
%% find trap position
trapx=mean(beadxp);
delx=beadxa-trapx; %VERY IMPORTANT: this is our delta x, our effective movement away from the trap position, which is what our spring force is given by
%% Note this is the first program for calibrating the data, so instead of finding the force using f=kx, k needs to be defined


myFittype = fittype('a*sin(b*x-p)','dependent',{'y'},'independent',{'x'},... %Fit both A and p
    'coefficients',{'a','p'},'problem',{'b'});
testfitstage=fit(Stagemove, Stagexn, myFittype,'problem',{wavelength});
%plot(testfitstage,Stagemove,Stagexn)
stagecoeff= coeffvalues(testfitstage);

testfitbead=fit(movetime, delx, myFittype, 'problem', {wavelength});
%plot(testfitbead,movetime,delx)
beadcoeff=coeffvalues(testfitbead);

phasecheck=abs(stagecoeff(1,2)-beadcoeff(1,2))-pi/2; %Double check this is close to 0, this is needed for the approximation that water is viscous and should be 90° out of phase

%% Now that we have our amplitudes, and we have double checked that our two data sets are approximately out of phase we can find our k value. 
%we use the equation k=6*pi*r*mu*v/x where v and x are the fit functions.
%They are completely out of phase and have the same wavelngth so the
%sinusoids are canceled (due to d/dx sin(x)=sin(x+pi/2)
%note here both amplitudes are in terms of um so the units cancel
r = 3.4e-06; %in m
mu= 1.0016e-03; %in N*s/m^2

k=6*pi*wavelength*r*mu*stagecoeff(1,2)/beadcoeff(1,2);



