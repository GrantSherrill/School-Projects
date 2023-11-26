clear;clc;close all;

addpath('matlab/210630')
addpath('matlab/210630 retracked')
% time1=importdata('19_28_20_times.txt');
bead1 = importdata('Resultsrack.csv');
% time2=importdata('19_27_41_times.txt');
% bead2 = importdata('Bretrack1.csv');
% time3=importdata('19_27_56_times.txt');
% bead3 = importdata('Bretrack2.csv');
% time8=importdata('19_28_09_times.txt');
% bead8 = importdata('Bretrack3.csv');
% time9=importdata('19_28_20_times.txt');
% bead9 = importdata('Bretrack4.csv');
% time4=importdata('19_28_35_times.txt');
% bead4 = importdata('Bretrack 5.csv');
% time5=importdata('19_28_45_times.txt');
% bead5 = importdata('Bretrack 6.csv');
% % time10=importdata('20_34_25_times.txt');
% bead10 = importdata('20_34_25_frames.tifResults.csv');
% 
% time6=importdata('20_35_20_times.txt');
% bead6 = importdata('20_35_20_frames.tifResults.csv');
% time7=importdata('20_35_33_times.txt');
% bead7 = importdata('20_35_33_frames.tifResults.csv');
% 





time1=bead1.data(:,1);
% time2=time2(:,1);
% time3=time3(:,1);
% time4=time4(:,1);
% time5=time5(:,1);
% % time6=time6(:,1);
% % time7=time7(:,1);
% time8=time8(:,1);
% time9=time9(:,1);
% time10=time10(:,1);



bead1=bead1.data(:,4);
plot(time1,bead1)
% bead2=bead2.data(:,4);
% bead3=bead3.data(:,4);
% bead4=bead4.data(:,4);
% bead5=bead5.data(:,4);
% % bead6=bead6.data(:,4);
% % bead7=bead7.data(:,4);
% bead8=bead8.data(:,4);
% bead9=bead9.data(:,4);
% % bead10=bead10.data(:,4);
% 
% 
% 
% subplot(4,2,1)
% plot(time1,bead1)
% subplot(4,2,2)
% plot(time1,bead2)
% subplot(4,2,3)
% plot(time1,bead3)
% subplot(4,2,4)
% plot(time1,bead8)
% subplot(4,2,5)
% plot(time1,bead9)
% subplot(4,2,6)
% plot(time1,bead4)
% subplot(4,2,7)
% plot(time1,bead5)
% % subplot(4,2,8)
% % plot(time10,bead10)
% % subplot(5,2,9)
% % plot(time6,bead6)
% % subplot(5,2,10)
% % plot(time7,bead7)
% 
% 











