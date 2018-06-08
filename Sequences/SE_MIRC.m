%% This file demonstrates the generation of a .seq file and corresponding toppe files for GE impl.
% Author: Sneha Potdar
% Date: 17th July 2017
%% File Location 
tic;
pulseq_dir = uigetdir();
pulseq_dir = [pulseq_dir,'/.'];
% cd('C:\Users\arush\Desktop\MIRC-GE\MATLAB\sequences'); 
% cd('/Users/smysr/Desktop/MIRC-GE/MATLAB/sequences');
addpath(genpath('.'));
addpath(genpath(pulseq_dir));
% addpath(genpath('C:\Users\arush\Desktop\MIRC-GE'));
% % sname = 'GRE_pulseq_17082017.seq';
% dirname = 'C:\Users\arush\Desktop\MIRC-GE\MATLAB\seqfiles';
% % dirname = '/Users/smysr/Desktop/MIRC-GE/MATLAB/seqfiles';
os = 'pc';
%% System limits
seq=mr.Sequence();
fov=256e-3;
Nx=256; Ny=256;
TE=100e-3;
TR=2000e-3;
dt_GE = 4e-6; %seconds
system = mr.opts('MaxGrad',33,'GradUnit','mT/m',...
        'MaxSlew',100,'SlewUnit','T/m/s','ADCdeadtime',10e-6,...
        'RFdeadTime',10e-6);

%% Implementation of RF 90 degree pulse and Gz slice selelction
deltak=1/fov;
kWidth = Nx*deltak;
readoutTime = Nx*dt_GE;  % Nx*dwelltime, here dwelltime is kept as 4e-6
[rf, gz] = mr.makeSincPulse(pi/2,system,'Duration',2.5e-3,...
    'SliceThickness',3e-3,'apodization',0.5,'timeBwProduct',4);

%% Implementation of Gx Frequency encoding, ADC
gx = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',readoutTime);
% Nx = ceil(gx.flatTime/dt_GE); %to ensure adc dwell time is considered
% adc=mr.makeAdc(Nx, 'Dwell', dt_GE, 'Delay',gx.riseTime) ; % This is to ensure dt of 4e-6, number of readout points increase through - harmless but for BW
adc = mr.makeAdc(Nx,system,'Duration',gx.flatTime,'Delay',gx.riseTime);

%% Implementation of GxPre, GzReph
% gxPre = mr.makeTrapezoid('x',system,'Area',gx.area/2,'Duration',readoutTime./2);
gxPre = mr.makeTrapezoid('x',system,'Area',-gx.area/2,'Duration',2.5e-3);
gzReph = mr.makeTrapezoid('z',system,'Area',-gz.area/2,'Duration',2.5e-3);

%% Implementation of RF180 degree, Gz180 slice selelction
[rf180, gz180] = mr.makeSincPulse(pi,system,'Duration',2.5e-3,...
  'SliceThickness',3e-3,'apodization',0.5,'timeBwProduct',4);

%% Calculating Delays
delayTE1=TE/2-mr.calcDuration(gzReph)-mr.calcDuration(rf)- mr.calcDuration(rf180)/2;
delayTE2=TE/2-mr.calcDuration(gx)./2-mr.calcDuration(rf180)/2;
delayTE3=TR-TE-mr.calcDuration(gx);

%%
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;
for i=1:Ny
    seq.addBlock(rf,gz);
%     gyPre = mr.makeTrapezoid('y',system,'Area',-(Ny/2-(i-1))*deltak,...
%         'Duration',readoutTime./2);
    gyPre = mr.makeTrapezoid('y',system,'Area',phaseAreas(i),'Duration',2.5e-3);
    seq.addBlock(gxPre,gyPre,gzReph);
    seq.addBlock(mr.makeDelay(delayTE1));
    seq.addBlock(rf180,gz180);
    seq.addBlock(mr.makeDelay(delayTE2));
    seq.addBlock(gx,adc);
    seq.addBlock(mr.makeDelay(delayTE3));
end
%%
% seq.plot('TimeRange',[0 2.*TR]);
% fname =[date,'SE_Mod_Pulseq.seq'];
% seq.write(fname);
toc;
%%  Call toppe functions to generate reqd .mod, modules.txt, scanloop.txt
% generate files for GE scanner
% mr_system = system;
% clear system; %this will be misinterpreted later else
% seq2ge_mod1(fname);
%%
%  seq2ge_mod1('SE_Python_22082017.seq');
  seq2ge_mod1('se_gpi_23082017.seq');

%%  simulate GE scan 
% system('tar xzf GEscan.tgz');
fprintf(1, 'Simulating... ');
ny = 256;
nBlocksPerTR = 7;%from the for loop above on the pulseq
d = readloop('scanloop.txt');
fprintf(1,'\n');
for ii = 1:1:ny
	fprintf(1, '\r%d of %d', ii, ny);
	scansim(1+(ii-1)*nBlocksPerTR, nBlocksPerTR+(ii-1)*nBlocksPerTR, d);
	pause(0.1);  % to allow display to update
end
fprintf(1,'\n');

%% Package the tgz
% system('tar czf GEscan.tgz modules.txt scanloop.txt *.wav seq2ge.m gre_demo.m');
switch os
    case 'pc'
        zip('SE_GPI_08092017', {'modules.txt','scanloop.txt','*.mod'});
        delete('*.mod','modules.txt','scanloop.txt');
%         system('zip GEscan_SE.tgz modules.txt scanloop.txt *.mod' );
%         system('del *.mod modules.txt scanloop.txt');
    case 'mac'
        system('tar czf GEscan_SE.tgz modules.txt scanloop.txt *.mod' );
        system('rm *.mod modules.txt scanloop.txt');
end
%%
% seq.plot();
% seq.write('Spin_Echo_MIRC.seq');
% seq.plot('TimeRange',[0 TR]);
% seq.write('//Users/sravan953/Documents/MIRC/Projects/pulseq-gpi/se_matlab.seq');
