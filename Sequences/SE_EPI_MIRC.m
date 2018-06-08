%% This file demonstrates the generation of a .seq file and corresponding toppe files for GE impl.
% Author: Sneha Potdar
% Date: 17th July 2017
%% File Location 
%cd('C:\Users\arush\Desktop\Tools\'); 
tic;
cd('C:\Users\arush\Desktop\MIRC-GE\MATLAB\sequences'); 
% cd('/Users/smysr/Desktop/MIRC-GE/MATLAB/sequences');
addpath(genpath('.'));
addpath(genpath('C:\Users\arush\Desktop\MIRC-GE'));
% sname = 'GRE_pulseq_17082017.seq';
dirname = 'C:\Users\arush\Desktop\MIRC-GE\MATLAB\seqfiles';
% dirname = '/Users/smysr/Desktop/MIRC-GE/MATLAB/seqfiles';
os = 'pc';
%% Set system limits
seq=mr.Sequence();              % Create a new sequence object
fov=220e-3; Nx=64; Ny=64;       % Define FOV and resolution
TE=70e-3;
TR=200e-3;
dt_GE = 4e-6; %seconds
system = mr.opts('MaxGrad',33,'GradUnit','mT/m',...
    'MaxSlew',110,'SlewUnit','T/m/s','ADCdeadtime',10e-6,...
    'RFdeadTime',10e-6,'ADCdeadTime',10e-6);  
% system = mr.opts('MaxGrad',33,'GradUnit','mT/m',...
%     'MaxSlew',100,'SlewUnit','T/m/s','ADCdeadtime',10e-6,...
%     'RFdeadTime',10e-6,'ADCdeadTime',10e-6);  
%% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,system,'Duration',2.5e-3,...
    'SliceThickness',3e-3,'apodization',0.5,'timeBwProduct',4);

%% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;
readoutTime = Nx*dt_GE;  % Nx*dwelltime, here dwelltime is kept as 4e-6
gx = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',readoutTime);
Nx = ceil(gx.flatTime/dt_GE); %to ensure adc dwell time is considered
adc=mr.makeAdc(Nx, 'Dwell', dt_GE, 'Delay',gx.riseTime) ; % This is to ensure dt of 4e-6, number of readout points increase through - harmless but for BW
% adc = mr.makeAdc(Nx,lims,'Duration',gx.flatTime,'Delay',gx.riseTime);

%%  Pre-phasing gradients
preTime=8e-4;
gxPre = mr.makeTrapezoid('x',system,'Area',-gx.area/2-deltak/2,'Duration',preTime);
gzReph = mr.makeTrapezoid('z',system,'Area',-gz.area/2,'Duration',preTime);
gyPre = mr.makeTrapezoid('y',system,'Area',-Ny/2*deltak,'Duration',preTime);

%% Phase blip in shortest possible time
dur = ceil(2*sqrt(deltak/system.maxSlew)/10e-6)*10e-6;
gy = mr.makeTrapezoid('y',system,'Area',deltak,'Duration',dur);

%% Refocusing pulse with spoiling gradients
rf180 = mr.makeBlockPulse(pi,system,'Duration',500e-6);
gzSpoil = mr.makeTrapezoid('z',system,'Area',gz.area*2,'Duration',3*preTime);

%% Calculate delay time
durationToCenter = (Nx/2+0.5)*mr.calcDuration(gx) + Ny/2*mr.calcDuration(gy);
delayTE1=TE/2 - mr.calcDuration(gz)/2 - preTime - mr.calcDuration(gzSpoil) - mr.calcDuration(rf180)/2;
delayTE2=TE/2 - mr.calcDuration(rf180)/2 - mr.calcDuration(gzSpoil) - durationToCenter;

%% Define sequence blocks
seq.addBlock(rf,gz);
seq.addBlock(gxPre,gyPre,gzReph);
seq.addBlock(mr.makeDelay(delayTE1));
seq.addBlock(gzSpoil);
seq.addBlock(rf180);
seq.addBlock(gzSpoil);
seq.addBlock(mr.makeDelay(delayTE2));
for i=1:Ny
    seq.addBlock(gx,adc);           % Read one line of k-space
    seq.addBlock(gy);               % Phase blip
    gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
end
seq.addBlock(mr.makeDelay(1));

%%
% seq.plot('TimeRange',[0 TR]);
% fname =[date,'SE-EPI_Mod_Pulseq.seq'];
% seq.write(fname);
toc;
%%  Call toppe functions to generate reqd .wav, cores.txt, scanloop.txt
% generate files for GE scanner
% mr_system = system;
% clear system; %this will be misinterpreted later else
% seq2ge_mod1(fname);

%%
%   seq2ge_mod1('SE_EPI_Python_01092017.seq');
    seq2ge_mod1('SE_EPI_GPI_01092017.seq');
%%  simulate GE scan 
% system('tar xzf GEscan.tgz');
fprintf(1, 'Simulating... ');
ny = 1;
nBlocksPerTR = 136;%from the for loop above on the pulseq
d = readloop('scanloop.txt');
fprintf(1,'\n');
for ii = 1:1:ny
	fprintf(1, '\r%d of %d', ii, ny);
	scansim(1+(ii-1)*nBlocksPerTR, nBlocksPerTR+(ii-1)*nBlocksPerTR, d);
	pause(0.1);  % to allow display to update
end
fprintf(1,'\n');

%% Package the tgz
% system('tar czf GEscan.tgz cores.txt scanloop.txt *.wav seq2ge.m se-epi_demo.m');
switch os
    case 'pc'
        zip('SE_EPI_GPI_08092017', {'modules.txt','scanloop.txt','*.mod'});
        delete('*.mod','modules.txt','scanloop.txt');
    case 'mac'
        system('tar czf GEscan_SE-EPI.tgz cores.txt scanloop.txt *.mod' );
        system('rm *.mod cores.txt scanloop.txt');
end
%%
% seq.plot();
% seq.write('SE-EPI_MIRC.seq');
% seq.write('/Users/sravan953/Documents/MIRC/Projects/pulseq-gpi/se_epi_matlab.seq');   % Output sequence for scanner
% seq.plot();             % Plot sequence waveforms
