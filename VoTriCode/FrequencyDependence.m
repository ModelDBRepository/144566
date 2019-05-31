% Code written by Claudia Clopath
% Please cite:	
% Clopath C, Vasilaki E, Büsing L and Gerstner W.
% "Connectivity reflects coding: A model of voltage-based
% spike-timing-dependent-plasticity with homeostasis"
% Nature Neuroscience, 13, 344-352, 2010
% and
% Clopath C, Gerstner W.
% "Voltage and Spike Timing interact in STDP - a Unified Model."
% Frontiers in Synaptic Neuroscience, doi:10.3389/fnsyn.2010.00025, 2010
 
% The plasticity model developed in Clopath et al is used to reproduce the 
% frequency dependency shown by Sjoestroem et al Neuron 2001. 

% The code uses the two files: VoTri.m (voltage-triplet plasticity rule developed in Clopath et al. 2010) and aEIF.m (neuron model)

clear all
par=[0.00014   0.00008  7 15  10];  % parameters of the plasticity model
rho=[0.1, 10, 20, 40, 50];          % frequency of the pairing for the simulation
DeltaTau=10;                        % time lag


%%%%%%%%%%%%%%%%%%%%% pre-post pairing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation of the model
nFreq=size(rho,2);
w=zeros(1,nFreq);
for i=1:nFreq
    w(i)=(VoTri(rho(i), DeltaTau, par)-0.5)/0.5;    
end


% Plot
figure; hold on
h=plot(rho,w*100+100);
set(h,'linewidth',3);
set(gca,'fontsize',20);
%%%%%%%%%%%%%%%%%%%%% post-pre pairing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation of the model
for i=1:nFreq
    w(i)=(VoTri(rho(i), -DeltaTau,par)-0.5)/0.5;
end

% Plot
h=plot(rho,w*100+100,'r');
set(h,'linewidth',2);
set(gca,'fontsize',20);
axis tight
xlabel('\rho [Hz]','fontsize',20);
ylabel('normalized weight [%]','fontsize',20);