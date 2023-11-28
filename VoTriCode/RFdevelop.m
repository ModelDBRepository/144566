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

% This code uses the voltage-triplet plasticity rule developed in Clopath
% et al. 2010 in a feedforward network to develop receptive fields (note
% that the amplitude for potentation and depression are multiplied by 10 to speed up the simulations).

% This code calls the function aEIF.m which is the neuron model

% clear all
% Simulation parameters
Ni        = 500;                                % number of inputs
n_neur    = 1;                                  % number of output neurons
EP        = 100;                                % number of time steps per epoch [ms]
NEP       = 1000;                               % number of epochs
dt        = 1;                                  % duration of time step [ms]
T         = EP*NEP;                             % number of time steps
tau_th    = 1.2*1000.0;                         % time constant for sliding threshold
theta_tar = 60;

% Plasticity model parameters
A2d      = 10*0.00014;                           % Depression Amplitude. 10 times faster to speed up the simulations
A3p      = 10*0.00008 ;                          % Potentiation Amplitude.  10 times faster to speed up the simulations
tau_p    = 7;                                    % time constant of the voltage low-pass for the potentiation term
tau_r    = 15;                                   % time constant for the presynaptic spike low-pass
tau_d    = 10;                                   % time constant of the voltage low-pass for the depression term

% Input parameters
sigma     = 10;                                 % standard deviation of synaptic input
nu_in_max = 0.015;                              % amplitude of synaptic input
nu_in_min = 0.0001;                             % baseline of synaptic input
P         = 10;                                 % number of input patterns
p         = zeros(1,T);                         % generate random pattern locations
for n=1:NEP
    p((n-1)*EP+1:n*EP) = floor(P*rand)+1;
end

% Initialisation of gaussian inputs
ind = 1:Ni;
gau = nu_in_min+nu_in_max*exp(-(ind-Ni/2).^2/(2*sigma^2));
gau(Ni+1:2*Ni) = gau;
for jj=1:P
    mup = 1+(jj-1)*Ni/P;
    pat(:,jj) = gau(mup:mup-1+Ni)';
end

% Initialisation of the variables
E_L      = -70.6;                                % resting potential
w0       = 1.0;                                  % initial condition for w
wmax     = 3.0;                                  % max synaptic weights
w       = rand(Ni,n_neur)*1.5+0.5;               % random weights to start with
u        = E_L*ones(n_neur, 1);                  % membrane potential
u_md     = E_L*ones(n_neur, 1);                  % low pas of u for the depression term
u_mp     = E_L*ones(n_neur, 1);                  % low pas of u for the potentiation term
uy       = 0*ones(n_neur, 1);                    % voltage thresholded
u_sig_md = 0*ones(n_neur, 1);                    % low-pass of u thresholded for the depression term
u_sig_mp = 0*ones(n_neur, 1);                    % low-pass of u thresholded for the potentiation term
r       = zeros(Ni,1);                           % presyn low pass
inp      = (rand(Ni,1)<pat(:,p(1)));             % generation of the inputs
C        = 281;                                  % membrane capacitance [pF]
theta    = 0.0*ones(n_neur, 1);                  % sliding threshold, relative to resting potential
u1s      = E_L*ones(n_neur, 1);                  % keeping voltage in memory for 2ms
u1ss     = E_L*ones(n_neur, 1);                  % keeping voltage in memory for 2ms
counter  = 0*ones(n_neur, 1);                    % initial values
w_tail   = 0*ones(n_neur, 1);                    % initial values
wad      = 0*ones(n_neur, 1);                    % initial values
V_T      = -50.4*ones(n_neur, 1);                % initial values
w_save = zeros(Ni, T/100);                       % save the weights for plotting

% Iterations
for t=2:T
    inp = (rand(Ni,1)<pat(:,p(t))); % input pattern
    I = w'*inp*C*4; % currents
    [u, wad,w_tail, counter, V_T] = aEIF(u, wad,w_tail, I,counter, V_T); % voltage, aEIF neuron model
    uy = ((u-E_L) > 25.3).*(u-E_L-25.3); % threshold of voltage
    r = (1-dt/tau_r)*r+dt/tau_r*inp; % low pass of presynaptic spike
    u_mp = u1ss*dt/tau_p +(1-dt/tau_p)*u_mp; % low pass of postsynaptic voltage for the potentiation term
    u_md = u1ss*dt/tau_d +(1-dt/tau_d)*u_md;% low pass of postsynaptic voltage for the depression term
    u_sig_md = ((u_md-E_L) > 0.).*(u_md-E_L-0.); % threshold the low pass
    u_sig_mp = ((u_mp-E_L) > 0.).*(u_mp-E_L-0.);% threshold the low pass
    theta = (1-dt/tau_th)*theta+dt/tau_th*((u1ss-E_L).^2); % homeostasis
    w = w + A3p*r*(u_sig_mp.*uy)'-A2d*inp*(u_sig_md.*(theta./theta_tar))'; % weight update
    w(find(w<0)) = 0; % weight lower bound
    w(find(w>wmax)) = wmax; % weight upper bound
    w_save(:,floor(t/100)+1) = w; % save the weights for plotting
    u1s = u; % save the voltage 2ms before for the update order
    u1ss = u1s; % save the voltage 2ms before for the update order
end

% Plot
figure; imagesc(w_save)
set(gcf,'PaperUnits','centimeter ')
set(gca,'FontSize',16,'FontName','Helvetica','linewidth',2)
ylabel('neuron index')
xlabel('time')
caxis([0 3])