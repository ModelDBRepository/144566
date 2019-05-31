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

function w = VoTri(rho, Dt, par)

% parameters of the plasticity rule
A_m=par(1);                   % amplitude for depression
A_p=par(2);                   % amplitude for potentiation
tau_p =par(3);                % time constant for voltage filter for potentiation
tetam = -70.6;                % low voltage threshold 
tetap = -45.3;                % high voltage threshold
tau_r = par(4);               % time constant low pass r [ms]
tau_d = par(5);               % low pass of the membrane pot [ms] for the depression term
 

% Protocol of Sjoestroem et al. Neuron 2001: pre and post spike trains
n = 5;                       % number of pairing 
f = 1000/rho;                % time between two pairs [ms]
l = n*f+abs(Dt)+1;           % length of the spike trains
x = zeros(1,l);              % presyn spike train
x(abs(Dt)+1:f:end-1) = 1;
y = zeros(1,l);              % postsyn spike train
y(abs(Dt)+1+Dt:f:end) = 1;

%Protocol: repetition at 0.1Hz
rho_low = 0.1;
n_rep = 15;
if rho == 0.1
    n_rep = 10;
end
x = repmat([x,zeros(1,1000/rho_low-length(x))],1,n_rep);
y = repmat([y,zeros(1,1000/rho_low-length(y))],1,n_rep);

% Parameters of the neuron model
E_L = -70.6;    %[mV]        % resting pot

% Initialization
l = length(x);               % simulation length
I_ext = zeros(1,l);          % extra current injected
I_s = y*1000000;             % current injected for a postsynaptic spike induction
w_tail = zeros(1,l);         % current for spike after-depolarisation
V_T = -50.4;                 % adaptive threshold
I_tot = zeros(1,l);          % total current
wad = zeros(1,l);            % adaptation
u = E_L*ones(1,l);           % membrane potential
u_md = E_L*ones(1,l);        % filtered membrane potential for the depression term
u_mp = E_L*ones(1,l);        % filtered membrane potential for the potentiation term 
r = zeros(size(x));          % low-pass the presynaptic spike-train, x
w = 0.5;                     % initialization of weights
counter = 0;                 % trick to count how long is the spike - here spikes are forced to be 2ms long

% Main loop over time
for t = 4:l
    I_tot(t) =  I_s(t)+I_ext(t); % current
    [u(t), wad(t), w_tail(t), counter, V_T] = aEIF(u(t-1), wad(t-1), w_tail(t-1),I_tot(t), counter, V_T); % voltage with the aEIF neuron model
    u_md(t+1) = u(t)/tau_d +(1-(1/tau_d))*u_md(t); % low-pass of voltage for depression
    u_mp(t+1) = u(t)/tau_p +(1-(1/tau_p))*u_mp(t);% low-pass of voltage for potentiation
    r(t+1) = x(t)/tau_r +(1-(1/tau_r))*r(t); % low-pass of presynaptic spike
    u_sig = (u(t) > tetap)*(u(t)-tetap); % voltage threhold
    u_md_sig = ((u_md(t-3)-tetam) > 0)*(u_md(t-3)-tetam); % low-pass threshold
    u_mp_sig = ((u_mp(t-3)-tetam) > 0)*(u_mp(t-3)-tetam); % low-pass threshold
    w=  w- A_m*x(t)*u_md_sig+ A_p*u_sig*r(t)*u_mp_sig; % weight updates
end