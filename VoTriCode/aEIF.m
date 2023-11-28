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

% This is the code to simulate the neuron model

function [u, w,w_tail,counter,V_T] = aEIF(u,w,w_tail,I,counter,V_T)

% AEIF: Simulate adex model with spike after depolarization current and
% adaptive threshold
%
% Preconditions:
%  u            Membrane potential at time t
%  w            Adaptation variable
%  w_tail       Current for spike afterdepolarization
%  V_T          Adaptive threshold 
%  I            Input Current
% counter       Counter to force the spike to be clamped at the high value
%               for 2ms


% Model parameters
th = 20;        % [mV] spike threshold
C = 281;        % [pF] membrane capacitance
g_L = 30;       % [nS] membrane conductance
E_L = -70.6;    % [mV] resting voltage
VT_rest = -50.4;% [mV] resetting voltage
Delta_T = 2;    % [mV] exponential parameters
tau_w = 144;    % [ms] time constant for adaptation variable w
a = 4;          % [nS] adaptation coupling constant
b = 0.0805;     % [nA] spike triggered adaptation
w_jump = 400;   % spike after depolarisation
tau_wtail = 40; % [ms] time constant for spike after depolarisation
tau_VT = 50;    % [ms] time constant for VT
VT_jump = 20;   % adaptive threshold

if counter ==2          % trick to force the spike to be 2ms long
    u = E_L+15+6.0984;  % resolution trick (simulation of the spike at a fine resolution - see below)
    w = w+b;
    w_tail = w_jump;
    counter = 0;
    V_T = VT_jump+VT_rest;
end

% Updates of the variables for the aEIF
udot = 1/C*(-g_L*(u-E_L) + g_L*Delta_T*exp((u-V_T)/Delta_T) - w +w_tail+ I);
wdot = 1/tau_w*(a*(u-E_L) - w);
u= u + udot;
w = w + wdot;
w_tail = w_tail-w_tail/tau_wtail;
V_T = VT_rest/tau_VT+(1-1/tau_VT)*V_T;

if counter == 1
    counter = 2;
    u = 29.4+3.462; % resolution trick (simulation of the spike at a fine resolution - see below)
    w = w-wdot;
end

if (u>th && counter ==0) % threshold condition for the spike
    u = 29.4;
    counter = 1;
end

% numerical trick for the aEIF model: I simulated, once and for all, the spike upswing and 
% integrated to know what is the integral of the spike, with high precision. Then I used this number in the simulation and 
% I clamped the spike for 2ms at the appropriated calculated value (this is to speed up the simulations 
% since network simulation is really time consuming). Since my time step in my simulation is 1ms, 
% I wait 3ms (the spike length (2ms) plus 1 time step (1ms)) before I read the filtered version of 
% the voltage (we want to read the value of the voltage trace before this spike). 
