function [Epulse] = genpulse(t, t0, carrier, chirp, duration)
% input parameters:
% t: time vector
% t0: time displacement
% lambda0: central wavelength (m)
% a: chirp parameter (unitless)
% tau: pulse temporal width

c = 2.9979 * 10^-5;       %speed of light[cm/fs]

Gaussian = @(t0, chirp, duration) exp(-((1-1i*chirp)/duration^2)*(t-t0).^2); % generic complex Gaussian (time = vector;  a,tau = scalar parameters) (see Eq. 22.1-11 in Saleh & Teich)
expPhase = @(t0) exp(-1i*2*pi*carrier*c*(t-t0)); % fast component of wave (depends on central wavelength)
osc = expPhase(t0);
env = Gaussian(t0, chirp, duration);
Epulse = osc.*env; % pulse electric field

end