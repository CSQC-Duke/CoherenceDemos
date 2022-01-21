clear

c = 2.9979 * 10^-5;                % speed of light[cm/fs]

%%%%%%%%%%%%%%%%%%%%%
%% Tunable parameters 
%%%%%%%%%%%%%%%%%%%%%

lambda0_main = 3000;                 % pulse 1: central wavelength [nm]
w0_main = 1/(lambda0_main/(10^7));   % pulse 1: carrier frequency [cm-1]
tau_p_main = 20;                     % pulse 1: temporal width of main pulse [fs]

lambda0_ref = 4000;                 % pulse 2: central wavelength [nm]
w0_ref = 1/(lambda0_ref/(10^7));    % pulse 2: carrier frequency [cm-1]
tau_p_ref = 20;                     % pulse 2: temporal width of main pulse [fs]

%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate main pulse 
%%%%%%%%%%%%%%%%%%%%%%%%

dt = (tau_p_main+tau_p_ref)/100;      % time resolution of time vector (fs)
trange = 10*(tau_p_main+tau_p_ref);   % range of main time vector (-min and max value)
trange_delay = 0.5*trange;            % range of delay time vector (-min and max value)
a_main = 0;                           % unitless number. When 0: no chirp
t0_main = 0;                          % time center for pulse [fs]
t= -trange:dt:trange;                               % time vector for pulse
Emain_t = genpulse(t, t0_main, w0_main, a_main, tau_p_main);  % create main pulse centered at time 0
dw_main = 1 / (c*dt*numel(t));                         % defining the frequency step [cm-1]
w_main = (-numel(t)*dw_main)/2:dw_main:...
    (numel(t)*dw_main)/2-dw_main;                        % frequency axis [cm-1%]
w_main = -w_main;
Emain_w = fftshift(fft(Emain_t));               % the spectrum of the main pulse

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up delay vectors 
%%%%%%%%%%%%%%%%%%%%%%%%%

t_delay = -trange_delay:dt:trange_delay;          % delay-time vector  [fs%]
dw_delay = 1 / (c*dt*numel(t_delay));                   % delay frequency step [cm-1%]
w_delay = (-numel(t_delay)*dw_delay)/2:dw_delay:...
    (numel(t_delay)*dw_delay)/2-dw_delay;               % frequency axis [cm-1%]
w_delay = -w_delay;

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate and sweep reference pulse 
%%%%%%%%%%%%%%%%%%%%%%%%%

%Field and intensity correlations are based off the equations in the
%corresponding sections of the following webpage: 
%en.wikipedia.org/wiki/Optical_autocorrelation

a_ref = 0;                          % unitless number. When 0: no chirp
I_FCint = zeros(1,numel(t_delay));
I_ICint = zeros(size(I_FCint));
FieldCorrelation = zeros(numel(t_delay),numel(t));
IntensityCorrelation = zeros(size(FieldCorrelation));

for i = 1:numel(t_delay)
    
    Eref_t = genpulse(t, t_delay(i), w0_ref, a_ref, tau_p_ref);  % delayed reference pulse (delayed using movable mirror)
    FieldCorrelation(i,:) = abs(Emain_t+Eref_t).^2;     % the total intensity on the detector vs. time
    I_FCint(i) = trapz(t, FieldCorrelation(i,:));       % the integrated detector intensity (since the detector is slow it basically integrates over the intensity of the pulses)
    
    IntensityCorrelation(i,:) = abs(Emain_t.*Eref_t).^2;
    I_ICint(i) = trapz(t, IntensityCorrelation(i,:));

end

I_FCint = I_FCint-50;       % remove background

%calculate frequency-domain correlation spectra
I_FCw = fftshift(fft(I_FCint));               % the spectrum of the field correlation 
I_ICw = fftshift(fft(I_ICint));               % the spectrum of the intensity correlation

Eref_t = genpulse(t, t0_main, w0_ref, a_ref, tau_p_ref); % delayed reference pulse (delayed using movable mirror)
Eref_w = fftshift(fft(Eref_t));               % the spectrum of the main pulse

%%%%%%%%%%%
%% Plotting
%%%%%%%%%%%

fig = figure;
subplot(3,1,1)
hold on
plot(t, real(Emain_t), 'LineWidth', 1, 'Color', 'r')
plot(t, real(Eref_t), 'LineWidth', 1, 'Color', 'g')
xlabel('t (fs)');
title('Main and reference pulse profiles')
legend('Pulse 1 field', 'Pulse 2 field')
box on
set(gca,'FontSize',14)
hold off

subplot(3,1,2)
hold on
plot(t_delay, I_FCint./max(abs(I_FCint)), 'LineWidth', 1, 'Color', 'b')
plot(t_delay, I_ICint./max(abs(I_ICint)), 'LineWidth', 1, 'Color', 'k');
xlabel('\Deltat (fs)');
legend('Field correlation', 'Intensity correlation')
title('Integrated detection')
box on
set(gca,'FontSize',14)
hold off

subplot(3,1,3);
hold on
plot(w_main./1000, abs(Emain_w)./max(abs(Emain_w)),'LineWidth', 1, 'Color', 'r')
plot(w_main./1000, abs(Eref_w)./max(abs(Eref_w)),'LineWidth', 1, 'Color', 'g')
plot(w_delay./1000, abs(I_FCw)./max(abs(I_FCw)),'LineWidth', 1, 'Color', 'b')
xlabel('\omega/2\pic (10^3 cm-1)');
title('Frequency domain')
legend('Pulse 1 Spectrum', 'Pulse 2 Spectrum', 'Field Correlation FFT')
xlim([(w0_main-3000)./1000, (w0_main+3000)./1000]);
set(gca,'FontSize',14)
box on

disp('done!');