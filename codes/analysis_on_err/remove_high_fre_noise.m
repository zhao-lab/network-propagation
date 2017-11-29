
%%

openExample('dsp/RemovingHighFrequencyNoiseFromAnECGSignalExample')

openExample('signal/SignalSmoothingExample')


%%
[humidity,time] = thingSpeakRead(12397,'NumPoints',8000,'Fields',3)

% Filter type
filtertype = 'FIR';
% Sampling Frequency
Fs = 1/60;
% Filter Order
N = 3;
% Passband frequency
Fpass = 1/(24*60*60);
% Stopband frequency
Fstop = 1/(2*60*60);
% Ripple Factor and Attenuation in stop band
Rp = 0.5;
Astop = 50;

% Design the filter using dsp.LowpassFilter
LPF = dsp.LowpassFilter('SampleRate',Fs,...
                        'FilterType',filtertype,...
                        'PassbandFrequency',Fpass,...
                        'StopbandFrequency',Fstop,...
                        'PassbandRipple',Rp,...
                        'StopbandAttenuation',Astop);
                    Output = step(LPF, humidity);
                    % Change the Channel ID and the '<write API key>' to send data to your channel
thingSpeakWrite(<channelID>,Output,'Timestamps',time,'WriteKey','<write API key>');
% Only use plot in MATLAB Visualizations
plot(time,humidity,time,Output);
ylabel('Relative Humidity');
legend('Raw Data', 'Filtered Data');
                    
                    
