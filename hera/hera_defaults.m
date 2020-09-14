%--------------------------------------------------------
%HeRA default settings
%--------------------------------------------------------

hera_def.pulsefile_extension = 'puls';
hera_def.prepeak = 5000;
hera_def.prepeakdelay = 2;          %Number of samples delay in peak detection of scanner
hera_def.samplerate = 50;
hera_def.barscalefactor = .57;      %scale factor for the progressbar (percentage > size in pixels)


%Settings for CWT
hera_def.hpd = 0.01;            %High pass filter cutoff for heartrate
hera_def.lpd = 2.5;             %Low pass filter cutoff for heartrate

%Settings for tracing
hera_def.avgHeartbeat = 1;      %Initial guess for average heart rate frequency (Hz)
hera_def.ballmobilityrange = .1;%Maximum tracing ball mobility (in Hz.) per sample (.1 default for 50 Hz)
hera_def.startCutout = 2;       %Cutout #of seconds for tracing at beginning of sample
hera_def.endCutout = 2;         %Cutout #of seconds for tracing at end of sample
hera_def.hpt = 0.005;           %High pass filter cutoff for trace
hera_def.lpt = 1;               %Low pass filter cutoff for trace

%Settings for power calculations
hera_def.lowerboundLF = .05;      %Lower bound frequency (Hz) for low frequency HRV
hera_def.upperboundLF = .15;      %Upper bound frequency (Hz) for low frequency HRV
hera_def.lowerboundHF = .15;      %Lower bound frequency (Hz) for high frequency HRV
hera_def.upperboundHF = .4;       %Upper bound frequency (Hz) for high frequency HRV



