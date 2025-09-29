clear all; clc

%% CREATE SEQUENCE

% set system limits (slew rate 130 and max_grad 30 work on Prisma)
sys = mr.opts('MaxGrad', 24, 'GradUnit', 'mT/m', ...
    'MaxSlew', 100, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 30e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'B0', 3);

seq = mr.Sequence();     % Create a new Pulseq sequence

%% PARAM

Ndummy=10;                      % number of dummy scans

sliceThickness=3e-3;           % Slice Thikness
Nspokes = 603;                 % Number of Spokes
Nx      = 200;                 % Samples per spoke
FOV     = 0.300;               % Field of view (m)

dwell = 1e-5;                  % Dwell Time [sec]
duration180 = 10e-3;           % 180° RF Duration [sec]
duration90 = 10e-3;            % 90° RF Duration [sec]
timeBwProduct90 = 4;           % Bandwidth
apodization = 0.4;             % Apodization

TE  = 50e-3;                   % TE [sec]
TR = 150e-3;                   % TR [sec]

adcDur  = Nx*dwell;
dk      = 1/FOV;
kWidth  = Nx*dk;
GA = ((pi/2)* (3-sqrt(5))) / pi * 180;
radian = mod((0:GA:(Nspokes-1)*GA)*pi/180,2*pi);
%% CREATE GRAD

gr = mr.makeTrapezoid('x','system',sys,  'FlatArea', kWidth, 'FlatTime', adcDur);               % Readout Gradient
grSpoil=mr.makeTrapezoid('x','Area',0.5*Nx*dk,'system',sys);                                    % Spoiler Gradient
grPre = mr.makeTrapezoid('x','system',sys, 'Area', (kWidth/2 + dk/2), 'Duration', adcDur/2);   % Prephasing Gradient

%% CREATE ADC

adc = mr.makeAdc(Nx, 'Dwell', dwell,  'Delay', gr.riseTime,'system', sys);

%% CREATE THE 90 DEGREES RF

[rf90, gz90,gz90_reph] = mr.makeSincPulse(deg2rad(90),'sys',sys,'duration', duration90,'timeBwProduct', timeBwProduct90,'apodization', apodization,...
    'sliceThickness', sliceThickness,'dwell', dwell,'use','excitation');

%% CREATE THE 180 DEGREES RF

phaseOffset = pi/2;
[rf180, gz180] = mr.makeSincPulse(deg2rad(180),'sys',sys,'duration', duration180,'timeBwProduct', timeBwProduct90,'apodization', apodization,...
    'sliceThickness', sliceThickness,'phaseOffset',phaseOffset,'dwell', dwell,'use','refocusing');
%% ADD BLOCKS

delayTE1=TE/2-(mr.calcDuration(gz90)-mr.calcRfCenter(rf90)-rf90.delay)-rf180.delay-mr.calcRfCenter(rf180)-mr.calcDuration(grPre) - mr.calcDuration(gz90_reph);
delayTE2=TE/2-mr.calcDuration(rf180)+rf180.delay+mr.calcRfCenter(rf180)-adc.delay-adcDur/2;
delayTR=TR-mr.calcDuration(gz90)-mr.calcDuration(grPre)-delayTE1-mr.calcDuration(rf180)-delayTE2-mr.calcDuration(gr);

assert(delayTE1>=0);
assert(delayTR>=0);

for i=(1-Ndummy):Nspokes

    % adding seq blocks
    seq.addBlock(rf90,gz90);
    seq.addBlock(gz90_reph);

    if i>0
        seq.addBlock(mr.rotate('z',radian(i),grPre));
    else
        seq.addBlock(mr.rotate('z',radian(1),grPre));
    end

    seq.addBlock(mr.makeDelay(delayTE1));
    seq.addBlock(rf180);
    seq.addBlock(mr.makeDelay(delayTE2));

    if (i>0)
        seq.addBlock(mr.rotate('z',radian(i),adc,gr));

    else
        seq.addBlock(mr.rotate('z',radian(1),gr));
    end

    seq.addBlock(mr.makeDelay(delayTR));
end


%% PLOT

seq.plot();

%% TIMINGS CHECK

% check whether the timing of the sequence is compatible with the scanner
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% PLOT TRAJ
% calculate k-space but only use it to check timing
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay',[0 0 0]*1e-6); % play with anisotropic trajectory delays -- zoom in to see the trouble ;-)


% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kz-axis
title('k-vector components as functions of time'); xlabel('time /s'); ylabel('k-component /m^-^1');
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('2D k-space trajectory'); xlabel('k_x /m^-^1'); ylabel('k_y /m^-^1');
%% ACCELLERATION FACTOR
nBin = 3;
ReqSpokes = (pi/2)*Nx*nBin;


fprintf('Sequence duration: %.2f s\n', seq.duration);
fprintf('Required spokes: %.0f\n', ReqSpokes);


%% SAVE SEQUENCE


fname = sprintf('2DRad_GA_nSpoke%d_SliceThikness%dmm_TE%dms_TR%.1dms',Nspokes,sliceThickness*1000,TE*1000,TR*1000);
seq.write(fname);
fprintf('Sequence saved as %s\n', fname);
