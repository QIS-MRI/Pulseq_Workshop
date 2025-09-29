function [spectrograms, spectrogram_rss, frequencies, times] = calculate_gradient_spectrum( ...
    seq, acoustic_resonances, use_derivative, frequency_oversampling, window_width, max_frequency, combine_mode, time_range, doPlot)

% MATLAB version of pypulseq Sequence.calculate_gradient_spectrum
% 
% Inputs:
%   seq                 - mr.Sequence object
%   acoustic_resonances - struct array with fields 'frequency' and 'bandwidth'
%   use_derivative      - logical, whether to use derivative of gradient waveform
%   frequency_oversampling - oversampling factor for FFT
%   window_width        - width of spectrogram window (s)
%   max_frequency       - maximum frequency (Hz)
%   combine_mode        - 'max', 'mean', 'rss', or 'none'
%   time_range          - [tmin, tmax] (s), or [] for full sequence
%   doPlot              - logical, whether to plot results
%
% Outputs:
%   spectrograms   - cell array of spectrograms per axis
%   spectrogram_rss- root-sum-of-squares combined spectrogram
%   frequencies    - frequency axis
%   times          - time axis (only relevant if combine_mode == 'none')

if nargin < 2 || isempty(acoustic_resonances), acoustic_resonances = []; end
if nargin < 3, use_derivative = false; end
if nargin < 4, frequency_oversampling = 3; end
if nargin < 5, window_width = 0.05; end
if nargin < 6, max_frequency = 2000; end
if nargin < 7, combine_mode = 'max'; end
if nargin < 8, time_range = []; end
if nargin < 9, doPlot = true; end

dt = seq.gradRasterTime;
%nwin = round(window_width / dt);
%nfft = round(frequency_oversampling * nwin);

% --- Get gradient waveforms ---
gw = seq.waveforms_and_times();
max_t = max(cellfun(@(g) g(1,end), gw(~cellfun('isempty',gw))));

% Time axis
if isempty(time_range)
    nt = ceil(max_t / dt);
    t = ((0:nt-1) + 0.5) * dt;
else
    tmax = min(time_range(2), max_t) - max(time_range(1),0);
    nt = ceil(tmax / dt);
    t = max(time_range(1),0) + ((0:nt-1)+0.5)*dt;
end

% Resample gradients
ng = 3;
g = zeros(ng, numel(t));
for i=1:ng
    if ~isempty(gw{i})
        g(i,:) = interp1(gw{i}(1,:), gw{i}(2,:), t, 'linear', 0);
    end
end

if use_derivative
    g = diff(g,1,2);
    t = t(1:end-1);
end

% --- Compute spectrograms ---
signal_length = size(g,2);

if signal_length <= 1
    % Sequence too short for spectrogram: return trivial outputs (prevents crash)
    warning('Sequence too short for spectral analysis (length %d). Returning trivial outputs.', signal_length);
    spectrograms = {0,0,0};
    spectrogram_rss = 0;
    frequencies = 0;
    times = 0;
    return;
end

% desired nwin from window_width (seconds)
nwin_desired = max(1, round(window_width / dt));   % at least 1
% cap to available samples and set a minimum reasonable window (e.g. 4)
nwin = min(nwin_desired, signal_length);
nwin = max(4, nwin);               % enforce a small minimum so hamming/fft behave well
nwin = min(nwin, signal_length);   % ensure still <= signal_length after min

% nfft: keep at least nwin for sensible FFT; respect frequency_oversampling
nfft = max(nwin, round(frequency_oversampling * nwin));

% overlap must be less than nwin
noverlap = min(floor(nwin/2), nwin-1);

% --- Compute spectrograms ---
spectrograms = cell(1,ng);
spectrogram_rss = 0;


for i=1:ng
    [sxx,freqs,times] = spectrogram(g(i,:), hamming(nwin), floor(nwin/2), nfft, 1/dt, 'yaxis');
    sxx = abs(sxx); % magnitude
    
    mask = freqs < max_frequency;
    sxx = sxx(mask,:);
    freqs = freqs(mask);
    
    spectrogram_rss = spectrogram_rss + sxx.^2;
    
    switch combine_mode
        case 'max'
            s = max(sxx,[],2);
        case 'mean'
            s = mean(sxx,2);
        case 'rss'
            s = sqrt(sum(sxx.^2,2));
        case 'none'
            s = sxx;
        otherwise
            error('Unknown combine_mode: %s', combine_mode);
    end
    spectrograms{i} = s;
end

spectrogram_rss = sqrt(spectrogram_rss);
if ~strcmp(combine_mode,'none')
    switch combine_mode
        case 'max'
            spectrogram_rss = max(spectrogram_rss,[],2);
        case 'mean'
            spectrogram_rss = mean(spectrogram_rss,2);
        case 'rss'
            spectrogram_rss = sqrt(sum(spectrogram_rss.^2,2));
    end
end

frequencies = freqs;

% --- Plotting ---
if doPlot
    if ~strcmp(combine_mode,'none')
        figure; hold on
        plot(frequencies, spectrograms{1});
        plot(frequencies, spectrograms{2});
        plot(frequencies, spectrograms{3});
        plot(frequencies, spectrogram_rss, 'k','LineWidth',1.5);
        xlabel('Frequency (Hz)'); ylabel('Amplitude');
        legend('X','Y','Z','RSS');
        title('Gradient Spectrum');
        
        for k=1:numel(acoustic_resonances)
            res = acoustic_resonances(k);
            xline(res.frequency,'r-');
            xline(res.frequency-res.bandwidth/2,'r--');
            xline(res.frequency+res.bandwidth/2,'r--');
        end
    else
        chans = {'X','Y','Z'};
        for i=1:ng
            figure;
            imagesc(times, frequencies, abs(spectrograms{i}));
            axis xy;
            xlabel('Time (s)'); ylabel('Frequency (Hz)');
            title(['Spectrum ' chans{i}]);
            
            for k=1:numel(acoustic_resonances)
                res = acoustic_resonances(k);
                yline(res.frequency,'r-');
                yline(res.frequency-res.bandwidth/2,'r--');
                yline(res.frequency+res.bandwidth/2,'r--');
            end
        end
    end
end

end
