% This file will run a fast pns and forbidden resonances check and print a
% png file. If the asc file from the scanner is not found, the pns check is
% skipped.  

% USAGE:
% Include the pulseq-master path with subfolders in your matlab path
% Run the safety check as: safety_check('my_sequence.seq')
% If you have the asc file, place it in the same folder as this file and
% change the name below to match that of your file.


function safety_check(seq_file,scanner)

    % Gradient system file (must be in the same folder)
    gradFile = 'MP_GPA_K2309_2250V_951A_AS82.asc';
    
    f = figure('Color','w');%("Position",[680   558   990   420]); %Adjust accordingly if needed
    t = tiledlayout(1,2);

    % Load sequence
    seq = mr.Sequence();
    seq.read(seq_file);

    % PNS check
    check_PNS(seq, gradFile,t);

    % Acoustic resonance check
    check_acoustic_resonances(seq, gradFile, t, scanner);

    % Add a global title with the hash (if found)
    title(t, sprintf('Hash: %s', seq.signatureValue), ...
        'FontWeight','bold','Interpreter','none');


    [~,name] = fileparts(seq_file);
    saveas(f,[name '_safety_check.png']);

    disp('Safety check finished')
end

function check_PNS(seq, gradFile,t)
    nexttile(t);
    
    if exist(gradFile, 'file')
        [pns_ok, pns_n, pns_comp, ~] = seq.calcPNS(gradFile,false);

        safe_plot(pns_comp'*100, seq.gradRasterTime);
    
        if pns_ok
            fprintf('PNS check passed: max PNS level = %.1f %%\n', 100*max(pns_n));
        else
            fprintf('PNS check failed: max PNS level = %.1f %%\n', 100*max(pns_n));
            error('PNS check failed!');
        end
    else
        warning('Gradient file not found. Skipping PNS check.');
        
        % Create an empty plot to keep layout consistent
        plot(0,0,'w'); axis off;
        title('PNS Check (Skipped)');
    end
end

function check_acoustic_resonances(seq, gradFile,t,scanner)

    if strcmp(scanner,'3T')
        DEFAULT_ACOUSTIC_RESONANCES = [ ...
            struct('frequency', 590,  'bandwidth', 100); ...
            struct('frequency', 1140, 'bandwidth', 220) ...
        ];
    elseif strcmp(scanner,'7T')
        DEFAULT_ACOUSTIC_RESONANCES = [ ...
            struct('frequency', 1100,  'bandwidth', 300); ...
            struct('frequency', 550, 'bandwidth', 100) ...
        ];
    end

    if exist(gradFile, 'file')
        asc = mr.Siemens.readasc(gradFile);
        freqs = asc.asGPAParameters.sGCParameters.aflAcousticResonanceFrequency;
        bws   = asc.asGPAParameters.sGCParameters.aflAcousticResonanceBandwidth;
        
        acoustic_resonances = [];
        for k = 1:numel(freqs)
            if freqs(k) ~= 0
                acoustic_resonances(end+1).frequency = freqs(k);
                acoustic_resonances(end).bandwidth = bws(k);
            end
        end
    else
        warning('Gradient file not found. Using default acoustic resonances.');
        acoustic_resonances = DEFAULT_ACOUSTIC_RESONANCES;
    end

    [spectrograms, spectrogram_rss, frequencies] = calculate_gradient_spectrum(seq, acoustic_resonances, false, 10, 0.02, 2000, 'max', [], false);
    
    nexttile(t);
    hold all;
    plot(frequencies, spectrograms{1},'DisplayName','x');
    plot(frequencies, spectrograms{2},'DisplayName','y');
    plot(frequencies, spectrograms{3},'DisplayName','z');
    plot(frequencies, spectrogram_rss,'r','LineWidth',1,'DisplayName','rss');

    % Add acoustic resonances
    for k=1:numel(acoustic_resonances)
        res = acoustic_resonances(k);
        xline(res.frequency,'k-','LineWidth',2,'HandleVisibility','off');
        xline(res.frequency-res.bandwidth/2,'k--','LineWidth',2,'HandleVisibility','off');
        xline(res.frequency+res.bandwidth/2,'k--','LineWidth',2,'HandleVisibility','off');
    end
    xlabel('Frequency (Hz)'); ylabel('Amplitude');
    box on;
    legend();
    title('Acoustic Resonances');
end

