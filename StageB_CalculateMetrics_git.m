% =========================================================================
% Stage B: Calculate Mean Sulci Metrics (Single Subject Batch Version) V4.4 FINAL
% =========================================================================
% Description:
% Loads pre-computed volumes from Stage A for a single subject.
% Calculates mean metrics based on reference code (localmax_mod) and
% local maxima overlap counts per sulcus.
% Saves mean+count backup per subject. Appends ONLY mean metrics (with string ID)
% to central file using fprintf.
% Assumes Stage A saved files with pre-permutation and this script uses
% loadRawImage with post-permutation.
%
% Called by: Batch script (e.g., qsub) which sets JOB_SUBJECT_ID env var.
% Requires: Outputs from StageA_GenerateVolumes.m
% =========================================================================

function StageB_CalculateMetrics() % Wrapped in function

    fprintf('--- Starting StageB_CalculateMetrics ---\n');
    startTime = tic;

    % --- Configuration ---
    config = struct();
    % Paths (VERIFY THESE)
    config.rootFolder = fullfile('/Users/yourname/Documents/NeuroData/FS_Outputs');
    config.outputSummaryFile = fullfile(config.rootFolder, 'summary_sulci_metrics_MEAN_v3.csv'); % CENTRAL file (Mean only, V2 to ensure clean start)
    config.imageSize = [256, 256, 256]; % MUST match Stage A
    config.maxExpectedLabels = 68; % Max columns for summary file header/padding (adjust if needed, e.g., 68)
    config.verbose = true; % Keep true for logs
    config.showPlots = false; % SET FALSE FOR BATCH JOBS
    config.overwriteSummary = false; % APPEND mode default for batch runs

    % Define log paths
    config.logFolder = fullfile(config.rootFolder, 'logs');
    config.errorLogFile = fullfile(config.logFolder, 'error_log.txt');
    if ~exist(config.logFolder, 'dir'); try mkdir(config.logFolder); catch; fprintf('Warning: Could not create log folder %s\n', config.logFolder); end; end

    % Paths needed by helpers if not running compiled
    if ~isdeployed
        config.helperScriptsPath = '/Users/yourname/Scripts_matlab';
        config.niftiToolsPath = '/Users/yourname/Scripts_matlab/Plugins_ML/NIfTI_20140122'; % Needed for niftiinfo in fallback
        addPaths(config); % Add paths if needed
    end

    % --- Get Subject ID ---
    datasetID = getenv('JOB_SUBJECT_ID');
    if isempty(datasetID); error('Environment variable JOB_SUBJECT_ID not set.'); end
    fprintf('\n--- Processing Dataset for Metrics: %s ---\n', datasetID);
    datasetFolder = fullfile(resolveTilde(config.rootFolder), datasetID);
    if ~exist(datasetFolder, 'dir'); error('Dataset folder not found: %s', datasetFolder); end
    
    paths = defineDatasetPaths(datasetFolder); % Define paths
    imSize = config.imageSize;
    sulciMetrics = []; % Initialize main output
    sulciOverlapCounts = []; % Initialize secondary output

    try
        % --- Initialize Central Summary File Header (less critical in append mode) ---
        initializeSummaryFile_StageB(config); % Ensure header exists/is correct

        % --- Step 1: Load Pre-Computed Volumes ---
        fprintf('  Loading pre-computed volumes...\n');

        labelFile = paths.sulciDilatedRaw; if ~exist(labelFile,'file'); error('Dilated labels file missing: %s', labelFile); end
        SL2 = loadRawImage(labelFile, imSize, 'int32', 'ieee-le'); % Loads & permutes

        localMaxFile = paths.localMaxRaw; if ~exist(localMaxFile,'file'); error('Local Maxima file missing: %s', localMaxFile); end
        VoxelMat3 = loadRawImage(localMaxFile, imSize, 'uint8', 'ieee-le'); % Loads & permutes
        VoxelMat3 = VoxelMat3/255;
        distMapFile = paths.distMapRaw; if ~exist(distMapFile,'file'); error('Distance Map file missing: %s', distMapFile); end
        %VoxelMat5 = loadRawImage(distMapFile, imSize, 'float', 'ieee-le'); % Loads & permutes
        VoxelMat5 = multibandread(distMapFile, [imSize(1), imSize(2), imSize(3)], 'float', 0, 'bsq', 'ieee-le');
        VoxelMat5=permute(VoxelMat5, [2,1,3]);
	fprintf('      Loaded Dilated Labels (SL2): Size=%s, NNZ=%d, Max=%d\n', mat2str(size(SL2)), nnz(SL2), max(SL2(:)));
        fprintf('      Loaded Local Max (VoxelMat3): Size=%s, NNZ=%d, Max=%d\n', mat2str(size(VoxelMat3)), nnz(VoxelMat3), max(VoxelMat3(:)));
        fprintf('      Loaded Dist Map (VoxelMat5): Size=%s, NNZ=%d\n', mat2str(size(VoxelMat5)), nnz(VoxelMat5));

        % --- Alignment Check ---
        if ~isequal(size(SL2), size(VoxelMat3), size(VoxelMat5)); error('Dimension mismatch after loading!'); end
        fprintf('      Dimensions consistent after loading.\n');
        % Visualization disabled for batch run

        % --- Step 2: Calculate `localmax_mod` ---
        fprintf('  Calculating modified local maxima map (localmax_mod)...\n');
        VoxelMat4 = VoxelMat3;
        localmax_mod = VoxelMat4 .* VoxelMat5.*2;
        fprintf('      Calculated localmax_mod map: NNZ=%d\n', nnz(localmax_mod));
        clear VoxelMat4 VoxelMat5; % Clear large inputs early (keep VoxelMat3 for count)

        % --- Step 3: Calculate Overlap Counts AND Mean Metrics ---
        fprintf('  Calculating overlap counts and mean metrics...\n');

        % --- Calculate Overlap Counts (SL2 vs VoxelMat3 > 0) ---
        maxActualLabel = max(SL2(:));
        if isempty(maxActualLabel) || maxActualLabel == 0; maxActualLabel = config.maxExpectedLabels; warning('No non-zero labels found in SL2 map for %s. Using default max label %d.', datasetID, config.maxExpectedLabels); end
        sulciOverlapCounts = zeros(1, maxActualLabel); % Initialize counts array
        localMaxMask = VoxelMat3 > 0;
        fprintf('      Counting overlap for labels 1 to %d...\n', maxActualLabel);
        for j = 1:maxActualLabel
            labelMask = (SL2 == j);
            if any(labelMask(:)); intersection = labelMask & localMaxMask; sulciOverlapCounts(j) = nnz(intersection); end
        end
        fprintf('      Finished overlap count loop. Total overlaps = %d.\n', sum(sulciOverlapCounts));
        clear localMaxMask labelMask intersection VoxelMat3; % Clear VoxelMat3 now

        % --- Calculate Mean ValueMap using infovol ---
        fprintf('  Calculating mean ValueMap using infovol method...\n');
        [x, y, z] = meshgrid(1:imSize(1), 1:imSize(2), 1:imSize(3));
        infovol = [x(:), y(:), z(:),(SL2(:)), localmax_mod(:)];
        clear x y z SL2 localmax_mod; % Clear large inputs

        
        fprintf('      Filtered infovol size [ %d x 5 ]\n', size(infovol, 1));

        % Initialize sulciMetrics (Means) array based on maxActualLabel
        sulciMetrics = zeros(1, maxActualLabel);

        if ~isempty(infovol)
            maxLabelInFiltered = ceil(max(infovol(:, 4)));
            if maxLabelInFiltered > maxActualLabel; sulciMetrics(maxLabelInFiltered) = 0; end % Ensure size if needed

            fprintf('      Calculating mean(localmax_mod) for labels present...\n');
            presentLabels = unique(infovol(:,4));
            for j = presentLabels'
                 if j > 0 && j <= length(sulciMetrics) % Bounds check
                    values_j = infovol(infovol(:, 4) == j, 5);
                    if ~isempty(values_j); sulciMetrics(j) = mean(nonzeros(values_j), 'omitnan'); end
                 else
                    fprintf('      Warning: Label %d from filtered infovol is out of bounds (maxLabel=%d).\n', j, maxActualLabel);
                 end
            end
            fprintf('      Finished mean metric calculation loop.\n');
        else
             warning('Infovol empty after filtering. Mean metrics will be all NaN for %s.', datasetID);
        end
        clear infovol indices_to_remove values_j presentLabels maxLabelInFiltered maxActualLabel;

        fprintf('      Final calculated Mean Metrics size: %s\n', mat2str(size(sulciMetrics)));
        fprintf('      Final calculated Overlap Counts size: %s\n', mat2str(size(sulciOverlapCounts)));

        % --- Step 4: Save Per-Subject Backup File (Means AND Counts) ---
        fprintf('  Saving per-subject backup file...\n');
        try
             maxLabelCalculated = max(length(sulciMetrics), length(sulciOverlapCounts));
             if maxLabelCalculated > 0 && (any(~isnan(sulciMetrics)) || any(sulciOverlapCounts > 0))
                 outputFilename = fullfile(paths.datasetFolder, sprintf('%s_sulci_metrics_counts_backup.txt', datasetID));
                 % Pad arrays to maxLabelCalculated if needed
                 if length(sulciMetrics) < maxLabelCalculated; sulciMetrics(maxLabelCalculated) = 0; end
                 if length(sulciOverlapCounts) < maxLabelCalculated; sulciOverlapCounts(maxLabelCalculated) = 0; end
                 % Prepare table
                 labels = (1:maxLabelCalculated)';
                 means = sulciMetrics(1:maxLabelCalculated)';
                 counts = sulciOverlapCounts(1:maxLabelCalculated)';
                 outputTable = table(labels, means, counts, 'VariableNames', {'Label', 'Mean_Value', 'LocalMax_Count'});
                 writetable(outputTable, outputFilename, 'Delimiter', '\t');
                 fprintf('    Per-subject metrics and counts saved to: %s\n', outputFilename);
             else
                 fprintf('    Skipping per-subject save as no valid metrics/counts were calculated.\n');
             end
             clear outputFilename maxLabelCalculated labels means counts outputTable;
        catch ME_write_subj
            warning('Failed to write per-subject backup file: %s', ME_write_subj.message);
            logError(ME_write_subj, config.errorLogFile, datasetID, 'Stage B - Save Subject Backup');
        end
        % --- End Per-Subject Backup ---

        % --- Step 5: Append Means (with String Subject ID) to Central Summary File ---
        if ~isempty(sulciMetrics) && any(~isnan(sulciMetrics)) % Check if we have any valid means
            % datasetID variable already holds the subject name string (e.g., 'EHBS0004')

            % Pad/truncate MEANS array consistently to match header
            maxExpectedLabels = config.maxExpectedLabels; % Use value from config
            if length(sulciMetrics) < maxExpectedLabels; sulciMetrics(maxExpectedLabels) = 0; % Pad with NaN
            elseif length(sulciMetrics) > maxExpectedLabels; sulciMetrics = sulciMetrics(1:maxExpectedLabels); end

            % Prepare output data - ONLY the means
            outputData = sulciMetrics(1:maxExpectedLabels); % Row vector of means

            % --- Use fprintf for mixed string/numeric output ---
            outputFilePath = resolveTilde(config.outputSummaryFile);
            fid = -1; % Initialize file ID
            try
                fid = fopen(outputFilePath, 'a'); % Open for appending ('a')
                if fid == -1; error('Cannot open summary file %s for appending.', outputFilePath); end
                % Create format spec: string, followed by N floats, newline
                formatSpec = ['%s', repmat(',%.8f', 1, length(outputData)), '\n'];
                % Write the subject ID string and the numeric data vector
                fprintf(fid, formatSpec, datasetID, outputData);
                fclose(fid); % Close the file
                fid = -1; % Reset file ID
                fprintf('  Mean metrics for %s appended to central summary file.\n', datasetID);
            catch ME_write_summary
                 warning('Failed to append to central summary file %s: %s', outputFilePath, ME_write_summary.message);
                 if fid ~= -1; try fclose(fid); catch; end; end % Ensure file is closed on error
                 logError(ME_write_summary, config.errorLogFile, datasetID, 'Stage B - Append Summary');
            end
            % --- End fprintf block ---

        else
             fprintf('  No valid mean metrics generated for %s, nothing appended to central summary.\n', datasetID);
        end
        % Clear variables for this iteration
        clear sulciMetrics sulciOverlapCounts outputData numericID numericIDstr outputFilePath fid formatSpec;


    catch ME % Catch errors for this dataset in Stage B
         totalTime = toc(startTime);
         fprintf(2,'\n!!! ERROR Stage B %s after %.2f sec: %s !!!\n', datasetID, totalTime, ME.message); % Use stderr
         logError(ME, config.errorLogFile, datasetID, 'Stage B');
         exit(1); % Exit with error status
    end

    totalTime = toc(startTime);
    fprintf('--- Stage B processing finished for %s (%.2f seconds) ---\n', datasetID, totalTime);
    exit(0); % Exit with success status

end % <<<< ****** END for function StageB_CalculateMetrics ******


% =========================================================================
% Stage B Helper Function Definitions
% =========================================================================
% --- I/O Function (Permute After Load) ---
function volume = loadRawImage(rawFile, imSize, dataType, byteOrder); if ~exist(rawFile,'file'); error('Raw file not found: %s', rawFile); end; fprintf('      Loading raw file: %s (Type: %s)\n', rawFile, dataType); try volume_temp = multibandread(rawFile, [imSize(1), imSize(2), imSize(3)], dataType, 0, 'bsq', byteOrder); fprintf('      Applying permute([2 1 3]) after loading raw file.\n'); volume = permute(volume_temp, [2, 1, 3]); clear volume_temp; catch ME; error('Failed load/permute raw %s: %s', rawFile, ME.message); end; end
% --- Setup Functions ---
function initializeSummaryFile_StageB(config) % Mean Only Header
    summaryFile = resolveTilde(config.outputSummaryFile);
    if config.overwriteSummary && exist(summaryFile, 'file'); fprintf('Deleting existing summary file: %s\n', summaryFile); try delete(summaryFile); catch ME_del; warning('Could not delete summary file: %s', ME_del.message); end; end
    if ~exist(summaryFile, 'file')
        fprintf('Summary file not found. Creating file with header (Means Only): %s\n', summaryFile);
        fid = fopen(summaryFile, 'a'); if fid == -1; error('Cannot open summary file: %s', summaryFile); end
        maxExpectedLabels = config.maxExpectedLabels; % Use config value
        headerLabels = sprintf('Sulcus_%d,', 1:maxExpectedLabels); % Mean only header
        header = ['DatasetID,' headerLabels(1:end-1)];
        fprintf(fid, '%s\n', header); fclose(fid);
     else
         fprintf('Appending mean metrics to existing summary file: %s\n', summaryFile);
    end
end
function logError(ME, logFilePath, datasetID, stage); resolvedLogPath = resolveTilde(logFilePath); try fid = fopen(resolvedLogPath, 'a'); if fid == -1; fprintf(2,'!!! Cannot open log file: %s !!!\n', resolvedLogPath); return; end; timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS'); fprintf(fid, 'Timestamp: %s\nDataset: %s\nStage: %s\n', timestamp, datasetID, stage); fprintf(fid, 'Error ID: %s\nError Msg: %s\n', ME.identifier, ME.message); if ~isempty(ME.stack); fprintf(fid, 'Stack Trace:\n'); for k = 1:length(ME.stack); fprintf(fid, '  File: %s, Name: %s, Line: %d\n', ME.stack(k).file, ME.stack(k).name, ME.stack(k).line); end; end; fprintf(fid, '----------------------------------------\n'); fclose(fid); catch logEx; fprintf(2,'!!! FAILED TO WRITE TO ERROR LOG FILE (%s): %s !!!\n', resolvedLogPath, logEx.message); end; end
function resolvedPath = resolveTilde(filePath); if ispc || isempty(filePath) || ~ischar(filePath); resolvedPath = filePath; return; end; if startsWith(filePath, '~'); try; [status, homeDir] = system('echo $HOME'); if status == 0 && ~isempty(strtrim(homeDir)); homeDir = strtrim(homeDir); if strcmp(filePath, '~'); resolvedPath = homeDir; elseif startsWith(filePath, '~/'); resolvedPath = fullfile(homeDir, filePath(3:end)); else; resolvedPath = fullfile(homeDir, filePath(2:end)); end; else; resolvedPath = filePath; end; catch; resolvedPath = filePath; end; else; resolvedPath = filePath; end; end
function base = basename(filepath); [~, name, ext] = fileparts(filepath); base = [name ext]; end
% --- Path Definition Function (Simplified for Stage B) ---
function paths = defineDatasetPaths(datasetFolder)
    paths = struct(); paths.datasetFolder = resolveTilde(datasetFolder);
    paths.distMapRaw = fullfile(paths.datasetFolder, 'Brain_distmap_anisotropic.raw');
    paths.localMaxRaw = fullfile(paths.datasetFolder, 'localmax_Brain_distmap_anisotropic.raw');
    paths.sulciDilatedRaw = fullfile(paths.datasetFolder, 'sulci_dil.raw');
    paths.combinedSulciNii = fullfile(paths.datasetFolder, 'Brain_sulci.nii'); % Needed for fallback maxLabel
end
% --- Visualization Function (Plots Disabled in Batch) ---
function visualizeMiddleSlice(volume, titleText, colorMapName); if nargin < 3; colorMapName = 'gray'; end; try; if isempty(volume) || ndims(volume) ~= 3; return; end; middleSliceIdx = round(size(volume, 3) / 2); if islogical(volume); volume = uint8(volume); elseif ~isfloat(volume); volume = single(volume); end; middleSlice = volume(:, :, middleSliceIdx); figHandle=figure('Name', titleText, 'NumberTitle', 'off','Visible','off'); imagesc(middleSlice'); axis image; axis off; colormap(gca, colorMapName); colorbar; title(sprintf('%s (Z=%d)', titleText, middleSliceIdx), 'Interpreter', 'none'); drawnow; try close(figHandle); catch; end; catch; fprintf('WARN: Failed to visualize %s\n',titleText); end; end
% --- Path Setup Function ---
function addPaths(config); try; if isfield(config,'niftiToolsPath') && exist(config.niftiToolsPath, 'dir'); addpath(config.niftiToolsPath); fprintf('  Added NIfTI Path: %s\n', config.niftiToolsPath); end; catch; end; try; if isfield(config,'helperScriptsPath') && exist(config.helperScriptsPath, 'dir'); addpath(config.helperScriptsPath); fprintf('  Added Helpers Path: %s\n', config.helperScriptsPath); end; catch; end; end

% --- NIfTI functions (Ensure these are available via addPaths) ---
% Requires niftiinfo(), niftiread()
