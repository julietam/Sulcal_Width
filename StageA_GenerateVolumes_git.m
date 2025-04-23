% =========================================================================
% Stage A: Generate Intermediate Volumes (Single Subject Batch Version) V4.1
% =========================================================================
% Description:
% Runs Stage 1 processing and label generation steps for a single subject.
% Uses the 'cd' workaround within callRegionGrowing24 helper function
% to address the tool's output naming issue with full paths.
% Uses consistent permutation: Permute [2 1 3] before saving raw files.
% =========================================================================

function StageA_GenerateVolumes() % Wrapped in a function

    fprintf('--- Starting StageA_GenerateVolumes ---\n');
    startTime = tic; % Start timer

    % --- Configuration ---
    config = struct();
    % Paths (VERIFY THESE)
    config.rootFolder = fullfile('/Users/yourname/Documents/NeuroData/FS_Outputs');
    config.niftiToolsPath = '/Users/yourname/Scripts_matlab/Plugins_ML/NIfTI_20140122';
    config.bwdistscPath = '/Users/yourname/Scripts_matlab/Plugins_ML/bwdistsc';
    config.helperScriptsPath = '/Users/yourname/Scripts_matlab';
    % *** VERIFY this executable path and name is correct ***
    config.rg24Path = '/Users/yourname//RG2020/RGmain/RegionGrowing24'; % Path to LOCAL MAXIMA tool
    config.mriConvertPath = 'mri_convert'; % Command/path for FreeSurfer's mri_convert

    % File Patterns and Names
    config.brainNiiFilename = 'brain.nii'; % Relative to mri subfolder (assumed)
    config.ribbonNiiFilename = 'ribbon.nii'; % Relative to mri subfolder (needed for defineDatasetPaths)
    config.labelPattern = '*.label.nii'; % Pattern for original label files

    % Processing Parameters
    config.imageSize = [256, 256, 256]; % Expected dimensions [X Y Z]
    config.closingRadius = 3; % Radius for morphological closing
    config.rg24Directions = 7; % Directions for RG24 local max (Verify)
    config.sulciDilateRadius = 2; % Radius for final label dilation

    % Execution Control
    config.verbose = true; % Print progress messages in logs
    config.showPlots = false; % SET TO FALSE FOR BATCH JOBS
    config.overwriteOutputs = true; % Typically true when rerunning batch jobs
    config.logFolder = fullfile(config.rootFolder, 'logs'); % Define log folder
    config.errorLogFile = fullfile(config.logFolder, 'error_log.txt');
    % --- End Configuration ---

    % --- Setup ---
    if ~isdeployed % Only add paths if running in standard MATLAB
        addPaths(config);
    end
    if ~exist(config.logFolder, 'dir'); try mkdir(config.logFolder); catch; end; end

    % --- Get Subject ID from Environment Variable ---
    datasetID = getenv('JOB_SUBJECT_ID');
    if isempty(datasetID); error('Environment variable JOB_SUBJECT_ID not set.'); end
    fprintf('\n--- Processing Dataset: %s ---\n', datasetID);

    % --- Setup for SINGLE Subject ---
    datasetFolder = fullfile(resolveTilde(config.rootFolder), datasetID);
    if ~exist(datasetFolder, 'dir'); error('Dataset folder not found: %s', datasetFolder); end
    paths = defineDatasetPaths(datasetFolder, config.brainNiiFilename, config.ribbonNiiFilename, datasetID);
    imSize = config.imageSize; % Use configured size initially

    try % Master try-catch for the whole dataset processing

        fprintf('--- Running Stage A Volume Generation Steps ---\n');

        % --- Step 1: Load Brain & Binarize ---
        fprintf('  Loading brain volume: %s\n', paths.brainNii);
        if ~exist(paths.brainNii, 'file'); error('Input brain NIfTI not found: %s', paths.brainNii); end
        voxelMat = niftiread(paths.brainNii);
        voxelInfo = niftiinfo(paths.brainNii);
        if ~isequal(voxelInfo.ImageSize(1:3), imSize); warning('Dataset %s: NIfTI Image size mismatch. Adapting.', datasetID); imSize = voxelInfo.ImageSize(1:3); fprintf('  Adapted imSize to [%d %d %d]\n', imSize(1), imSize(2), imSize(3)); end
        voxelSize = voxelInfo.PixelDimensions(1:3);
        fprintf('  Binarizing image...\n');
        voxelMatBinary = imbinarize(voxelMat);
        if nnz(voxelMatBinary) == 0; error('Input brain mask %s is empty!', paths.brainNii); end

        % --- Step 2: Distance Map ---
        fprintf('  Computing anisotropic distance map...\n');
        distMap = computeDistanceMap(voxelMatBinary, voxelSize);
        fprintf('  Saving distance map (Raw): %s\n', paths.distMapRaw);
        saveVolume(distMap, paths.distMapRaw, 'float'); % Permutes before save
        fprintf('  Saving distance map (NIfTI): %s\n', paths.distMapNii);
        saveNiftiWithInfo(distMap, paths.distMapNii, voxelInfo);

        % --- Step 3: Brain Closure ---
        fprintf('  Generating morphological closure...\n');
        brainClosure = generateBrainClosure(voxelMatBinary, config.closingRadius);
        brainClosure = double(brainClosure);
        if nnz(brainClosure) == 0; warning('Dataset %s: Brain closure result is all zeros!', datasetID); else fprintf('      Brain closure NNZ: %d\n', nnz(brainClosure)); end
        fprintf('  Saving brain closure (Raw): %s\n', paths.brainClosureRaw);
        saveVolume(brainClosure, paths.brainClosureRaw, 'double'); % Permutes before save
        fprintf('  Saving brain closure (NIfTI): %s\n', paths.brainClosureNii);
        saveNiftiWithInfo(brainClosure, paths.brainClosureNii, voxelInfo);

        % --- Step 4: Masked Distance Map (DistMap * BrainClosure) ---
        fprintf('  Calculating Masked Distance Map (DistMap * BrainClosure)...\n');
        maskedDistMap = distMap .* brainClosure;
        fprintf('      Masked DistMap NNZ=%d\n', nnz(maskedDistMap));
        paths.maskedDistMapRaw = fullfile(paths.datasetFolder, 'Brain_masked_distmap_anisotropic.raw');
        fprintf('  Saving Masked Distance Map (Raw): %s\n', paths.maskedDistMapRaw);
        saveVolume(maskedDistMap, paths.maskedDistMapRaw, 'float'); % Permutes before save

        % --- Step 5: Local Maxima using RG24 ---
        fprintf('  Calling RegionGrowing24 on MASKED distance map...\n');
        callRegionGrowing24(paths.maskedDistMapRaw, paths.localMaxRaw, imSize, config); % Use the 'cd' workaround version defined below
        if ~exist(paths.localMaxRaw, 'file'); error('RG24 failed to create expected output file: %s', paths.localMaxRaw); end; % Final check still useful
        fprintf('  Local maxima file step completed: %s\n', paths.localMaxRaw);
        % Load to check type and save Nii - use loadRawImage (permutes after load)
        localMaxVolRaw = loadRawImage(paths.localMaxRaw, imSize, 'uint8', 'ieee-le');
        fprintf('  Saving local maxima map (NIfTI): %s\n', paths.localMaxNii);
        saveNiftiWithInfo(localMaxVolRaw, paths.localMaxNii, voxelInfo);

        % --- Step 6: Combine Sulci Labels ---
        fprintf('  Combining sulci labels from: %s\n', paths.labelsFolder);
        originalLabelPattern = config.labelPattern;
        if ~isempty(dir(fullfile(paths.labelsFolder, [originalLabelPattern '.gz'])))
             fprintf('    Unzipping labels...\n');
             gunzip(fullfile(paths.labelsFolder, [originalLabelPattern '.gz']), paths.labelsFolder);
        end
        combinedSulci = combineSulciLabels(paths.labelsFolder, originalLabelPattern, paths.brainNii, config.mriConvertPath);
        if isempty(combinedSulci) || all(combinedSulci(:)==0); warning('Dataset %s: Combined sulci map is empty/zeros.', datasetID); combinedSulci = zeros(imSize, 'int32'); else fprintf('      Combined map NNZ=%d, Max=%d\n', nnz(combinedSulci), max(combinedSulci(:))); end
        fprintf('  Saving initial combined sulci map (Raw - INT32): %s\n', paths.combinedSulciRaw);
        saveraw_int(combinedSulci, paths.combinedSulciRaw, 'int32'); % Permutes before save
        fprintf('  Saving initial combined sulci map (NIfTI): %s\n', paths.combinedSulciNii);
        saveNiftiWithInfo(int32(combinedSulci), paths.combinedSulciNii, voxelInfo);

        % --- Step 7: Dilate Combined Labels ---
        fprintf('  Dilating INITIAL combined sulci labels (radius %d)...\n', config.sulciDilateRadius);
        SL2 = dilateSulciLabels(combinedSulci, config.sulciDilateRadius);
        fprintf('  Saving dilated initial labels (Raw - INT32): %s\n', paths.sulciDilatedRaw);
        saveVolume(SL2, paths.sulciDilatedRaw, 'int32'); % Permutes before save
        fprintf('  Saving dilated initial labels (NIfTI): %s\n', paths.sulciDilatedNii);
        saveNiftiWithInfo(int32(SL2), paths.sulciDilatedNii, voxelInfo);

        totalTime = toc(startTime);
        fprintf('--- Stage A Volume Generation for %s Completed Successfully (%.2f seconds) ---\n', datasetID, totalTime);

    catch ME % Catch errors for this dataset
         totalTime = toc(startTime);
         fprintf(2,'\n!!! ERROR in Stage A for dataset %s after %.2f sec: %s !!!\n', datasetID, totalTime, ME.message); % Use stderr
         logError(ME, config.errorLogFile, datasetID, 'Stage A');
         exit(1); % Exit MATLAB with non-zero status on error
    end

    exit(0); % Explicitly exit MATLAB with success status

end % End of main function wrapper StageA_GenerateVolumes

% =========================================================================
% Helper Function Definitions for Stage A
% =========================================================================

% --- I/O Functions (Permute Before Save, Permute After Load) ---
function saveVolume(volume, fileID, dataType); fprintf('      Permuting [2 1 3] before saving raw file: %s\n', fileID); volumeToSave = permute(volume, [2, 1, 3]); try fid = fopen(fileID, 'w', 'ieee-le'); if fid == -1; error('Cannot open %s', fileID); end; fwrite(fid, volumeToSave, dataType); fclose(fid); catch ME; if fid ~= -1; fclose(fid); end; error('Failed write raw %s: %s', fileID, ME.message); end; end
function saveraw_int(volume, fileID, dataType); fprintf('      Permuting [2 1 3] before saving raw integer file: %s\n', fileID); volumeToSave = permute(volume, [2, 1, 3]); try fid = fopen(fileID, 'w', 'ieee-le'); if fid == -1; error('Cannot open %s', fileID); end; fwrite(fid, volumeToSave, dataType); fclose(fid); catch ME; if fid ~= -1; fclose(fid); end; error('Failed write int raw %s: %s', fileID, ME.message); end; end
function volume = loadRawImage(rawFile, imSize, dataType, byteOrder); if ~exist(rawFile,'file'); error('Raw file not found: %s', rawFile); end; fprintf('      Loading raw file: %s (Type: %s)\n', rawFile, dataType); try volume_temp = multibandread(rawFile, [imSize(1), imSize(2), imSize(3)], dataType, 0, 'bsq', byteOrder); fprintf('      Applying permute([2 1 3]) after loading raw file.\n'); volume = permute(volume_temp, [2, 1, 3]); clear volume_temp; catch ME; error('Failed load/permute raw %s: %s', rawFile, ME.message); end; end

% --- Setup Functions ---
function addPaths(config); try; if isfield(config,'niftiToolsPath') && exist(config.niftiToolsPath, 'dir'); addpath(config.niftiToolsPath); fprintf('  Added NIfTI Path: %s\n', config.niftiToolsPath); end; catch; end; try; if isfield(config,'bwdistscPath') && exist(config.bwdistscPath, 'dir'); addpath(config.bwdistscPath); fprintf('  Added bwdistsc Path: %s\n', config.bwdistscPath); end; catch; end; try; if isfield(config,'helperScriptsPath') && exist(config.helperScriptsPath, 'dir'); addpath(config.helperScriptsPath); fprintf('  Added Helpers Path: %s\n', config.helperScriptsPath); end; catch; end; end % Made more robust
function logError(ME, logFilePath, datasetID, stage); resolvedLogPath = resolveTilde(logFilePath); try fid = fopen(resolvedLogPath, 'a'); if fid == -1; fprintf(2,'!!! Cannot open log file: %s !!!\n', resolvedLogPath); return; end; timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS'); fprintf(fid, 'Timestamp: %s\nDataset: %s\nStage: %s\n', timestamp, datasetID, stage); fprintf(fid, 'Error ID: %s\nError Msg: %s\n', ME.identifier, ME.message); if ~isempty(ME.stack); fprintf(fid, 'Stack Trace:\n'); for k = 1:length(ME.stack); fprintf(fid, '  File: %s, Name: %s, Line: %d\n', ME.stack(k).file, ME.stack(k).name, ME.stack(k).line); end; end; fprintf(fid, '----------------------------------------\n'); fclose(fid); catch logEx; fprintf(2,'!!! FAILED TO WRITE TO ERROR LOG FILE (%s): %s !!!\n', resolvedLogPath, logEx.message); end; end
function resolvedPath = resolveTilde(filePath); if ispc || isempty(filePath) || ~ischar(filePath); resolvedPath = filePath; return; end; if startsWith(filePath, '~'); try; [status, homeDir] = system('echo $HOME'); if status == 0 && ~isempty(strtrim(homeDir)); homeDir = strtrim(homeDir); if strcmp(filePath, '~'); resolvedPath = homeDir; elseif startsWith(filePath, '~/'); resolvedPath = fullfile(homeDir, filePath(3:end)); else; resolvedPath = fullfile(homeDir, filePath(2:end)); end; else; resolvedPath = filePath; end; catch; resolvedPath = filePath; end; else; resolvedPath = filePath; end; end
function base = basename(filepath); [~, name, ext] = fileparts(filepath); base = [name ext]; end
function paths = defineDatasetPaths(datasetFolder, brainNiiFilename, ribbonNiiFilename, datasetID)
    paths = struct();
    paths.datasetFolder = resolveTilde(datasetFolder); paths.mriFolder = fullfile(paths.datasetFolder, 'mri'); paths.labelsFolder = fullfile(paths.datasetFolder, 'labels');
    paths.brainNii = fullfile(paths.datasetFolder, brainNiiFilename); % Corrected path to mriFolder
    paths.distMapRaw = fullfile(paths.datasetFolder, 'Brain_distmap_anisotropic.raw');
    paths.maskedDistMapRaw = fullfile(paths.datasetFolder, 'Brain_masked_distmap_anisotropic.raw');
    paths.localMaxRaw = fullfile(paths.datasetFolder, 'localmax_Brain_distmap_anisotropic.raw');
    paths.distMapNii = fullfile(paths.datasetFolder, 'Brain_distmap_anisotropic.nii');
    paths.brainClosureRaw = fullfile(paths.datasetFolder, 'brain_closing.raw');
    paths.brainClosureNii = fullfile(paths.datasetFolder, 'brain_closing.nii');
    paths.localMaxNii = fullfile(paths.datasetFolder, 'localmax_Brain_distmap_anisotropic.nii');
    paths.combinedSulciRaw = fullfile(paths.datasetFolder, 'Brain_sulci.raw');
    paths.combinedSulciNii = fullfile(paths.datasetFolder, 'Brain_sulci.nii');
    paths.sulciDilatedRaw = fullfile(paths.datasetFolder, 'sulci_dil.raw');
    paths.sulciDilatedNii = fullfile(paths.datasetFolder, 'sulci_dil.nii');
    paths.ribbonNii = fullfile(paths.mriFolder, ribbonNiiFilename);
end

% --- Stage 1 Helpers ---
function distMap = computeDistanceMap(voxelMatBinary, voxelSize); if ~exist('bwdistsc', 'file'); error('bwdistsc function not found. Ensure path is added in config.'); end; distMap = bwdistsc(voxelMatBinary, voxelSize); end
function brainClosure = generateBrainClosure(voxelMatBinary, radius); SE = strel('sphere', radius); brainClosure = imclose(voxelMatBinary, SE); end

% --- callRegionGrowing24 using 'cd' WORKAROUND ---
function callRegionGrowing24(inputFile, outputFile, imSize, config)
    resolvedRg24Path = resolveTilde(config.rg24Path); % Path to executable
    if ~exist(resolvedRg24Path, 'file'); error('RG24 executable not found: %s', config.rg24Path); end

    inputFileFullPath = resolveTilde(inputFile); % Full path to input
    outputFileFullPath = resolveTilde(outputFile); % Full path to desired final output

    if ~exist(inputFileFullPath, 'file'); error('Input file for RG24 not found: %s', inputFileFullPath); end

    [dataDir, inputBaseName, inputExt] = fileparts(inputFileFullPath); % Get directory and input filename
    inputRelative = [inputBaseName, inputExt]; % Relative input filename
    % Assume tool's default output name is based on relative input name
    outputBaseName = ['localmax_' inputBaseName '.raw']; %<<<<<<<<< This is the name RG24 TRIES to create by default
    expectedOutputFullPathInDir = fullfile(dataDir, outputBaseName); % Full path to where it likely tries to write

    currentDir = pwd; % Store current directory
    cleanupObj = onCleanup(@() cd(currentDir)); % Ensure we change back even on error

    try
        fprintf('  Changing directory to: %s\n', dataDir);
        cd(dataDir); % Change to the data directory

        % Delete existing expected output file in target directory BEFORE running
         if exist(outputBaseName, 'file')
             fprintf('  Deleting existing output file before RG24 run: %s\n', outputBaseName);
             delete(outputBaseName);
         end
         % Also delete the final target file just in case rename fails later
         if exist(outputFileFullPath, 'file') && ~strcmp(expectedOutputFullPathInDir, outputFileFullPath)
             fprintf('  Deleting existing final target file before RG24 run: %s\n', outputFileFullPath);
             delete(outputFileFullPath);
         end


        dimsStr = sprintf('%d %d %d', imSize(1), imSize(2), imSize(3));
        % Construct command using ONLY RELATIVE input filename
        cmd = sprintf('"%s" "%s" -local_maxF -dims %s -num_directions %d', ...
                      resolvedRg24Path, inputRelative, dimsStr, config.rg24Directions);
        % NOTE: No -o flag used here!

        fprintf('  Executing RG24 (in data dir): %s\n', cmd);
        [status, cmdout] = system(cmd);
        if config.verbose || status ~= 0; fprintf('  RG24 Output:\n%s\n', cmdout); end

        cd(currentDir); % Change back immediately
        fprintf('  Restored directory to: %s\n', currentDir);

        if status ~= 0
             error('RegionGrowing24 failed with status %d.', status);
        end

        % Verify the default output file name was created IN dataDir
        if ~exist(expectedOutputFullPathInDir, 'file')
             error('Expected default RG24 output file %s was not found after execution in %s.', outputBaseName, dataDir);
        else
             fprintf('  Verified default RG24 output file exists: %s\n', expectedOutputFullPathInDir);
             % Rename IF the automatically generated name differs from the desired final name from paths struct
             if ~strcmp(expectedOutputFullPathInDir, outputFileFullPath)
                 fprintf('  Renaming %s to %s\n', expectedOutputFullPathInDir, outputFileFullPath);
                 try
                    movefile(expectedOutputFullPathInDir, outputFileFullPath);
                 catch ME_move
                     error('Failed to move RG24 output file %s to %s: %s', expectedOutputFullPathInDir, outputFileFullPath, ME_move.message);
                 end
             end
        end

        % Final check for the desired file name defined in paths struct
         if ~exist(outputFileFullPath, 'file')
             error('Final RG24 output file %s not found after execution/renaming.', outputFileFullPath);
         else
              fprintf('  Confirmed final RG24 output file %s exists.\n', outputFileFullPath);
         end

    catch ME
        cd(currentDir); % Ensure change back on error
        rethrow(ME);
    end
end % End of function callRegionGrowing24 (cd version)


% --- Stage 2 Helpers (needed by Stage A) ---
function combinedSulci = combineSulciLabels(sulciFolder, labelPattern, referenceNii, mriConvertPath)
    fprintf('    Running combineSulciLabels (Internal Reslicing)...\n'); combinedSulci = []; filesAddedCounter = 0;
    try
        resolvedReferenceNii = resolveTilde(referenceNii); if ~exist(resolvedReferenceNii,'file'); error('Reference Nii missing: %s', resolvedReferenceNii); end; refInfo = niftiinfo(resolvedReferenceNii); imSize = refInfo.ImageSize; labels = dir(fullfile(sulciFolder, labelPattern)); if isempty(labels); warning('No labels matching "%s" found in %s.', labelPattern, sulciFolder); combinedSulci = zeros(imSize, 'int32'); return; end; % Return zeros if no labels
        fprintf('    Found %d files matching "%s".\n', length(labels), labelPattern); combinedSulci = zeros(imSize, 'int32'); fprintf('    Initialized combinedSulci size: %s\n', mat2str(imSize)); mriConvertCmd = mriConvertPath; tempResliceDir = fullfile(sulciFolder, 'temp_resliced_labels'); if ~exist(tempResliceDir, 'dir'); mkdir(tempResliceDir); end; cleanupObj = onCleanup(@() rmdir(tempResliceDir, 's'));
        for k = 1:length(labels)
            baseFileName = labels(k).name; fullFileName = fullfile(labels(k).folder, baseFileName); fprintf('    Label %d/%d: %s\n', k, length(labels), baseFileName);
            try; originalLabelData = niftiread(fullFileName); if nnz(originalLabelData) == 0; continue; end; catch; fprintf(' WARN: Cannot read %s\n',baseFileName); continue; end;
            [~, name_noext, ~] = fileparts(baseFileName); fullFileName_res = fullfile(tempResliceDir, sprintf('%s_res.nii', name_noext));
            cmd_res = sprintf('%s --reslice_like "%s" "%s" "%s" -rt nearest', mriConvertCmd, resolvedReferenceNii, fullFileName, fullFileName_res);
            [status, ~] = system(cmd_res); % Suppress system output unless error?
            if status ~= 0 || ~exist(fullFileName_res, 'file'); fprintf(' WARN: mri_convert failed for %s\n', baseFileName); continue; end;
            try
                niftiResData = niftiread(fullFileName_res);
                if ~isequal(size(niftiResData), refInfo.ImageSize); fprintf(' WARN: Resliced size mismatch %s\n', baseFileName); continue; end;
                if nnz(niftiResData) > 0; combinedSulci(niftiResData > 0) = int32(k); filesAddedCounter = filesAddedCounter + 1; end % Overwrite logic
            catch ME_load
                fprintf(' WARN: Failed load/process resliced %s: %s\n', fullFileName_res, ME_load.message); continue;
            end
        end
    catch ME_combine
        error('Error in combineSulciLabels: %s', ME_combine.message);
    end
    fprintf('    Combined labels from %d files. NNZ=%d\n', filesAddedCounter, nnz(combinedSulci)); if filesAddedCounter == 0; warning('combineSulciLabels: Final map empty.'); end
end

function dilatedSulci = dilateSulciLabels(sulciLabels, radius); if isempty(sulciLabels); dilatedSulci = []; return; end; if radius > 0; SE = strel('sphere', radius); dilatedSulci = imdilate(sulciLabels, SE); else dilatedSulci = sulciLabels; end; end
function saveNiftiWithInfo(volume, niftiFile, refNiiInfo); try; matlabType = class(volume); niftiType = ''; switch matlabType; case 'uint8'; niftiType = 'uint8'; bitpix = 8; case 'int16'; niftiType = 'int16'; bitpix = 16; case 'int32'; niftiType = 'int32'; bitpix = 32; case 'single'; niftiType = 'single'; bitpix = 32; case 'double'; niftiType = 'double'; bitpix = 64; case 'logical'; volume = uint8(volume); niftiType = 'uint8'; bitpix = 8; otherwise; error('Unsupported NIfTI type: %s', matlabType); end; newInfo = refNiiInfo; newInfo.Filename = niftiFile; newInfo.ImageSize = size(volume); try; newInfo.PixelDimensions = refNiiInfo.PixelDimensions(1:length(newInfo.ImageSize)); catch; newInfo.PixelDimensions = [1 1 1]; end; newInfo.Datatype = niftiType; newInfo.BitsPerPixel = bitpix; newInfo.Description = sprintf('Generated by MATLAB %s', datestr(now)); niftiwrite(volume, niftiFile, newInfo, 'Compressed', true); catch ME; error('Failed write NIfTI %s: %s', niftiFile, ME.message); end; end
function visualizeMiddleSlice(volume, titleText, colorMapName); if nargin < 3; colorMapName = 'gray'; end; try; if isempty(volume) || ndims(volume) ~= 3; return; end; middleSliceIdx = round(size(volume, 3) / 2); if islogical(volume); volume = uint8(volume); elseif ~isfloat(volume); volume = single(volume); end; middleSlice = volume(:, :, middleSliceIdx); figHandle=figure('Name', titleText, 'NumberTitle', 'off','Visible','off'); imagesc(middleSlice'); axis image; axis off; colormap(gca, colorMapName); colorbar; title(sprintf('%s (Z=%d)', titleText, middleSliceIdx), 'Interpreter', 'none'); drawnow; try close(figHandle); catch; end; catch; fprintf('WARN: Failed to visualize %s\n',titleText); end; end % Make invisible/close

% --- NIfTI functions (Ensure these are available or paths added) ---
% function info = niftiinfo(filename); ... end
% function data = niftiread(filename_or_info); ... end
% function niftiwrite(vol, filename, info, varargin); ... end
