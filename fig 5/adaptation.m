% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Simulate the model for several adaptive local-resource response rates.

% Generate a reproducible set of species-specific resource-use traits.
masterSeed = 1;
traitStream = RandStream('Threefry', 'Seed', masterSeed);
traitStream.Substream = 1;

% Consumer-resource model parameters.
model.S = 50;
model.k = 1;
model.gamma = 0.01;
model.Rstar = 5;
model.N0 = 1;
model.R0 = model.Rstar;
model.T0 = model.Rstar;
model.extinctionThreshold = 1e-3;
model.alpha = 1 + 4 .* rand(traitStream, 1, model.S);

% Fixed diffusion rate and the three adaptation rates to compare.
beta = 1;
etaValues = [1e-3, 0.5, 1.5];

% Simulation and output settings.
tMax = 3000;
outputInterval = 1;
explicitSafetyFactor = 0.5;

% Set a conservative base timestep cap. Each trajectory reduces this cap
% further at every step using its current resource and population rates.
maximumDt = min(0.001, explicitSafetyFactor ./ ...
    ((model.S + 1) .* beta));

% Store every trajectory at the same fixed physical times.
sampleTimes = (0:outputInterval:tMax).';
if sampleTimes(end) < tMax
    sampleTimes(end + 1, 1) = tMax;
end

nEta = numel(etaValues);
nSamples = numel(sampleTimes);
abundanceHistory = zeros(nSamples, model.S, nEta);
thetaHistory = zeros(nSamples, model.S, nEta);
finalAbundances = zeros(nEta, model.S);
minimumDtUsed = zeros(nEta, 1);
maximumDtUsed = zeros(nEta, 1);
numSteps = zeros(nEta, 1);

% Run every eta value from identical initial conditions.
for etaIndex = 1:nEta
    eta = etaValues(etaIndex);
    [currentAbundanceHistory, currentThetaHistory, currentFinalN, ...
        currentMinimumDt, currentMaximumDt, currentNumSteps] = ...
        simulateAdaptiveTrajectory( ...
        eta, beta, sampleTimes, maximumDt, explicitSafetyFactor, model);

    abundanceHistory(:, :, etaIndex) = currentAbundanceHistory;
    thetaHistory(:, :, etaIndex) = currentThetaHistory;
    finalAbundances(etaIndex, :) = currentFinalN;
    minimumDtUsed(etaIndex) = currentMinimumDt;
    maximumDtUsed(etaIndex) = currentMaximumDt;
    numSteps(etaIndex) = currentNumSteps;
end

% A species' theta dynamics are plotted only if that species survives to
% tMax in the corresponding eta trajectory.
survivingSpecies = finalAbundances > model.extinctionThreshold;

% Locate this source script so the plain Excel workbook is saved beside it,
% including when MATLAB runs a temporary Editor_* copy of the script.
candidateDirectories = {};
try
    activeFilename = matlab.desktop.editor.getActiveFilename;
    if ~isempty(activeFilename)
        candidateDirectories{end + 1} = fileparts(activeFilename);
    end
catch
    % MATLAB may be running without the desktop Editor.
end
scriptFullPath = mfilename('fullpath');
if ~isempty(scriptFullPath)
    candidateDirectories{end + 1} = fileparts(scriptFullPath);
end
resolvedScriptPath = which('adaptation');
if ~isempty(resolvedScriptPath)
    candidateDirectories{end + 1} = fileparts(resolvedScriptPath);
end
candidateDirectories{end + 1} = pwd;
candidateDirectories = unique(candidateDirectories, 'stable');

scriptDirectory = '';
for directoryIndex = 1:numel(candidateDirectories)
    candidateScriptPath = fullfile( ...
        candidateDirectories{directoryIndex}, 'adaptation.m');
    if isfile(candidateScriptPath)
        scriptDirectory = candidateDirectories{directoryIndex};
        break
    end
end
if isempty(scriptDirectory)
    error('Could not locate the folder containing adaptation.m.')
end

% Export one plain worksheet per eta. Abundances and theta values share the
% same time grid and are separated by one blank column for readability.
outputWorkbookPath = fullfile( ...
    scriptDirectory, 'adaptation_dynamics.xlsx');
abundanceHeaders = [{'time'}, arrayfun(@(speciesIndex) ...
    sprintf('abundance_%d', speciesIndex), 1:model.S, ...
    'UniformOutput', false)];
thetaHeaders = [{'time'}, arrayfun(@(speciesIndex) ...
    sprintf('theta_%d', speciesIndex), 1:model.S, ...
    'UniformOutput', false)];
thetaStartColumn = model.S + 3;
thetaHeaderRange = sprintf('%s1', ...
    excelColumnName(thetaStartColumn));
thetaDataRange = sprintf('%s2', ...
    excelColumnName(thetaStartColumn));

for etaIndex = 1:nEta
    sheetName = sprintf('eta=%.3g', etaValues(etaIndex));

    % Overwrite only this eta worksheet so rerunning the script refreshes
    % its data without adding spreadsheet formatting or Excel tables.
    writecell(abundanceHeaders, outputWorkbookPath, ...
        'Sheet', sheetName, 'Range', 'A1', ...
        'WriteMode', 'overwritesheet')
    writematrix([sampleTimes, abundanceHistory(:, :, etaIndex)], ...
        outputWorkbookPath, 'Sheet', sheetName, 'Range', 'A2')
    writecell(thetaHeaders, outputWorkbookPath, ...
        'Sheet', sheetName, 'Range', thetaHeaderRange)
    writematrix([sampleTimes, thetaHistory(:, :, etaIndex)], ...
        outputWorkbookPath, 'Sheet', sheetName, ...
        'Range', thetaDataRange)
end

% Apply the shared figure style.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
figureFontSize = 60;
curveLineWidth = 3;
axesLineWidth = 1.5;
figurePosition = [50 80 1800 850];
speciesColors = lines(model.S);

% Use a common abundance range so changes across eta panels are directly
% comparable. Display only reference ticks reached by the plotted range.
maximumAbundance = max(abundanceHistory, [], 'all');
abundanceUpperLimit = maximumAbundance;
candidateAbundanceTicks = [1e-2, 1e0, 1e2];
abundanceTicks = candidateAbundanceTicks( ...
    candidateAbundanceTicks <= abundanceUpperLimit);

% Give every eta case its own figure, with abundance and theta side by side.
for etaIndex = 1:nEta
    figureHandle = figure(etaIndex);
    clf(figureHandle)
    set(figureHandle, 'Color', 'w', 'Position', figurePosition, ...
        'Renderer', 'painters')
    figureLayout = tiledlayout(1, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    abundanceAxes = nexttile(figureLayout, 1);
    plotAbundancePanel(abundanceAxes, sampleTimes, ...
        abundanceHistory(:, :, etaIndex), etaValues(etaIndex), ...
        speciesColors, model.extinctionThreshold, abundanceUpperLimit, ...
        abundanceTicks, tMax, figureFontSize, curveLineWidth, ...
        axesLineWidth)

    thetaAxes = nexttile(figureLayout, 2);
    plotThetaPanel(thetaAxes, sampleTimes, ...
        thetaHistory(:, :, etaIndex), ...
        survivingSpecies(etaIndex, :), etaValues(etaIndex), ...
        speciesColors, tMax, figureFontSize, curveLineWidth, ...
        axesLineWidth)

    % Export the line plots as true vector objects inside the SVG file.
    svgOutputPath = fullfile(scriptDirectory, ...
        sprintf('adaptation_eta_%.3g.svg', etaValues(etaIndex)));
    print(figureHandle, svgOutputPath, '-dsvg', '-painters')
end

% Report integration details and the final survivor count for every eta.
for etaIndex = 1:nEta
    fprintf('eta = %.6g:\n', etaValues(etaIndex))
    fprintf('  Simulation completed to t = %.6g.\n', tMax)
    fprintf('  Adaptive Euler steps: %d.\n', numSteps(etaIndex))
    fprintf('  Timestep range used: [%.6g, %.6g].\n', ...
        minimumDtUsed(etaIndex), maximumDtUsed(etaIndex))
    fprintf('  Species above %.1e at tMax: %d of %d.\n', ...
        model.extinctionThreshold, ...
        nnz(survivingSpecies(etaIndex, :)), model.S)
end
fprintf('Plain trajectory data written to %s.\n', outputWorkbookPath)

function plotAbundancePanel(axesHandle, sampleTimes, abundanceValues, ...
    eta, speciesColors, extinctionThreshold, abundanceUpperLimit, ...
    abundanceTicks, tMax, figureFontSize, curveLineWidth, axesLineWidth)

% Plot every species abundance on a logarithmic y-axis in one square panel.
hold(axesHandle, 'on')
abundanceForPlot = abundanceValues;
abundanceForPlot(abundanceForPlot <= 0) = NaN;

for speciesIndex = 1:size(abundanceForPlot, 2)
    semilogy(axesHandle, sampleTimes, ...
        abundanceForPlot(:, speciesIndex), ...
        'Color', speciesColors(speciesIndex, :), ...
        'LineWidth', curveLineWidth)
end

xlabel(axesHandle, '$t$')
ylabel(axesHandle, '$N_i(t)$')
title(axesHandle, sprintf('$\\eta = %.3g$', eta))
xlim(axesHandle, [0 tMax])
xticks(axesHandle, [0, tMax ./ 2, tMax])
ylim(axesHandle, [extinctionThreshold, abundanceUpperLimit])
yticks(axesHandle, abundanceTicks)
set(axesHandle, 'YScale', 'log', 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth)
pbaspect(axesHandle, [1 1 1])
box(axesHandle, 'on')
end

function plotThetaPanel(axesHandle, sampleTimes, thetaValues, ...
    survivingSpecies, eta, speciesColors, tMax, figureFontSize, ...
    curveLineWidth, axesLineWidth)

% Plot theta only for species that survive to the final simulation time.
hold(axesHandle, 'on')
survivingIndices = find(survivingSpecies);
for speciesIndex = survivingIndices
    plot(axesHandle, sampleTimes, thetaValues(:, speciesIndex), ...
        'Color', speciesColors(speciesIndex, :), ...
        'LineWidth', curveLineWidth)
end

xlabel(axesHandle, '$t$')
ylabel(axesHandle, '$\theta_i(t)$')
title(axesHandle, sprintf('$\\eta = %.3g$', eta))
xlim(axesHandle, [0 tMax])
xticks(axesHandle, [0, tMax ./ 2, tMax])
ylim(axesHandle, [0 1])
set(axesHandle, 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth)
pbaspect(axesHandle, [1 1 1])
box(axesHandle, 'on')
end

function columnName = excelColumnName(columnNumber)

% Convert a positive Excel column number to its A1-style letter label.
if ~isscalar(columnNumber) || columnNumber < 1 || ...
        columnNumber ~= floor(columnNumber)
    error('Excel column numbers must be positive integers.')
end

columnName = '';
while columnNumber > 0
    remainder = mod(columnNumber - 1, 26);
    columnName = [char(65 + remainder), columnName]; %#ok<AGROW>
    columnNumber = floor((columnNumber - 1) ./ 26);
end
end

function [abundanceHistory, thetaHistory, N, minimumDtUsed, ...
    maximumDtUsed, numSteps] = simulateAdaptiveTrajectory( ...
    eta, beta, sampleTimes, maximumDt, explicitSafetyFactor, model)

% Simulate one deterministic trajectory from the common initial state.
N = model.N0 .* ones(1, model.S);
R = model.R0;
T = model.T0 .* ones(1, model.S);
theta = zeros(1, model.S);

nSamples = numel(sampleTimes);
tMax = sampleTimes(end);
abundanceHistory = zeros(nSamples, model.S);
thetaHistory = zeros(nSamples, model.S);
abundanceHistory(1, :) = N;
thetaHistory(1, :) = theta;

sampleIndex = 2;
currentTime = 0;
timeTolerance = 10 .* eps(max(1, tMax));
minimumDtUsed = inf;
maximumDtUsed = 0;
numSteps = 0;

while currentTime < tMax
    % Calculate the current per-capita growth rate for every species.
    resourceUse = theta .* T + (1 - theta) .* R;
    perCapitaGrowth = model.alpha .* resourceUse - model.k;

    % Calculate the deterministic consumer-resource rates.
    dN = perCapitaGrowth .* N;
    dT = beta .* (R - T) - ...
        model.alpha .* model.gamma .* theta .* T .* N;
    dR = model.Rstar - R - beta .* ...
        (model.S .* R - sum(T)) - model.gamma .* R .* ...
        sum(model.alpha .* (1 - theta) .* N);

    % Bound the explicit step using resource removal, demographic change,
    % and theta adaptation. Extinct species do not contribute to the two
    % species-level rate bounds because their N_i and theta_i are frozen.
    survivingNow = N > model.extinctionThreshold;
    localRemovalRate = beta + ...
        model.alpha .* model.gamma .* theta .* N;
    sharedRemovalRate = 1 + model.S .* beta + ...
        model.gamma .* sum(model.alpha .* (1 - theta) .* N);

    if any(survivingNow)
        demographicRate = max(abs(perCapitaGrowth(survivingNow)));
        adaptationRate = max(eta .* ...
            abs(perCapitaGrowth(survivingNow)));
    else
        demographicRate = 0;
        adaptationRate = 0;
    end

    maximumCurrentRate = max([localRemovalRate, sharedRemovalRate, ...
        demographicRate, adaptationRate]);

    if ~isfinite(maximumCurrentRate) || maximumCurrentRate <= 0
        error(['An invalid adaptive rate occurred at t = %.6g ' ...
            'for eta = %.6g.'], currentTime, eta)
    end

    dt = min(maximumDt, ...
        explicitSafetyFactor ./ maximumCurrentRate);
    dt = min(dt, tMax - currentTime);

    % Land exactly on each requested output time.
    if sampleIndex <= nSamples
        dt = min(dt, sampleTimes(sampleIndex) - currentTime);
    end

    if ~isfinite(dt) || dt <= 0
        error(['An invalid timestep occurred at t = %.6g ' ...
            'for eta = %.6g.'], currentTime, eta)
    end

    % Advance abundances and resources explicitly.
    NNew = N + dt .* dN;
    TNew = T + dt .* dT;
    RNew = R + dt .* dR;

    % Force extinct populations to zero immediately and permanently.
    NNew(NNew < model.extinctionThreshold) = 0;

    % Advance theta semi-implicitly with growth held fixed during this step:
    % dtheta_i/dt = eta[(1-theta_i)max(-g_i,0)
    %                  - theta_i max(g_i,0)].
    switchLocal = eta .* max(-perCapitaGrowth, 0);
    switchGlobal = eta .* max(perCapitaGrowth, 0);
    thetaNew = theta;

    negativeGrowth = perCapitaGrowth < 0 & NNew > 0;
    positiveGrowth = perCapitaGrowth > 0 & NNew > 0;

    thetaNew(negativeGrowth) = ...
        (theta(negativeGrowth) + dt .* switchLocal(negativeGrowth)) ...
        ./ (1 + dt .* switchLocal(negativeGrowth));
    thetaNew(positiveGrowth) = theta(positiveGrowth) ...
        ./ (1 + dt .* switchGlobal(positiveGrowth));

    % Extinct populations retain their last theta value, and roundoff is
    % clipped so every numerical theta remains within its biological bounds.
    thetaNew(NNew == 0) = theta(NNew == 0);
    thetaNew = min(max(thetaNew, 0), 1);

    % Commit the completed step.
    N = NNew;
    T = TNew;
    R = RNew;
    theta = thetaNew;
    currentTime = currentTime + dt;
    numSteps = numSteps + 1;
    minimumDtUsed = min(minimumDtUsed, dt);
    maximumDtUsed = max(maximumDtUsed, dt);

    % Stop rather than silently retain an invalid trajectory.
    if any(~isfinite(N), 'all') || ~isfinite(R) || ...
            any(~isfinite(T), 'all') || any(~isfinite(theta), 'all')
        error(['A nonfinite state occurred at t = %.6g ' ...
            'for eta = %.6g.'], currentTime, eta)
    end
    if R < 0 || any(T < 0, 'all')
        error(['A negative resource state occurred at t = %.6g ' ...
            'for eta = %.6g.'], currentTime, eta)
    end

    % Record a requested output time once it is reached to floating-point
    % precision. The step restriction above prevents skipped samples.
    if sampleIndex <= nSamples && ...
            currentTime >= sampleTimes(sampleIndex) - timeTolerance
        abundanceHistory(sampleIndex, :) = N;
        thetaHistory(sampleIndex, :) = theta;
        sampleIndex = sampleIndex + 1;
    end
end
end
