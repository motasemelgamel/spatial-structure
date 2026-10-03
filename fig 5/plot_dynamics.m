% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Reproduce the saved abundance and theta dynamics without rerunning the
% stochastic consumer-resource simulation.
inputFileName = 'adaptation_dynamics.xlsx';

% Locate the workbook beside this script, including when MATLAB runs a
% temporary Editor_* copy from another working directory.
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
resolvedScriptPath = which('plot_dynamics');
if ~isempty(resolvedScriptPath)
    candidateDirectories{end + 1} = fileparts(resolvedScriptPath);
end
candidateDirectories{end + 1} = pwd;
candidateDirectories = unique(candidateDirectories, 'stable');

scriptDirectory = '';
for directoryIndex = 1:numel(candidateDirectories)
    candidateWorkbookPath = fullfile( ...
        candidateDirectories{directoryIndex}, inputFileName);
    candidateScriptPath = fullfile( ...
        candidateDirectories{directoryIndex}, 'plot_dynamics.m');
    if isfile(candidateWorkbookPath) && isfile(candidateScriptPath)
        scriptDirectory = candidateDirectories{directoryIndex};
        break
    end
end
if isempty(scriptDirectory)
    error('Could not locate %s beside plot_dynamics.m.', inputFileName)
end
inputWorkbookPath = fullfile(scriptDirectory, inputFileName);

% Read the two plain data sheets and identify their species columns.
abundanceTable = readtable(inputWorkbookPath, 'Sheet', 'Abundances', ...
    'VariableNamingRule', 'preserve');
thetaTable = readtable(inputWorkbookPath, 'Sheet', 'Theta', ...
    'VariableNamingRule', 'preserve');

abundanceNames = string(abundanceTable.Properties.VariableNames);
thetaNames = string(thetaTable.Properties.VariableNames);
abundanceColumns = find(startsWith(abundanceNames, 'species_'));
thetaColumns = find(startsWith(thetaNames, 'theta_'));

if isempty(abundanceColumns) || ...
        numel(abundanceColumns) ~= numel(thetaColumns)
    error('The workbook must contain matching species_* and theta_* columns.')
end

timeValues = abundanceTable{:, 'time'};
thetaTime = thetaTable{:, 'time'};
abundanceHistory = abundanceTable{:, abundanceColumns};
thetaHistory = thetaTable{:, thetaColumns};

% Read the saved parameter values used in the panel titles and extinction
% filter. They are repeated down their columns to keep the workbook flat.
etaValues = abundanceTable{:, 'eta'};
sigmaValues = abundanceTable{:, 'sigma'};
extinctionThresholdValues = abundanceTable{:, 'extinction_threshold'};
exampleEta = etaValues(1);
exampleSigma = sigmaValues(1);
extinctionThreshold = extinctionThresholdValues(1);

% Validate the saved arrays before plotting them.
timeValues = timeValues(:);
thetaTime = thetaTime(:);
nSpecies = numel(abundanceColumns);
timeTolerance = 100 .* eps(max(1, max(abs(timeValues), [], 'omitnan')));
if numel(timeValues) ~= numel(thetaTime) || ...
        any(abs(timeValues - thetaTime) > timeTolerance) || ...
        any(diff(timeValues) < 0) || ...
        ~isequal(size(abundanceHistory), size(thetaHistory)) || ...
        any(~isfinite(timeValues)) || ...
        any(~isfinite(abundanceHistory), 'all') || ...
        any(~isfinite(thetaHistory), 'all') || ...
        any(abundanceHistory < 0, 'all') || ...
        any(thetaHistory < 0 | thetaHistory > 1, 'all')
    error('The saved abundance or theta trajectories are invalid.')
end
if any(abs(etaValues - exampleEta) > timeTolerance) || ...
        any(abs(sigmaValues - exampleSigma) > timeTolerance) || ...
        any(abs(extinctionThresholdValues - extinctionThreshold) > ...
        timeTolerance) || ~isfinite(exampleEta) || ...
        ~isfinite(exampleSigma) || ~isfinite(extinctionThreshold) || ...
        extinctionThreshold <= 0
    error('The saved dynamics metadata are invalid or inconsistent.')
end

% Apply the same plotting format and species colors as the original script.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
figureFontSize = 60;
curveLineWidth = 3;
axesLineWidth = 1.5;
speciesColors = lines(nSpecies);

tMinimum = min(timeValues);
tMaximum = max(timeValues);
panelTitle = sprintf('$\\eta = %.3g,\\; \\sigma = %.3g$', ...
    exampleEta, exampleSigma);
survivingSpecies = abundanceHistory(end, :) > extinctionThreshold;

positiveAbundances = abundanceHistory(abundanceHistory > 0);
if isempty(positiveAbundances)
    error('The saved trajectory contains no positive abundances to plot.')
end
abundanceUpperLimit = max(positiveAbundances);
candidateAbundanceTicks = [1e-2, 1e0, 1e2];
abundanceTicks = candidateAbundanceTicks( ...
    candidateAbundanceTicks <= abundanceUpperLimit);

figureHandle = figure(1);
clf(figureHandle)
set(figureHandle, 'Color', 'w', 'Position', [50 80 1900 900], ...
    'Renderer', 'painters')
layoutHandle = tiledlayout(figureHandle, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Left panel: abundance dynamics for every species.
abundanceAxes = nexttile(layoutHandle, 1);
hold(abundanceAxes, 'on')
abundanceForPlot = abundanceHistory;
abundanceForPlot(abundanceForPlot <= 0) = NaN;
for speciesIndex = 1:nSpecies
    semilogy(abundanceAxes, timeValues, ...
        abundanceForPlot(:, speciesIndex), ...
        'Color', speciesColors(speciesIndex, :), ...
        'LineWidth', curveLineWidth)
end
xlabel(abundanceAxes, 'Time')
ylabel(abundanceAxes, '$N_i(t)$')
title(abundanceAxes, panelTitle)
xlim(abundanceAxes, [tMinimum tMaximum])
xticks(abundanceAxes, ...
    [tMinimum, (tMinimum + tMaximum) ./ 2, tMaximum])
ylim(abundanceAxes, [extinctionThreshold abundanceUpperLimit])
yticks(abundanceAxes, abundanceTicks)
set(abundanceAxes, 'YScale', 'log', 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth)
pbaspect(abundanceAxes, [1 1 1])
box(abundanceAxes, 'on')

% Right panel: theta dynamics only for species surviving at the final time.
% The same speciesColors row is used in both panels, so N_i and theta_i have
% identical colors for every plotted species.
thetaAxes = nexttile(layoutHandle, 2);
hold(thetaAxes, 'on')
survivingIndices = find(survivingSpecies);
for speciesIndex = survivingIndices
    plot(thetaAxes, timeValues, thetaHistory(:, speciesIndex), ...
        'Color', speciesColors(speciesIndex, :), ...
        'LineWidth', curveLineWidth)
end
xlabel(thetaAxes, 'Time')
ylabel(thetaAxes, '$\theta_i(t)$')
title(thetaAxes, panelTitle)
xlim(thetaAxes, [tMinimum tMaximum])
xticks(thetaAxes, [tMinimum, (tMinimum + tMaximum) ./ 2, tMaximum])
ylim(thetaAxes, [0 1])
set(thetaAxes, 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth)
pbaspect(thetaAxes, [1 1 1])
box(thetaAxes, 'on')

fprintf(['Plotted %d saved species trajectories from %s; %d species ' ...
    'survive at the final time.\n'], nSpecies, inputWorkbookPath, ...
    nnz(survivingSpecies))
