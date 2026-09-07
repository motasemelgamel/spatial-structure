% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Plot the stored two-species results without rerunning the simulation.

% Locate susceptibility.xlsx beside the saved plotting script. Additional
% candidates allow normal execution from the Editor and MATLAB path.
inputFileName = 'susceptibility.xlsx';
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
resolvedScriptPath = which('plot_two_species_susceptibility');
if ~isempty(resolvedScriptPath)
    candidateDirectories{end + 1} = fileparts(resolvedScriptPath);
end
candidateDirectories{end + 1} = pwd;
candidateDirectories = unique(candidateDirectories, 'stable');

scriptDirectory = '';
for directoryIndex = 1:numel(candidateDirectories)
    candidateDirectory = candidateDirectories{directoryIndex};
    if isfile(fullfile(candidateDirectory, inputFileName))
        scriptDirectory = candidateDirectory;
        break
    end
end

if isempty(scriptDirectory)
    error('Could not locate %s. Checked: %s', inputFileName, ...
        strjoin(candidateDirectories, ', '))
end
workbookPath = fullfile(scriptDirectory, inputFileName);

% Read the three plotted matrices and their beta and theta coordinates.
[betaValues, thetaValues, susceptibilityForPlot] = ...
    readHeatmapSheet(workbookPath, 'Figure1_susceptibility');
[betaForN1, thetaForN1, logN1ForPlot] = ...
    readHeatmapSheet(workbookPath, 'Figure2_N1');
[betaForN2, thetaForN2, logN2ForPlot] = ...
    readHeatmapSheet(workbookPath, 'Figure2_N2');

% Require identical parameter grids across all three workbook sheets.
assertSameGrid(betaValues, betaForN1, 'beta', 'Figure2_N1')
assertSameGrid(betaValues, betaForN2, 'beta', 'Figure2_N2')
assertSameGrid(thetaValues, thetaForN1, 'theta', 'Figure2_N1')
assertSameGrid(thetaValues, thetaForN2, 'theta', 'Figure2_N2')

if any(betaValues <= 0)
    error('All beta values must be positive before applying log10(beta).')
end
logBetaValues = log10(betaValues);

% Use one abundance color scale so the two species remain comparable.
abundanceColorLimits = [min([logN1ForPlot(:); logN2ForPlot(:)]), ...
    max([logN1ForPlot(:); logN2ForPlot(:)])];
if ~all(isfinite(abundanceColorLimits))
    error('The stored abundance heatmaps contain nonfinite values.')
end
if abundanceColorLimits(1) == abundanceColorLimits(2)
    abundanceColorLimits = abundanceColorLimits + [-0.5, 0.5];
end

% Apply the same presentation format as the simulation script.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
figureFontSize = 75;
titleFontSize = 60;
colorbarFontSize = 60;
axesLineWidth = 1.5;

% Figure 1: normalized susceptibility coefficient.
figure(1)
clf
set(gcf, 'Color', 'w', 'Position', [180 100 900 850])
susceptibilityAxes = axes;
formatHeatmapAxes(susceptibilityAxes, thetaValues, logBetaValues, ...
    susceptibilityForPlot, figureFontSize, axesLineWidth)
title(susceptibilityAxes, ...
    '$\log_{10}(|\chi_{11}|/\min |\chi_{11}|)$', ...
    'FontSize', titleFontSize)
formatHeatmapColorbar(susceptibilityAxes, colorbarFontSize, axesLineWidth)

% Figure 2: equilibrium abundance of each species.
figure(2)
clf
set(gcf, 'Color', 'w', 'Position', [50 80 1900 900])
abundanceLayout = tiledlayout(1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

firstAbundanceAxes = nexttile(abundanceLayout, 1);
formatHeatmapAxes(firstAbundanceAxes, thetaValues, logBetaValues, ...
    logN1ForPlot, figureFontSize, axesLineWidth)
clim(firstAbundanceAxes, abundanceColorLimits)
title(firstAbundanceAxes, '$\log_{10}(N_1)$', ...
    'FontSize', titleFontSize)

secondAbundanceAxes = nexttile(abundanceLayout, 2);
formatHeatmapAxes(secondAbundanceAxes, thetaValues, logBetaValues, ...
    logN2ForPlot, figureFontSize, axesLineWidth)
clim(secondAbundanceAxes, abundanceColorLimits)
title(secondAbundanceAxes, '$\log_{10}(N_2)$', ...
    'FontSize', titleFontSize)
abundanceColorbar = formatHeatmapColorbar( ...
    secondAbundanceAxes, colorbarFontSize, axesLineWidth);
abundanceColorbar.Layout.Tile = 'east';

fprintf('Plotted the stored susceptibility results from %s.\n', workbookPath)

function [betaValues, thetaValues, heatmapValues] = ...
    readHeatmapSheet(workbookPath, sheetName)

% Read one beta column followed by heatmap columns named theta=<value>.
dataTable = readtable(workbookPath, 'Sheet', sheetName, ...
    'VariableNamingRule', 'preserve');
variableNames = string(dataTable.Properties.VariableNames);

if width(dataTable) < 2 || ~strcmpi(variableNames(1), 'beta') || ...
        any(~startsWith(variableNames(2:end), 'theta='))
    error(['Sheet %s must contain beta followed by theta=<value> ' ...
        'columns.'], sheetName)
end

betaValues = dataTable{:, 1};
thetaValues = str2double(extractAfter(variableNames(2:end), 'theta='));
heatmapValues = dataTable{:, 2:end};

if ~isnumeric(betaValues) || ~isnumeric(heatmapValues) || ...
        any(~isfinite(betaValues), 'all') || ...
        any(~isfinite(thetaValues), 'all') || ...
        any(~isfinite(heatmapValues), 'all')
    error('Sheet %s must contain only finite numeric plot data.', sheetName)
end
betaValues = betaValues(:);
thetaValues = thetaValues(:).';

if size(heatmapValues, 1) ~= numel(betaValues) || ...
        size(heatmapValues, 2) ~= numel(thetaValues)
    error('The heatmap dimensions are inconsistent in sheet %s.', sheetName)
end
if any(diff(betaValues) <= 0) || any(diff(thetaValues) <= 0)
    error('Beta and theta must be strictly increasing in sheet %s.', sheetName)
end
end

function assertSameGrid(referenceValues, candidateValues, ...
    parameterName, sheetName)

% Compare parameter grids up to floating-point roundoff.
referenceValues = referenceValues(:);
candidateValues = candidateValues(:);
if numel(referenceValues) ~= numel(candidateValues)
    error('%s-grid length differs in sheet %s.', parameterName, sheetName)
end
comparisonScale = max(1, max(abs(referenceValues), [], 'omitnan'));
tolerance = 100 .* eps(comparisonScale);
if any(abs(candidateValues - referenceValues) > tolerance)
    error('%s values differ in sheet %s.', parameterName, sheetName)
end
end

function formatHeatmapAxes(axesHandle, thetaValues, logBetaValues, ...
    heatmapValues, figureFontSize, axesLineWidth)

% Plot uniformly sampled log10(beta) coordinates without axis distortion.
imagesc(axesHandle, thetaValues, logBetaValues, heatmapValues)
set(axesHandle, 'YDir', 'normal', 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth, 'TickDir', 'out', 'Layer', 'top')
xlabel(axesHandle, '$\theta$')
ylabel(axesHandle, '$\beta$')
xlim(axesHandle, [min(thetaValues), max(thetaValues)])
ylim(axesHandle, [min(logBetaValues), max(logBetaValues)])
xticks(axesHandle, [0.05, 0.5, 0.95])

minimumLogBeta = floor(min(logBetaValues));
maximumLogBeta = ceil(max(logBetaValues));
logBetaTicks = minimumLogBeta:2:maximumLogBeta;
if logBetaTicks(end) < maximumLogBeta
    logBetaTicks(end + 1) = maximumLogBeta;
end
yticks(axesHandle, logBetaTicks)
yticklabels(axesHandle, compose('$10^{%d}$', logBetaTicks))
colormap(axesHandle, parula(256))
pbaspect(axesHandle, [1 1 1])
box(axesHandle, 'on')
end

function colorbarHandle = formatHeatmapColorbar( ...
    axesHandle, colorbarFontSize, axesLineWidth)

% Apply the established typography and line weight to one colorbar.
colorbarHandle = colorbar(axesHandle);
colorbarHandle.FontSize = colorbarFontSize;
colorbarHandle.LineWidth = axesLineWidth;
colorbarHandle.TickLabelInterpreter = 'latex';
end
