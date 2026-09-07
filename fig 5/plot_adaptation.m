% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Reproduce the adaptation figures using only the saved Excel trajectories.
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
resolvedScriptPath = which('plot_adaptation');
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
        candidateDirectories{directoryIndex}, 'plot_adaptation.m');
    if isfile(candidateWorkbookPath) && isfile(candidateScriptPath)
        scriptDirectory = candidateDirectories{directoryIndex};
        break
    end
end
if isempty(scriptDirectory)
    error('Could not locate %s beside plot_adaptation.m.', inputFileName)
end
inputWorkbookPath = fullfile(scriptDirectory, inputFileName);

% Read and order the eta worksheets by their numeric adaptation rates.
sheetNames = string(sheetnames(inputWorkbookPath));
etaValues = str2double(extractAfter(sheetNames, 'eta='));
if isempty(sheetNames) || any(~isfinite(etaValues))
    error('Every worksheet name must follow the form eta=<numeric value>.')
end
[etaValues, etaOrder] = sort(etaValues);
sheetNames = sheetNames(etaOrder);
nEta = numel(etaValues);

timeByEta = cell(nEta, 1);
abundanceByEta = cell(nEta, 1);
thetaByEta = cell(nEta, 1);
survivingSpeciesByEta = cell(nEta, 1);
nSpecies = [];
globalAbundanceMaximum = 0;

% Extract the abundance and theta blocks from every eta worksheet.
for etaIndex = 1:nEta
    sheetName = sheetNames(etaIndex);
    dataTable = readtable(inputWorkbookPath, ...
        'Sheet', char(sheetName), 'VariableNamingRule', 'preserve');
    variableNames = string(dataTable.Properties.VariableNames);

    abundanceColumns = find(startsWith(variableNames, 'abundance_'));
    thetaColumns = find(startsWith(variableNames, 'theta_'));
    if isempty(abundanceColumns) || isempty(thetaColumns) || ...
            numel(abundanceColumns) ~= numel(thetaColumns)
        error(['Sheet %s must contain matching abundance_* and theta_* ' ...
            'columns.'], sheetName)
    end

    currentSpeciesCount = numel(abundanceColumns);
    if etaIndex == 1
        nSpecies = currentSpeciesCount;
    elseif currentSpeciesCount ~= nSpecies
        error('The species count is inconsistent in sheet %s.', sheetName)
    end

    % The saved layout places abundance time first and theta time directly
    % before theta_1, with a blank separator between the two data blocks.
    abundanceTime = dataTable{:, 1};
    thetaTimeColumn = thetaColumns(1) - 1;
    thetaTime = dataTable{:, thetaTimeColumn};
    abundanceValues = dataTable{:, abundanceColumns};
    thetaValues = dataTable{:, thetaColumns};

    if ~isnumeric(abundanceTime) || ~isnumeric(thetaTime) || ...
            ~isnumeric(abundanceValues) || ~isnumeric(thetaValues)
        error('Sheet %s contains nonnumeric trajectory data.', sheetName)
    end
    abundanceTime = abundanceTime(:);
    thetaTime = thetaTime(:);
    if any(~isfinite(abundanceTime)) || any(~isfinite(thetaTime)) || ...
            numel(abundanceTime) ~= numel(thetaTime) || ...
            any(abs(abundanceTime - thetaTime) > ...
            100 .* eps(max(1, max(abs(abundanceTime), [], 'omitnan'))))
        error('Abundance and theta time grids differ in sheet %s.', sheetName)
    end
    if any(diff(abundanceTime) < 0) || ...
            any(~isfinite(abundanceValues), 'all') || ...
            any(~isfinite(thetaValues), 'all') || ...
            any(abundanceValues < 0, 'all') || ...
            any(thetaValues < 0 | thetaValues > 1, 'all')
        error('Sheet %s contains an invalid saved trajectory.', sheetName)
    end

    timeByEta{etaIndex} = abundanceTime;
    abundanceByEta{etaIndex} = abundanceValues;
    thetaByEta{etaIndex} = thetaValues;
    survivingSpeciesByEta{etaIndex} = ...
        abundanceValues(end, :) > 1e-3;
    globalAbundanceMaximum = max(globalAbundanceMaximum, ...
        max(abundanceValues, [], 'all'));
end

% Match the formatting and biological plotting rules in adaptation.m.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
figureFontSize = 100;
curveLineWidth = 3;
axesLineWidth = 1.5;
figurePosition = [50 80 1800 850];
extinctionThreshold = 1e-3;
speciesColors = lines(nSpecies);
candidateAbundanceTicks = [1e-2, 1e0, 1e2];
abundanceTicks = candidateAbundanceTicks( ...
    candidateAbundanceTicks <= globalAbundanceMaximum);

% Reproduce one two-panel figure for every saved eta value.
for etaIndex = 1:nEta
    timeValues = timeByEta{etaIndex};
    tMinimum = min(timeValues);
    tMaximum = max(timeValues);

    figureHandle = figure(etaIndex);
    clf(figureHandle)
    set(figureHandle, 'Color', 'w', 'Position', figurePosition, ...
        'Renderer', 'painters')
    figureLayout = tiledlayout(1, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    abundanceAxes = nexttile(figureLayout, 1);
    plotAbundancePanel(abundanceAxes, timeValues, ...
        abundanceByEta{etaIndex}, etaValues(etaIndex), ...
        speciesColors, extinctionThreshold, globalAbundanceMaximum, ...
        abundanceTicks, tMinimum, tMaximum, figureFontSize, ...
        curveLineWidth, axesLineWidth)

    thetaAxes = nexttile(figureLayout, 2);
    plotThetaPanel(thetaAxes, timeValues, thetaByEta{etaIndex}, ...
        survivingSpeciesByEta{etaIndex}, etaValues(etaIndex), ...
        speciesColors, tMinimum, tMaximum, figureFontSize, ...
        curveLineWidth, axesLineWidth)

    % Match adaptation.m by saving each completed figure as vector SVG.
    svgOutputPath = fullfile(scriptDirectory, ...
        sprintf('adaptation_eta_%.3g.svg', etaValues(etaIndex)));
    print(figureHandle, svgOutputPath, '-dsvg', '-painters')
end

fprintf('Plotted %d eta cases containing %d species from %s.\n', ...
    nEta, nSpecies, inputWorkbookPath)

function plotAbundancePanel(axesHandle, timeValues, abundanceValues, ...
    eta, speciesColors, extinctionThreshold, abundanceUpperLimit, ...
    abundanceTicks, tMinimum, tMaximum, figureFontSize, ...
    curveLineWidth, axesLineWidth)

% Plot every species abundance on the original logarithmic y-axis.
hold(axesHandle, 'on')
abundanceForPlot = abundanceValues;
abundanceForPlot(abundanceForPlot <= 0) = NaN;

for speciesIndex = 1:size(abundanceForPlot, 2)
    semilogy(axesHandle, timeValues, ...
        abundanceForPlot(:, speciesIndex), ...
        'Color', speciesColors(speciesIndex, :), ...
        'LineWidth', curveLineWidth)
end

xlabel(axesHandle, 'Time')
ylabel(axesHandle, 'Species abundance')
title(axesHandle, sprintf('$\\eta = %.3g$', eta))
xlim(axesHandle, [tMinimum tMaximum])
xticks(axesHandle, [tMinimum, (tMinimum + tMaximum) ./ 2, tMaximum])
ylim(axesHandle, [extinctionThreshold abundanceUpperLimit])
yticks(axesHandle, abundanceTicks)
set(axesHandle, 'YScale', 'log', 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth)
pbaspect(axesHandle, [1 1 1])
box(axesHandle, 'on')
end

function plotThetaPanel(axesHandle, timeValues, thetaValues, ...
    survivingSpecies, eta, speciesColors, tMinimum, tMaximum, ...
    figureFontSize, curveLineWidth, axesLineWidth)

% Plot theta only for species surviving at the saved final time point.
hold(axesHandle, 'on')
survivingIndices = find(survivingSpecies);
for speciesIndex = survivingIndices
    plot(axesHandle, timeValues, thetaValues(:, speciesIndex), ...
        'Color', speciesColors(speciesIndex, :), ...
        'LineWidth', curveLineWidth)
end

xlabel(axesHandle, 'Time')
ylabel(axesHandle, '$\theta_i(t)$')
title(axesHandle, sprintf('$\\eta = %.3g$', eta))
xlim(axesHandle, [tMinimum tMaximum])
xticks(axesHandle, [tMinimum, (tMinimum + tMaximum) ./ 2, tMaximum])
ylim(axesHandle, [0 1])
set(axesHandle, 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth)
pbaspect(axesHandle, [1 1 1])
box(axesHandle, 'on')
end
