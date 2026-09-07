% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Locate both workbooks robustly even when MATLAB executes an Editor_* copy.
inputFileNames = { ...
    'example_abundance_trajectories.xlsx'; ...
    'burn_surviving_fraction.xlsx'};
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
resolvedScriptPath = which('plot_Sfig');
if ~isempty(resolvedScriptPath)
    candidateDirectories{end + 1} = fileparts(resolvedScriptPath);
end
candidateDirectories{end + 1} = pwd;
candidateDirectories = unique(candidateDirectories, 'stable');

scriptDirectory = '';
for directoryIndex = 1:numel(candidateDirectories)
    candidateDirectory = candidateDirectories{directoryIndex};
    containsAllFiles = all(cellfun(@(fileName) ...
        isfile(fullfile(candidateDirectory, fileName)), inputFileNames));
    if containsAllFiles
        scriptDirectory = candidateDirectory;
        break
    end
end

if isempty(scriptDirectory)
    error(['Could not locate both required workbooks. Checked: %s'], ...
        strjoin(candidateDirectories, ', '))
end
abundanceWorkbookPath = fullfile(scriptDirectory, inputFileNames{1});
burnWorkbookPath = fullfile(scriptDirectory, inputFileNames{2});

% Read every worksheet and retain its time, phase, and abundance data.
sheetNames = string(sheetnames(abundanceWorkbookPath));
nPanels = numel(sheetNames);
if nPanels == 0
    error('The workbook contains no worksheets to plot.')
end
if nPanels ~= 5
    error(['A 2-by-3 layout requires exactly five abundance sheets, ' ...
        'but the workbook contains %d.'], nPanels)
end

timeByPanel = cell(nPanels, 1);
phaseByPanel = cell(nPanels, 1);
abundanceByPanel = cell(nPanels, 1);
parameterSymbols = strings(nPanels, 1);
parameterValues = nan(nPanels, 1);
betaValues = nan(nPanels, 1);

globalPositiveMinimum = inf;
globalPositiveMaximum = 0;
globalTimeMinimum = inf;
globalTimeMaximum = -inf;
nSpecies = [];

for panelIndex = 1:nPanels
    sheetName = sheetNames(panelIndex);
    dataTable = readtable(abundanceWorkbookPath, ...
        'Sheet', char(sheetName), ...
        'VariableNamingRule', 'preserve');
    variableNames = string(dataTable.Properties.VariableNames);

    if width(dataTable) < 3 || ...
            ~strcmpi(variableNames(1), 'time') || ...
            ~strcmpi(variableNames(2), 'phase') || ...
            any(~startsWith(variableNames(3:end), 'species_'))
        error(['Sheet %s must contain time, phase, and species_* ' ...
            'columns in that order.'], sheetName)
    end

    timeValues = dataTable{:, 1};
    phaseValues = string(dataTable{:, 2});
    abundanceValues = dataTable{:, 3:end};
    if ~isnumeric(timeValues) || ~isnumeric(abundanceValues)
        error('Time and species abundances must be numeric in sheet %s.', ...
            sheetName)
    end
    timeValues = timeValues(:);

    if numel(phaseValues) ~= numel(timeValues) || ...
            size(abundanceValues, 1) ~= numel(timeValues)
        error('Column lengths are inconsistent in sheet %s.', sheetName)
    end
    if any(diff(timeValues) < 0)
        error('Time must be nondecreasing in sheet %s.', sheetName)
    end

    if panelIndex == 1
        nSpecies = size(abundanceValues, 2);
    elseif size(abundanceValues, 2) ~= nSpecies
        error('Species count is inconsistent in sheet %s.', sheetName)
    end

    [parameterSymbols(panelIndex), parameterValues(panelIndex), ...
        betaValues(panelIndex)] = parseSheetParameters(sheetName);

    timeByPanel{panelIndex} = timeValues;
    phaseByPanel{panelIndex} = phaseValues;
    abundanceByPanel{panelIndex} = abundanceValues;

    positiveValues = abundanceValues( ...
        isfinite(abundanceValues) & abundanceValues > 0);
    if ~isempty(positiveValues)
        globalPositiveMinimum = min(globalPositiveMinimum, ...
            min(positiveValues));
        globalPositiveMaximum = max(globalPositiveMaximum, ...
            max(positiveValues));
    end
    globalTimeMinimum = min(globalTimeMinimum, min(timeValues));
    globalTimeMaximum = max(globalTimeMaximum, max(timeValues));
end

% Read the deterministic burn-in survival curves formerly shown as Figure 2
% by plot_fig.m. Each theta value is stored in a separate worksheet.
burnSheetNames = string(sheetnames(burnWorkbookPath));
burnThetaValues = str2double(extractAfter(burnSheetNames, 'theta='));
if isempty(burnSheetNames) || any(~isfinite(burnThetaValues))
    error(['Every burn-in worksheet name must follow the form ' ...
        'theta=<numeric value>.'])
end
[burnThetaValues, burnThetaOrder] = sort(burnThetaValues);
burnSheetNames = burnSheetNames(burnThetaOrder);
nBurnTheta = numel(burnThetaValues);
burnBetaValues = [];
burnSurvivingFraction = [];

for thetaIndex = 1:nBurnTheta
    [currentBetaValues, currentBurnFraction] = readBurnSheet( ...
        burnWorkbookPath, burnSheetNames(thetaIndex));

    if thetaIndex == 1
        burnBetaValues = currentBetaValues;
        burnSurvivingFraction = nan( ...
            numel(burnBetaValues), nBurnTheta);
    else
        assertSameBeta(burnBetaValues, currentBetaValues, ...
            sprintf('%s in burn-in data', burnSheetNames(thetaIndex)))
    end
    burnSurvivingFraction(:, thetaIndex) = currentBurnFraction;
end

if ~isfinite(globalPositiveMinimum) || globalPositiveMaximum <= 0
    error('The workbook contains no positive species abundances to plot.')
end

% Use common decade limits so abundance magnitudes are comparable by panel.
abundanceLowerLimit = 10 .^ floor(log10(globalPositiveMinimum));
abundanceUpperLimit = 10 .^ ceil(log10(globalPositiveMaximum));
if abundanceUpperLimit <= abundanceLowerLimit
    abundanceUpperLimit = 10 .* abundanceLowerLimit;
end

% Apply the same scientific figure formatting used in plot_fig.m.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
figureFontSize = 60;
titleFontSize = 42;
axesLineWidth = 1.5;
speciesLineWidth = 3;
speciesColors = lines(nSpecies);
burnColors = lines(nBurnTheta);
burnLegendLabels = compose('$\\theta = %.2g$', burnThetaValues);

burnPanelMinimum = min( ...
    burnSurvivingFraction, [], 'all', 'omitnan');
if ~isfinite(burnPanelMinimum) || burnPanelMinimum >= 1
    burnPanelMinimum = 0.9;
end

% Combine the burn-in richness panel and five abundance panels in a 2-by-3
% supplementary figure.
totalPanels = nPanels + 1;
nColumns = 3;
nRows = ceil(totalPanels ./ nColumns);
figureHandle = figure(1);
clf(figureHandle)
set(figureHandle, 'Color', 'w', 'Position', [50 40 2200 1400], ...
    'Renderer', 'painters')
layoutHandle = tiledlayout(nRows, nColumns, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: surviving fraction at the end of deterministic burn-in.
burnAxes = nexttile(layoutHandle, 1);
hold(burnAxes, 'on')
burnFigureHandles = gobjects(nBurnTheta, 1);
for thetaIndex = 1:nBurnTheta
    burnFigureHandles(thetaIndex) = plot( ...
        burnAxes, burnBetaValues, ...
        burnSurvivingFraction(:, thetaIndex), ...
        'Color', burnColors(thetaIndex, :), 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerSize', 15, ...
        'MarkerFaceColor', burnColors(thetaIndex, :), ...
        'MarkerEdgeColor', burnColors(thetaIndex, :), ...
        'LineWidth', 3);
end
set(burnAxes, 'XScale', 'log', 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth)
xlabel(burnAxes, '$\beta$')
ylabel(burnAxes, 'Fraction of surviving species')
xlim(burnAxes, [min(burnBetaValues) max(burnBetaValues)])
ylim(burnAxes, [burnPanelMinimum 1])
xticks(burnAxes, [1e-2 1 1e2])
pbaspect(burnAxes, [1 1 1])
box(burnAxes, 'on')
title(burnAxes, 'Before noise', 'FontSize', titleFontSize)
legend(burnAxes, burnFigureHandles, burnLegendLabels, ...
    'Location', 'best', 'Orientation', 'vertical', 'FontSize', 42)

for panelIndex = 1:nPanels
    axesHandle = nexttile(layoutHandle, panelIndex + 1);
    colororder(axesHandle, speciesColors)
    semilogy(axesHandle, timeByPanel{panelIndex}, ...
        abundanceByPanel{panelIndex}, 'LineWidth', speciesLineWidth)
    hold(axesHandle, 'on')

    % Mark the transition from deterministic burn-in to noisy observation.
    observationIndex = find( ...
        phaseByPanel{panelIndex} == "observation", 1, 'first');
    if ~isempty(observationIndex)
        if observationIndex > 1
            transitionTime = timeByPanel{panelIndex}(observationIndex - 1);
        else
            transitionTime = timeByPanel{panelIndex}(observationIndex);
        end
        xline(axesHandle, transitionTime, '--', ...
            'Color', [0.15 0.15 0.15], 'LineWidth', axesLineWidth, ...
            'HandleVisibility', 'off');
    end

    set(axesHandle, 'YScale', 'log', 'FontSize', figureFontSize, ...
        'LineWidth', axesLineWidth)
    xlabel(axesHandle, 'Time')
    ylabel(axesHandle, 'Species abundance')
    xlim(axesHandle, [globalTimeMinimum globalTimeMaximum])
    ylim(axesHandle, [abundanceLowerLimit abundanceUpperLimit])
    pbaspect(axesHandle, [1 1 1])
    box(axesHandle, 'on')

    titleText = sprintf( ...
        '$%s = %.2g,\\; \\beta = %.3g$', ...
        parameterSymbols(panelIndex), parameterValues(panelIndex), ...
        betaValues(panelIndex));
    title(axesHandle, titleText, 'FontSize', titleFontSize)
end

% Save the completed figure as editable vector graphics rather than a
% raster image embedded inside the SVG file.
svgOutputPath = fullfile(scriptDirectory, 'plot_Sfig.svg');
print(figureHandle, svgOutputPath, '-dsvg', '-painters')

fprintf(['Plotted one burn-in panel and %d abundance panels containing ' ...
    '%d species from %s.\nSaved the vector figure to %s.\n'], ...
    nPanels, nSpecies, scriptDirectory, svgOutputPath)

function [betaValues, burnFraction] = ...
    readBurnSheet(workbookPath, sheetName)

% Read the deterministic burn-in surviving fraction for one theta value.
dataTable = readtable(workbookPath, 'Sheet', char(sheetName), ...
    'VariableNamingRule', 'preserve');
variableNames = string(dataTable.Properties.VariableNames);
if width(dataTable) ~= 2 || ~strcmpi(variableNames(1), 'beta') || ...
        ~strcmpi(variableNames(2), 'S_burn_over_S0')
    error(['Sheet %s in %s must contain exactly beta and ' ...
        'S_burn_over_S0.'], sheetName, workbookPath)
end

betaValues = dataTable{:, 1};
burnFraction = dataTable{:, 2};
if ~isnumeric(betaValues) || ~isnumeric(burnFraction) || ...
        any(~isfinite(betaValues), 'all') || ...
        any(~isfinite(burnFraction), 'all')
    error('Sheet %s in %s must contain finite numeric data.', ...
        sheetName, workbookPath)
end
betaValues = betaValues(:);
burnFraction = burnFraction(:);
end

function assertSameBeta(referenceBeta, candidateBeta, contextText)

% Require identical beta grids up to floating-point roundoff.
candidateBeta = candidateBeta(:);
if numel(candidateBeta) ~= numel(referenceBeta)
    error('Beta-grid length mismatch for %s.', contextText)
end
comparisonScale = max(1, max(abs(referenceBeta), [], 'omitnan'));
tolerance = 100 .* eps(comparisonScale);
if any(abs(candidateBeta - referenceBeta) > tolerance)
    error('Beta values do not match for %s.', contextText)
end
end

function [parameterSymbol, parameterValue, betaValue] = ...
    parseSheetParameters(sheetName)

% Parse names such as theta=0.5, beta=3.29 without hard-coding values.
tokens = regexp(char(sheetName), ...
    '^(theta|eta)=([^,]+),\s*beta=(.+)$', 'tokens', 'once');
if isempty(tokens)
    error(['Sheet %s must follow either theta=<value>, beta=<value> ' ...
        'or eta=<value>, beta=<value>.'], sheetName)
end

parameterName = tokens{1};
parameterValue = str2double(tokens{2});
betaValue = str2double(tokens{3});
if ~isfinite(parameterValue) || ~isfinite(betaValue)
    error('Sheet %s contains a nonnumeric parameter value.', sheetName)
end

if strcmp(parameterName, 'theta')
    parameterSymbol = "\theta";
else
    parameterSymbol = "\eta";
end
end
