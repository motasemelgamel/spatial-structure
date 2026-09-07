% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Reproduce the eta and sigma sweeps using only the saved Excel statistics.
inputFileName = 'fraction_of_species_with_noise_and_adaptation.xlsx';

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
resolvedScriptPath = which( ...
    'plot_fraction_of_species_with_noise_and_adaptation');
if ~isempty(resolvedScriptPath)
    candidateDirectories{end + 1} = fileparts(resolvedScriptPath);
end
candidateDirectories{end + 1} = pwd;
candidateDirectories = unique(candidateDirectories, 'stable');

scriptDirectory = '';
for directoryIndex = 1:numel(candidateDirectories)
    candidateWorkbookPath = fullfile( ...
        candidateDirectories{directoryIndex}, inputFileName);
    candidateScriptPath = fullfile(candidateDirectories{directoryIndex}, ...
        'plot_fraction_of_species_with_noise_and_adaptation.m');
    if isfile(candidateWorkbookPath) && isfile(candidateScriptPath)
        scriptDirectory = candidateDirectories{directoryIndex};
        break
    end
end
if isempty(scriptDirectory)
    error(['Could not locate %s beside ' ...
        'plot_fraction_of_species_with_noise_and_adaptation.m.'], ...
        inputFileName)
end
inputWorkbookPath = fullfile(scriptDirectory, inputFileName);

% Read the two plain sweep sheets. Each mean column is immediately followed
% by its corresponding sample-standard-deviation column.
etaTable = readtable(inputWorkbookPath, 'Sheet', 'Eta sweep', ...
    'VariableNamingRule', 'preserve');
sigmaTable = readtable(inputWorkbookPath, 'Sheet', 'Sigma sweep', ...
    'VariableNamingRule', 'preserve');

etaNames = string(etaTable.Properties.VariableNames);
sigmaNames = string(sigmaTable.Properties.VariableNames);
etaMeanColumns = find(startsWith(etaNames, 'mean_sigma_'));
etaStdColumns = find(startsWith(etaNames, 'std_sigma_'));
sigmaMeanColumns = find(startsWith(sigmaNames, 'mean_eta_'));
sigmaStdColumns = find(startsWith(sigmaNames, 'std_eta_'));

if isempty(etaMeanColumns) || ...
        numel(etaMeanColumns) ~= numel(etaStdColumns) || ...
        isempty(sigmaMeanColumns) || ...
        numel(sigmaMeanColumns) ~= numel(sigmaStdColumns)
    error('The sweep workbook has incomplete mean or standard-deviation data.')
end

etaSweepValues = etaTable{:, 'eta'};
etaSweepMean = etaTable{:, etaMeanColumns};
etaSweepStd = etaTable{:, etaStdColumns};
sigmaForEtaSweep = parseParameterValues( ...
    etaNames(etaMeanColumns), 'mean_sigma_');

sigmaSweepValues = sigmaTable{:, 'sigma'};
sigmaSweepMean = sigmaTable{:, sigmaMeanColumns};
sigmaSweepStd = sigmaTable{:, sigmaStdColumns};
etaForSigmaSweep = parseParameterValues( ...
    sigmaNames(sigmaMeanColumns), 'mean_eta_');

% Confirm that the paired mean and standard-deviation headers describe the
% same fixed parameter values in the same column order.
sigmaFromStdHeaders = parseParameterValues( ...
    etaNames(etaStdColumns), 'std_sigma_');
etaFromStdHeaders = parseParameterValues( ...
    sigmaNames(sigmaStdColumns), 'std_eta_');
if ~isequal(sigmaForEtaSweep, sigmaFromStdHeaders) || ...
        ~isequal(etaForSigmaSweep, etaFromStdHeaders)
    error('Mean and standard-deviation headers are paired inconsistently.')
end

validateSweep(etaSweepValues, etaSweepMean, etaSweepStd, 'eta')
validateSweep(sigmaSweepValues, sigmaSweepMean, sigmaSweepStd, 'sigma')

% Apply the same two-panel plotting format used by the original script.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
figureFontSize = 60;
legendFontSize = 36;
markerSize = 8;
curveLineWidth = 3;
axesLineWidth = 1.5;
curveColors = lines(max(numel(sigmaForEtaSweep), ...
    numel(etaForSigmaSweep)));

sweepFigure = figure(1);
clf(sweepFigure)
set(sweepFigure, 'Color', 'w', 'Position', [50 80 1900 900], ...
    'Renderer', 'painters')
sweepLayout = tiledlayout(sweepFigure, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

etaSweepAxes = nexttile(sweepLayout, 1);
plotOverlaySweep(etaSweepAxes, etaSweepValues, etaSweepMean, etaSweepStd, ...
    sigmaForEtaSweep, 'sigma', '$\eta$', [0, 0.5, 1], ...
    curveColors, figureFontSize, legendFontSize, markerSize, ...
    curveLineWidth, axesLineWidth)

sigmaSweepAxes = nexttile(sweepLayout, 2);
plotOverlaySweep(sigmaSweepAxes, sigmaSweepValues, sigmaSweepMean, ...
    sigmaSweepStd, etaForSigmaSweep, 'eta', '$\sigma$', [0, 0.5, 1], ...
    curveColors, figureFontSize, legendFontSize, markerSize, ...
    curveLineWidth, axesLineWidth)

fprintf('Plotted the saved eta and sigma sweeps from %s.\n', ...
    inputWorkbookPath)


function parameterValues = parseParameterValues(variableNames, prefix)

% Convert plain Excel header tokens such as 0p1 back to numeric values.
parameterTokens = extractAfter(variableNames, prefix);
parameterValues = str2double(strrep(parameterTokens, 'p', '.'));
parameterValues = parameterValues(:).';
if any(~isfinite(parameterValues))
    error('Could not parse the fixed parameter values from the headers.')
end
end


function validateSweep(xValues, meanValues, standardDeviations, sweepName)

% Check dimensions, ordering, and the allowed range of survival statistics.
xValues = xValues(:);
if isempty(xValues) || any(~isfinite(xValues)) || ...
        any(diff(xValues) <= 0) || ...
        size(meanValues, 1) ~= numel(xValues) || ...
        ~isequal(size(meanValues), size(standardDeviations)) || ...
        any(~isfinite(meanValues), 'all') || ...
        any(~isfinite(standardDeviations), 'all') || ...
        any(meanValues < 0 | meanValues > 1, 'all') || ...
        any(standardDeviations < 0, 'all')
    error('The saved %s sweep data are invalid.', sweepName)
end
end


function plotOverlaySweep(axesHandle, xValues, ...
    meanValues, standardDeviations, fixedValues, fixedParameterName, ...
    xAxisLabel, xTicks, curveColors, figureFontSize, legendFontSize, ...
    markerSize, curveLineWidth, axesLineWidth)

% Plot saved means with one-sample-standard-deviation error bars.
expectedSize = [numel(xValues), numel(fixedValues)];
if ~isequal(size(meanValues), expectedSize) || ...
        ~isequal(size(standardDeviations), expectedSize)
    error('The saved sweep statistics have inconsistent dimensions.')
end

hold(axesHandle, 'on')
curveHandles = gobjects(numel(fixedValues), 1);
legendLabels = cell(numel(fixedValues), 1);
for curveIndex = 1:numel(fixedValues)
    curveHandles(curveIndex) = errorbar(axesHandle, xValues, ...
        meanValues(:, curveIndex), standardDeviations(:, curveIndex), ...
        'Color', curveColors(curveIndex, :), 'LineStyle', '-', ...
        'Marker', 'o', 'MarkerSize', markerSize, ...
        'MarkerFaceColor', curveColors(curveIndex, :), ...
        'MarkerEdgeColor', curveColors(curveIndex, :), ...
        'LineWidth', curveLineWidth, 'CapSize', 8);

    if strcmp(fixedParameterName, 'sigma')
        legendLabels{curveIndex} = sprintf( ...
            '$\\sigma = %.3g$', fixedValues(curveIndex));
    elseif strcmp(fixedParameterName, 'eta')
        legendLabels{curveIndex} = sprintf( ...
            '$\\eta = %.3g$', fixedValues(curveIndex));
    else
        error('Unknown fixed parameter name: %s.', fixedParameterName)
    end
end

survivingFractionMinimum = min( ...
    meanValues - standardDeviations, [], 'all', 'omitnan');
if ~isfinite(survivingFractionMinimum) || survivingFractionMinimum >= 1
    survivingFractionMinimum = 0.9;
end

xlabel(axesHandle, xAxisLabel)
ylabel(axesHandle, 'Fraction of surviving species')
xlim(axesHandle, [min(xValues), max(xValues)])
ylim(axesHandle, [survivingFractionMinimum, 1])
xticks(axesHandle, xTicks)
set(axesHandle, 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth)
pbaspect(axesHandle, [1 1 1])
box(axesHandle, 'on')
legend(axesHandle, curveHandles, legendLabels, ...
    'Location', 'best', 'FontSize', legendFontSize, 'Box', 'on')
end
