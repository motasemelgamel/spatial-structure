% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Reproduce the adaptation-noise heatmap using only the saved Excel data.
inputFileName = 'survival_heatmap.xlsx';

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
resolvedScriptPath = which('plot_heatmap_wtih_noise');
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
        'plot_heatmap_wtih_noise.m');
    if isfile(candidateWorkbookPath) && isfile(candidateScriptPath)
        scriptDirectory = candidateDirectories{directoryIndex};
        break
    end
end
if isempty(scriptDirectory)
    error('Could not locate %s beside plot_heatmap_wtih_noise.m.', ...
        inputFileName)
end
inputWorkbookPath = fullfile(scriptDirectory, inputFileName);

% The sheet stores eta across the first row and sigma down the first column,
% so the numeric block is the surviving fraction with rows=sigma, cols=eta.
sheetName = 'Figure data';
gridData = readmatrix(inputWorkbookPath, 'Sheet', sheetName);
if size(gridData, 1) < 2 || size(gridData, 2) < 2
    error('Sheet %s must contain an eta row, a sigma column, and a grid.', ...
        sheetName)
end

etaValues = gridData(1, 2:end);
sigmaValues = gridData(2:end, 1).';
survivingFraction = gridData(2:end, 2:end);

if ~isequal(size(survivingFraction), ...
        [numel(sigmaValues), numel(etaValues)]) || ...
        isempty(survivingFraction) || ...
        any(~isfinite(etaValues)) || any(~isfinite(sigmaValues)) || ...
        any(~isfinite(survivingFraction), 'all')
    error('The saved heatmap data are invalid or incomplete.')
end
if any(etaValues <= 0) || any(sigmaValues <= 0)
    error('Eta and sigma must be positive for logarithmic coordinates.')
end
if any(diff(etaValues) <= 0) || any(diff(sigmaValues) <= 0)
    error('Eta and sigma must both be strictly increasing.')
end
if any(survivingFraction < 0, 'all') || any(survivingFraction > 1, 'all')
    error('The surviving fraction must lie between 0 and 1.')
end

% Plot one square panel over the complete logarithmic parameter range. The
% color range ends at the observed maximum rounded up to 0.05.
plotAdaptationHeatmapFullRange(etaValues, sigmaValues, ...
    survivingFraction, 'Fraction of surviving species', 1, 0.05);

fprintf('Plotted a %dx%d eta-sigma grid from %s.\n', ...
    numel(sigmaValues), numel(etaValues), inputWorkbookPath)


function figureHandle = plotAdaptationHeatmapFullRange( ...
    etaValues, sigmaValues, heatmapValues, colorbarLabel, ...
    figureNumber, colorStep)
%PLOTADAPTATIONHEATMAPFULLRANGE Plot the stored logarithmic parameter sweep.

% Check that the result matrix follows rows=sigma and columns=eta.
expectedSize = [numel(sigmaValues), numel(etaValues)];
if ~isequal(size(heatmapValues), expectedSize)
    error(['heatmapValues must have %d sigma rows and %d eta columns.'], ...
        expectedSize(1), expectedSize(2))
end
if ~isscalar(colorStep) || ~isfinite(colorStep) || colorStep <= 0
    error('colorStep must be a finite positive scalar.')
end

% End the color scale at a rounded observed maximum to reveal variation.
finiteValues = heatmapValues(isfinite(heatmapValues));
if isempty(finiteValues)
    error('heatmapValues contains no finite values to plot.')
end
observedMaximum = max(finiteValues);
upperColorLimit = max(colorStep, ...
    colorStep .* ceil(observedMaximum ./ colorStep));

% Keep the colorbar readable when the rounded range contains many levels.
allColorTicks = 0:colorStep:upperColorLimit;
if numel(allColorTicks) > 7
    tickStride = ceil((numel(allColorTicks) - 1) ./ 5);
    colorTicks = allColorTicks(1:tickStride:end);
    if colorTicks(end) < upperColorLimit
        colorTicks(end + 1) = upperColorLimit;
    end
else
    colorTicks = allColorTicks;
end

% Use log10 coordinates so every image cell correctly represents the
% logarithmically spaced parameter grid.
logEtaValues = log10(etaValues);
logSigmaValues = log10(sigmaValues);
etaTickValues = [1e-3, 1e-1, 10];
sigmaTickValues = [1e-3, 1e-1, 10];

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
figureFontSize = 60;
colorbarFontSize = 42;
axesLineWidth = 1.5;

% The only output figure shows the complete requested parameter range.
figureHandle = figure(figureNumber);
clf(figureHandle)
% The window is wide enough to keep the whole colorbar label on the canvas.
set(figureHandle, 'Color', 'w', 'Position', [180 100 1090 850])
axesHandle = axes(figureHandle);
imagesc(axesHandle, logEtaValues, logSigmaValues, heatmapValues)
set(axesHandle, 'YDir', 'normal', 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth, 'TickDir', 'out', 'Layer', 'top')
xlabel(axesHandle, '$\eta$', 'Interpreter', 'latex')
ylabel(axesHandle, '$\sigma$', 'Interpreter', 'latex')
xlim(axesHandle, [min(logEtaValues), max(logEtaValues)])
ylim(axesHandle, [min(logSigmaValues), max(logSigmaValues)])
xticks(axesHandle, log10(etaTickValues))
yticks(axesHandle, log10(sigmaTickValues))
xticklabels(axesHandle, {'$10^{-3}$', '$10^{-1}$', '$10^{1}$'})
yticklabels(axesHandle, {'$10^{-3}$', '$10^{-1}$', '$10^{1}$'})
clim(axesHandle, [0, upperColorLimit])
colormap(axesHandle, parula(256))
pbaspect(axesHandle, [1 1 1])
box(axesHandle, 'on')

% Label the color scale explicitly and retain the established formatting.
colorbarHandle = colorbar(axesHandle);
colorbarHandle.FontSize = colorbarFontSize;
colorbarHandle.LineWidth = axesLineWidth;
colorbarHandle.Ticks = colorTicks;
% defaultAxesTickLabelInterpreter does not reach colorbars, so set it here.
colorbarHandle.TickLabelInterpreter = 'latex';
colorbarHandle.Label.String = colorbarLabel;
colorbarHandle.Label.Interpreter = 'latex';
colorbarHandle.Label.FontSize = colorbarFontSize;
end
