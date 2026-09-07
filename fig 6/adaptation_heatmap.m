% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Calculate the surviving fraction across adaptation and noise strengths.

% Generate one reproducible set of species-specific resource-use traits.
% Every parameter pair uses this community and the same initial conditions.
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

% Fixed diffusion rate and logarithmically spaced heatmap coordinates.
beta = 1;
etaValues = logspace(-3, 1, 100);
sigmaValues = logspace(-3, 1, 100);

% Simulation settings retained from the stochastic adaptive model.
tMax = 5000;
explicitSafetyFactor = 0.5;
maximumDt = min(0.001, explicitSafetyFactor ./ ...
    ((model.S + 1) .* beta));

nEta = numel(etaValues);
nSigma = numel(sigmaValues);
nParameterPairs = nEta .* nSigma;

% Store only final metrics and timestep diagnostics for each parameter pair.
survivingFractionFlat = zeros(nParameterPairs, 1);
minimumDtUsedFlat = zeros(nParameterPairs, 1);
maximumDtUsedFlat = zeros(nParameterPairs, 1);
numStepsFlat = zeros(nParameterPairs, 1);

% Simulate every independent parameter pair concurrently. Each pair gets a
% distinct reproducible noise substream, while species traits remain fixed.
parfor pairIndex = 1:nParameterPairs
    [sigmaIndex, etaIndex] = ind2sub([nSigma, nEta], pairIndex);
    eta = etaValues(etaIndex);
    sigma = sigmaValues(sigmaIndex);

    pairNoiseStream = RandStream('Threefry', 'Seed', masterSeed);
    pairNoiseStream.Substream = pairIndex + 1;

    [survivingFractionFlat(pairIndex), ...
        minimumDtUsedFlat(pairIndex), ...
        maximumDtUsedFlat(pairIndex), numStepsFlat(pairIndex)] = ...
        simulateFinalSurvivingFraction( ...
        eta, sigma, beta, tMax, maximumDt, ...
        explicitSafetyFactor, model, pairNoiseStream);
end

% Rows correspond to sigma and columns correspond to eta for imagesc.
survivingFraction = reshape( ...
    survivingFractionFlat, [nSigma, nEta]);

% Plot one square panel over the complete logarithmic parameter range. The
% color range ends at the observed maximum rounded up to 0.05.
plotAdaptationHeatmapFullRange(etaValues, sigmaValues, ...
    survivingFraction, 'Fraction of surviving species', 1, 0.05)

% Report global integration diagnostics after all workers finish.
fprintf('Completed %d eta-sigma parameter pairs.\n', nParameterPairs)
fprintf('Overall timestep range: [%.6g, %.6g].\n', ...
    min(minimumDtUsedFlat), max(maximumDtUsedFlat))
fprintf('Total adaptive Euler steps: %d.\n', sum(numStepsFlat))

function [survivingFraction, minimumDtUsed, maximumDtUsed, numSteps] = ...
    simulateFinalSurvivingFraction( ...
    eta, sigma, beta, tMax, maximumDt, ...
    explicitSafetyFactor, model, noiseStream)

% Simulate one stochastic parameter pair and retain only final richness.
N = model.N0 .* ones(1, model.S);
R = model.R0;
T = model.T0 .* ones(1, model.S);
theta = zeros(1, model.S);

currentTime = 0;
minimumDtUsed = inf;
maximumDtUsed = 0;
numSteps = 0;

while currentTime < tMax
    % Calculate each species' deterministic per-capita growth signal.
    resourceUse = theta .* T + (1 - theta) .* R;
    perCapitaGrowth = model.alpha .* resourceUse - model.k;

    % Calculate the deterministic consumer-resource rates.
    dN = perCapitaGrowth .* N;
    dT = beta .* (R - T) - ...
        model.alpha .* model.gamma .* theta .* T .* N;
    dR = model.Rstar - R - beta .* ...
        (model.S .* R - sum(T)) - model.gamma .* R .* ...
        sum(model.alpha .* (1 - theta) .* N);

    % Bound the timestep using current resource, demographic, adaptation,
    % and stochastic rates. Extinct species have frozen N_i and theta_i.
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

    noiseRate = sigma .^ 2;
    maximumCurrentRate = max([localRemovalRate, sharedRemovalRate, ...
        demographicRate, adaptationRate, noiseRate]);

    if ~isfinite(maximumCurrentRate) || maximumCurrentRate <= 0
        error(['An invalid adaptive rate occurred at t = %.6g for ' ...
            'eta = %.6g and sigma = %.6g.'], currentTime, eta, sigma)
    end

    dt = min(maximumDt, ...
        explicitSafetyFactor ./ maximumCurrentRate);
    dt = min(dt, tMax - currentTime);

    if ~isfinite(dt) || dt <= 0
        error(['An invalid timestep occurred at t = %.6g for ' ...
            'eta = %.6g and sigma = %.6g.'], currentTime, eta, sigma)
    end

    % Advance abundances with Euler-Maruyama. Skip the random-number draw
    % along the deterministic sigma=0 boundary.
    if sigma > 0
        brownianIncrement = sqrt(dt) .* ...
            randn(noiseStream, 1, model.S);
        NNew = N + dt .* dN + sigma .* N .* brownianIncrement;
    else
        NNew = N + dt .* dN;
    end

    % Resources remain deterministic and are advanced explicitly.
    TNew = T + dt .* dT;
    RNew = R + dt .* dR;

    % Force populations below the extinction threshold permanently to zero.
    NNew(NNew < model.extinctionThreshold) = 0;

    % Advance theta semi-implicitly using the deterministic growth signal.
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

    % Extinct species retain their final theta, and roundoff is bounded.
    thetaNew(NNew == 0) = theta(NNew == 0);
    thetaNew = min(max(thetaNew, 0), 1);

    % Commit the completed stochastic adaptive step.
    N = NNew;
    T = TNew;
    R = RNew;
    theta = thetaNew;
    currentTime = currentTime + dt;
    numSteps = numSteps + 1;
    minimumDtUsed = min(minimumDtUsed, dt);
    maximumDtUsed = max(maximumDtUsed, dt);

    % Stop instead of retaining a numerically invalid trajectory.
    if any(~isfinite(N), 'all') || ~isfinite(R) || ...
            any(~isfinite(T), 'all') || any(~isfinite(theta), 'all')
        error(['A nonfinite state occurred at t = %.6g for eta = %.6g ' ...
            'and sigma = %.6g.'], currentTime, eta, sigma)
    end
    if R < 0 || any(T < 0, 'all')
        error(['A negative resource state occurred at t = %.6g for ' ...
            'eta = %.6g and sigma = %.6g.'], currentTime, eta, sigma)
    end
end

survivingFraction = nnz(N > model.extinctionThreshold) ./ model.S;
end


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
