% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Calculate replicate-averaged survival for complementary adaptation-rate
% and multiplicative-noise sweeps.

% Use separate reproducible seeds for community traits and stochastic paths.
masterSeed = 1;
noiseSeed = 2;
traitStream = RandStream('Threefry', 'Seed', masterSeed);

% Consumer-resource model parameters retained from the averaged sweep.
model.S = 50;
model.k = 1;
model.gamma = 0.01;
model.Rstar = 5;
model.N0 = 1;
model.R0 = model.Rstar;
model.T0 = model.Rstar;
model.extinctionThreshold = 1e-3;

% Fixed diffusion rate and the requested one-dimensional parameter sweeps.
beta = 1;
nSweepPoints = 100;
etaSweepValues = linspace(0, 1, nSweepPoints);
sigmaForEtaSweep = [1e-2, 1e-1, 0.5, 1e0];
sigmaSweepValues = linspace(0, 1, nSweepPoints);
etaForSigmaSweep = [1e-2, 1e-1, 1e0, 1e1];

% Sample 20 independent communities once and reuse the same communities at
% every parameter value. This paired design reduces community-sampling noise
% when comparing curves, while each trajectory has an independent noise path.
nReplicates = 20;
alphaByReplicate = zeros(nReplicates, model.S);
for replicateIndex = 1:nReplicates
    traitStream.Substream = replicateIndex;
    alphaByReplicate(replicateIndex, :) = ...
        1 + 4 .* rand(traitStream, 1, model.S);
end

% Simulation settings retained from the stochastic adaptive model.
tMax = 5000;
explicitSafetyFactor = 0.5;
maximumDt = min(0.001, explicitSafetyFactor ./ ...
    ((model.S + 1) .* beta));

% Sweep 1: vary eta from 0 to 1 for each requested fixed sigma value.
nEtaSweep = numel(etaSweepValues);
nFixedSigma = numel(sigmaForEtaSweep);
nEtaSweepSimulations = nEtaSweep .* nFixedSigma .* nReplicates;

etaSweepFractionFlat = zeros(nEtaSweepSimulations, 1);
etaSweepMinimumDtFlat = zeros(nEtaSweepSimulations, 1);
etaSweepMaximumDtFlat = zeros(nEtaSweepSimulations, 1);
etaSweepNumStepsFlat = zeros(nEtaSweepSimulations, 1);

% Run every eta-sigma-community combination concurrently. The substream is
% determined only by the array index, so results do not depend on scheduling.
parfor simulationIndex = 1:nEtaSweepSimulations
    [etaIndex, sigmaCaseIndex, replicateIndex] = ind2sub( ...
        [nEtaSweep, nFixedSigma, nReplicates], simulationIndex);
    eta = etaSweepValues(etaIndex);
    sigma = sigmaForEtaSweep(sigmaCaseIndex);

    replicateModel = model;
    replicateModel.alpha = alphaByReplicate(replicateIndex, :);
    noiseStream = RandStream('Threefry', 'Seed', noiseSeed);
    noiseStream.Substream = simulationIndex;

    [etaSweepFractionFlat(simulationIndex), ...
        etaSweepMinimumDtFlat(simulationIndex), ...
        etaSweepMaximumDtFlat(simulationIndex), ...
        etaSweepNumStepsFlat(simulationIndex)] = ...
        simulateFinalSurvivingFraction( ...
        eta, sigma, beta, tMax, maximumDt, ...
        explicitSafetyFactor, replicateModel, noiseStream);
end

% Recover eta-by-sigma-by-replicate arrays and calculate ensemble statistics.
etaSweepFractionByReplicate = reshape(etaSweepFractionFlat, ...
    [nEtaSweep, nFixedSigma, nReplicates]);
etaSweepMean = mean(etaSweepFractionByReplicate, 3, 'omitnan');
etaSweepStd = std(etaSweepFractionByReplicate, 0, 3, 'omitnan');

% Sweep 2: vary sigma from 0 to 1 for each requested fixed eta value.
nSigmaSweep = numel(sigmaSweepValues);
nFixedEta = numel(etaForSigmaSweep);
nSigmaSweepSimulations = nSigmaSweep .* nFixedEta .* nReplicates;

sigmaSweepFractionFlat = zeros(nSigmaSweepSimulations, 1);
sigmaSweepMinimumDtFlat = zeros(nSigmaSweepSimulations, 1);
sigmaSweepMaximumDtFlat = zeros(nSigmaSweepSimulations, 1);
sigmaSweepNumStepsFlat = zeros(nSigmaSweepSimulations, 1);

parfor simulationIndex = 1:nSigmaSweepSimulations
    [sigmaIndex, etaCaseIndex, replicateIndex] = ind2sub( ...
        [nSigmaSweep, nFixedEta, nReplicates], simulationIndex);
    sigma = sigmaSweepValues(sigmaIndex);
    eta = etaForSigmaSweep(etaCaseIndex);

    replicateModel = model;
    replicateModel.alpha = alphaByReplicate(replicateIndex, :);
    noiseStream = RandStream('Threefry', 'Seed', noiseSeed);
    noiseStream.Substream = nEtaSweepSimulations + simulationIndex;

    [sigmaSweepFractionFlat(simulationIndex), ...
        sigmaSweepMinimumDtFlat(simulationIndex), ...
        sigmaSweepMaximumDtFlat(simulationIndex), ...
        sigmaSweepNumStepsFlat(simulationIndex)] = ...
        simulateFinalSurvivingFraction( ...
        eta, sigma, beta, tMax, maximumDt, ...
        explicitSafetyFactor, replicateModel, noiseStream);
end

% Recover sigma-by-eta-by-replicate arrays and calculate ensemble statistics.
sigmaSweepFractionByReplicate = reshape(sigmaSweepFractionFlat, ...
    [nSigmaSweep, nFixedEta, nReplicates]);
sigmaSweepMean = mean(sigmaSweepFractionByReplicate, 3, 'omitnan');
sigmaSweepStd = std(sigmaSweepFractionByReplicate, 0, 3, 'omitnan');

% Record one additional stochastic trajectory for examining abundance and
% theta dynamics. Use community 1 and a noise substream not used by either
% parameter sweep.
exampleEta = 0.25;
exampleSigma = 0.1;
exampleOutputInterval = 1;
exampleTime = (0:exampleOutputInterval:tMax).';
if exampleTime(end) < tMax
    exampleTime(end + 1, 1) = tMax;
end

exampleModel = model;
exampleModel.alpha = alphaByReplicate(1, :);
exampleNoiseStream = RandStream('Threefry', 'Seed', noiseSeed);
exampleNoiseStream.Substream = ...
    nEtaSweepSimulations + nSigmaSweepSimulations + 1;

[exampleAbundanceHistory, exampleThetaHistory, ...
    exampleFinalAbundances, exampleMinimumDtUsed, ...
    exampleMaximumDtUsed, exampleNumSteps] = ...
    simulateRecordedTrajectory( ...
    exampleEta, exampleSigma, beta, exampleTime, maximumDt, ...
    explicitSafetyFactor, exampleModel, exampleNoiseStream);
exampleSurvivingSpecies = ...
    exampleFinalAbundances > model.extinctionThreshold;

% Apply the shared presentation-scale figure formatting.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
figureFontSize = 60;
legendFontSize = 36;
markerSize = 8;
curveLineWidth = 3;
axesLineWidth = 1.5;
curveColors = lines(max(nFixedSigma, nFixedEta));

% Figure 1: place the eta and sigma sweeps in two side-by-side panels.
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

% Figure 2: one stored stochastic abundance and theta trajectory.
plotExampleTrajectory(2, exampleTime, exampleAbundanceHistory, ...
    exampleThetaHistory, exampleSurvivingSpecies, exampleEta, ...
    exampleSigma, model.extinctionThreshold, figureFontSize, ...
    curveLineWidth, axesLineWidth)

% Report the computational scale and integration diagnostics.
allMinimumDt = [etaSweepMinimumDtFlat; sigmaSweepMinimumDtFlat];
allMaximumDt = [etaSweepMaximumDtFlat; sigmaSweepMaximumDtFlat];
allNumSteps = [etaSweepNumStepsFlat; sigmaSweepNumStepsFlat];
totalSimulations = nEtaSweepSimulations + nSigmaSweepSimulations;

fprintf(['Completed %d stochastic trajectories across both sweeps ' ...
    'using %d communities.\n'], totalSimulations, nReplicates)
fprintf('Eta-sweep trajectories: %d.\n', nEtaSweepSimulations)
fprintf('Sigma-sweep trajectories: %d.\n', nSigmaSweepSimulations)
fprintf('Overall timestep range: [%.6g, %.6g].\n', ...
    min(allMinimumDt), max(allMaximumDt))
fprintf('Total adaptive Euler steps: %.0f.\n', sum(allNumSteps))
fprintf(['Stored example trajectory: eta = %.6g, sigma = %.6g, ' ...
    'community = 1.\n'], exampleEta, exampleSigma)
fprintf('Example timestep range: [%.6g, %.6g].\n', ...
    exampleMinimumDtUsed, exampleMaximumDtUsed)
fprintf('Example adaptive Euler steps: %.0f.\n', exampleNumSteps)
fprintf('Example surviving species at tMax: %d of %d.\n', ...
    nnz(exampleSurvivingSpecies), model.S)

function plotOverlaySweep(axesHandle, xValues, ...
    meanValues, standardDeviations, fixedValues, fixedParameterName, ...
    xAxisLabel, xTicks, curveColors, figureFontSize, legendFontSize, ...
    markerSize, curveLineWidth, axesLineWidth)

% Plot replicate means with one-standard-deviation error bars in one panel.
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

function figureHandle = plotExampleTrajectory(figureNumber, timeValues, ...
    abundanceHistory, thetaHistory, survivingSpecies, eta, sigma, ...
    extinctionThreshold, figureFontSize, curveLineWidth, axesLineWidth)

% Plot the stored example as two square side-by-side panels, matching the
% panel geometry and spacing used in optimal_beta_final.m.
nSpecies = size(abundanceHistory, 2);
if ~isequal(size(thetaHistory), size(abundanceHistory)) || ...
        size(abundanceHistory, 1) ~= numel(timeValues) || ...
        numel(survivingSpecies) ~= nSpecies
    error('The stored example trajectory has inconsistent dimensions.')
end

speciesColors = lines(nSpecies);
tMinimum = min(timeValues);
tMaximum = max(timeValues);
panelTitle = sprintf('$\\eta = %.3g,\\; \\sigma = %.3g$', eta, sigma);

positiveAbundances = abundanceHistory(abundanceHistory > 0);
if isempty(positiveAbundances)
    error('The stored example contains no positive abundances to plot.')
end
abundanceUpperLimit = max(positiveAbundances);
candidateAbundanceTicks = [1e-2, 1e0, 1e2];
abundanceTicks = candidateAbundanceTicks( ...
    candidateAbundanceTicks <= abundanceUpperLimit);

figureHandle = figure(figureNumber);
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

% Right panel: theta dynamics only for species surviving at tMax.
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
end

function [survivingFraction, minimumDtUsed, maximumDtUsed, numSteps] = ...
    simulateFinalSurvivingFraction( ...
    eta, sigma, beta, tMax, maximumDt, ...
    explicitSafetyFactor, model, noiseStream)

% Simulate one stochastic trajectory and retain only final species richness.
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

    % Advance abundances with Euler-Maruyama. Skip random draws at sigma=0.
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

% Normalize final richness by the initial number of species.
survivingFraction = nnz(N > model.extinctionThreshold) ./ model.S;
end

function [abundanceHistory, thetaHistory, N, minimumDtUsed, ...
    maximumDtUsed, numSteps] = simulateRecordedTrajectory( ...
    eta, sigma, beta, sampleTimes, maximumDt, ...
    explicitSafetyFactor, model, noiseStream)

% Simulate one stochastic trajectory while retaining abundance and theta at
% fixed output times. The numerical update matches the parameter sweeps.
sampleTimes = sampleTimes(:);
if isempty(sampleTimes) || sampleTimes(1) ~= 0 || ...
        any(diff(sampleTimes) <= 0)
    error('sampleTimes must start at zero and increase strictly.')
end

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
    % Calculate the current deterministic per-capita growth signal.
    resourceUse = theta .* T + (1 - theta) .* R;
    perCapitaGrowth = model.alpha .* resourceUse - model.k;

    % Calculate deterministic abundance and resource rates.
    dN = perCapitaGrowth .* N;
    dT = beta .* (R - T) - ...
        model.alpha .* model.gamma .* theta .* T .* N;
    dR = model.Rstar - R - beta .* ...
        (model.S .* R - sum(T)) - model.gamma .* R .* ...
        sum(model.alpha .* (1 - theta) .* N);

    % Apply the same adaptive timestep rule used in the two sweeps.
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
        error(['An invalid recorded-trajectory rate occurred at ' ...
            't = %.6g.'], currentTime)
    end

    dt = min(maximumDt, ...
        explicitSafetyFactor ./ maximumCurrentRate);
    dt = min(dt, tMax - currentTime);

    % Land exactly on the next requested output time.
    if sampleIndex <= nSamples
        dt = min(dt, sampleTimes(sampleIndex) - currentTime);
    end
    if ~isfinite(dt) || dt <= 0
        error(['An invalid recorded-trajectory timestep occurred at ' ...
            't = %.6g.'], currentTime)
    end

    % Advance abundances with Euler-Maruyama.
    if sigma > 0
        brownianIncrement = sqrt(dt) .* ...
            randn(noiseStream, 1, model.S);
        NNew = N + dt .* dN + sigma .* N .* brownianIncrement;
    else
        NNew = N + dt .* dN;
    end

    % Advance deterministic resources explicitly.
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
    thetaNew(NNew == 0) = theta(NNew == 0);
    thetaNew = min(max(thetaNew, 0), 1);

    % Commit the completed step and update diagnostics.
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
        error(['A nonfinite recorded state occurred at t = %.6g.'], ...
            currentTime)
    end
    if R < 0 || any(T < 0, 'all')
        error(['A negative recorded resource state occurred at ' ...
            't = %.6g.'], currentTime)
    end

    % Record the state when the next requested physical time is reached.
    if sampleIndex <= nSamples && ...
            currentTime >= sampleTimes(sampleIndex) - timeTolerance
        abundanceHistory(sampleIndex, :) = N;
        thetaHistory(sampleIndex, :) = theta;
        sampleIndex = sampleIndex + 1;
    end
end

if sampleIndex <= nSamples
    error('The recorded trajectory did not reach every requested time.')
end
end
