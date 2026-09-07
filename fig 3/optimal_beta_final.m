% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all

% Model and simulation settings shared by every parameter pair.
masterSeed = 1;
sigma = 1;
nTrajectories = 20;
observationTime = 2000;

model.S = 50;
model.k = 1;
model.gamma = 0.01;
model.Rstar = 5;
model.N0 = 1;
model.R0 = model.Rstar;
model.T0 = model.Rstar;
model.extinctionThreshold = 1e-3;

traitStream = RandStream('Threefry', 'Seed', masterSeed);
traitStream.Substream = 1;
model.alpha = 1 + 4 * rand(traitStream, 1, model.S);

thetaValues = [0 0.3 0.5 0.8 1];
betaValues = logspace(-3, 3, 30);

% The base cap resolves exchange; each step is reduced for current consumption.
explicitSafetyFactor = 0.5;
maximumDtValues = min(0.001, explicitSafetyFactor ./ ...
    ((model.S + 1) .* betaValues));

% Use a fixed deterministic burn-in without testing for steady state.
burnTime = 3000;

nBeta = numel(betaValues);
nTheta = numel(thetaValues);
nParameterPairs = nBeta * nTheta;

% Deterministic burn-in is calculated once for each parameter pair.
burnN = zeros(nParameterPairs, model.S);
burnR = zeros(nParameterPairs, 1);
burnT = zeros(nParameterPairs, model.S);
burnRichness = zeros(nParameterPairs, 1);
burnNonfiniteTime = nan(nParameterPairs, 1);
burnNegativeResourceTime = nan(nParameterPairs, 1);

parfor pair = 1:nParameterPairs
    [p, g] = ind2sub([nBeta, nTheta], pair);
    [Nburn, Rburn, Tburn, nonfiniteTime, ...
        negativeResourceTime] = simulateBurnIn( ...
        betaValues(p), thetaValues(g), maximumDtValues(p), ...
        explicitSafetyFactor, burnTime, model);

    burnN(pair, :) = Nburn;
    burnR(pair) = Rburn;
    burnT(pair, :) = Tburn;
    burnRichness(pair) = nnz(Nburn);
    burnNonfiniteTime(pair) = nonfiniteTime;
    burnNegativeResourceTime(pair) = negativeResourceTime;
end

nonfiniteBurnPairs = find(~isnan(burnNonfiniteTime));
if ~isempty(nonfiniteBurnPairs)
    [p, g] = ind2sub([nBeta, nTheta], nonfiniteBurnPairs);
    nonfiniteBurnReport = table( ...
        reshape(betaValues(p), [], 1), reshape(thetaValues(g), [], 1), ...
        burnNonfiniteTime(nonfiniteBurnPairs), ...
        'VariableNames', {'beta', 'theta', 'firstNonfiniteTime'});
    disp(nonfiniteBurnReport)
end

negativeBurnPairs = find(~isnan(burnNegativeResourceTime));
if ~isempty(negativeBurnPairs)
    [p, g] = ind2sub([nBeta, nTheta], negativeBurnPairs);
    negativeBurnResourceReport = table( ...
        reshape(betaValues(p), [], 1), reshape(thetaValues(g), [], 1), ...
        burnNegativeResourceTime(negativeBurnPairs), ...
        'VariableNames', {'beta', 'theta', 'firstNegativeResourceTime'});
    disp(negativeBurnResourceReport)
end

if ~isempty(nonfiniteBurnPairs) || ~isempty(negativeBurnPairs)
    error('Invalid resource states occurred during burn-in; observations were not run.')
end

% Simulate all trajectories together for each parameter pair.
survivingFractionS0Flat = nan(nParameterPairs, nTrajectories);
relativeRichnessChangeFlat = nan(nParameterPairs, nTrajectories);
nonfiniteTimeFlat = nan(nParameterPairs, nTrajectories);
negativeResourceTimeFlat = nan(nParameterPairs, nTrajectories);
observationStreams = parallel.pool.Constant( ...
    RandStream('Threefry', 'Seed', masterSeed));

parfor pair = 1:nParameterPairs
    [p, g] = ind2sub([nBeta, nTheta], pair);
    Sburn = burnRichness(pair);

    if Sburn == 0
        % No species remain to simulate after burn-in.
        survivingFractionS0Flat(pair, :) = 0;
    else
        % One fixed substream generates 20 independent trajectory rows.
        stream = observationStreams.Value;
        stream.Substream = pair + 1;
        [finalRichness, firstNonfiniteTimes, ...
            firstNegativeResourceTimes] = ...
            simulateObservationEnsemble(burnN(pair, :), burnR(pair), ...
            burnT(pair, :), betaValues(p), thetaValues(g), ...
            sigma, observationTime, maximumDtValues(p), ...
            explicitSafetyFactor, nTrajectories, model, stream);

        survivingFractionS0Flat(pair, :) = (finalRichness / model.S).';
        relativeRichnessChangeFlat(pair, :) = ...
            ((finalRichness - Sburn) / Sburn).';
        nonfiniteTimeFlat(pair, :) = firstNonfiniteTimes.';
        negativeResourceTimeFlat(pair, :) = firstNegativeResourceTimes.';
    end
end
delete(observationStreams)

survivingFractionS0 = reshape(survivingFractionS0Flat, ...
    [nBeta, nTheta, nTrajectories]);
relativeRichnessChange = reshape(relativeRichnessChangeFlat, ...
    [nBeta, nTheta, nTrajectories]);

% Report invalid trajectories explicitly and exclude them from all metrics.
nonfiniteIndices = find(~isnan(nonfiniteTimeFlat));
if ~isempty(nonfiniteIndices)
    [pair, trajectory] = ind2sub( ...
        [nParameterPairs, nTrajectories], nonfiniteIndices);
    [p, g] = ind2sub([nBeta, nTheta], pair);
    nonfiniteTrajectoryReport = table( ...
        reshape(betaValues(p), [], 1), reshape(thetaValues(g), [], 1), ...
        trajectory, nonfiniteTimeFlat(nonfiniteIndices), ...
        'VariableNames', {'beta', 'theta', 'trajectory', ...
        'firstNonfiniteTime'});
    warning('%d nonfinite trajectories were excluded from the results.', ...
        numel(nonfiniteIndices))
    disp(nonfiniteTrajectoryReport)
end

negativeResourceIndices = find(~isnan(negativeResourceTimeFlat));
if ~isempty(negativeResourceIndices)
    [pair, trajectory] = ind2sub( ...
        [nParameterPairs, nTrajectories], negativeResourceIndices);
    [p, g] = ind2sub([nBeta, nTheta], pair);
    negativeResourceTrajectoryReport = table( ...
        reshape(betaValues(p), [], 1), reshape(thetaValues(g), [], 1), ...
        trajectory, negativeResourceTimeFlat(negativeResourceIndices), ...
        'VariableNames', {'beta', 'theta', 'trajectory', ...
        'firstNegativeResourceTime'});
    warning('%d negative-resource trajectories were excluded from the results.', ...
        numel(negativeResourceIndices))
    disp(negativeResourceTrajectoryReport)
end

if isempty(nonfiniteIndices) && isempty(negativeResourceIndices)
    fprintf('Resource-state check: no invalid states detected.\n')
end

%% Export trajectory-level results and selected abundance examples

% Write all generated files to one folder in the job's working directory.
exportDirectory = fullfile(pwd, 'optimal_beta_final_outputs');
if ~isfolder(exportDirectory)
    mkdir(exportDirectory)
end

survivalFile = fullfile( ...
    exportDirectory, 'surviving_fraction_trajectories.xlsx');
richnessChangeFile = fullfile( ...
    exportDirectory, 'richness_change_trajectories.xlsx');
exampleAbundanceFile = fullfile( ...
    exportDirectory, 'example_abundance_trajectories.xlsx');

% Each theta sheet contains beta and all trajectory-level panel values.
trajectoryColumnNames = [{'beta'}, ...
    cellstr(compose("trajectory_%02d", 1:nTrajectories))];
for g = 1:nTheta
    thetaSheet = sprintf('theta=%.2g', thetaValues(g));
    if g == 1
        sheetWriteMode = 'replacefile';
    else
        sheetWriteMode = 'overwritesheet';
    end

    survivalValues = reshape( ...
        survivingFractionS0(:, g, :), [nBeta, nTrajectories]);
    survivalTable = array2table( ...
        [betaValues(:), survivalValues], ...
        'VariableNames', trajectoryColumnNames);
    writetable(survivalTable, survivalFile, ...
        'Sheet', thetaSheet, 'WriteMode', sheetWriteMode)

    richnessChangeValues = reshape( ...
        relativeRichnessChange(:, g, :), [nBeta, nTrajectories]);
    richnessChangeTable = array2table( ...
        [betaValues(:), richnessChangeValues], ...
        'VariableNames', trajectoryColumnNames);
    writetable(richnessChangeTable, richnessChangeFile, ...
        'Sheet', thetaSheet, 'WriteMode', sheetWriteMode)
end

% Select one reproducible beta in [1, 10] for every theta.
eligibleBetaIndices = find(betaValues >= 1 & betaValues <= 10);
if isempty(eligibleBetaIndices)
    error('No simulated beta values lie between 1 and 10.')
end
exampleSelectionStream = RandStream('Threefry', 'Seed', masterSeed);
exampleSelectionStream.Substream = nParameterPairs + 2;
selectionIndices = randi(exampleSelectionStream, ...
    numel(eligibleBetaIndices), nTheta, 1);
exampleBetaIndices = eligibleBetaIndices(selectionIndices);
exampleBetaValues = betaValues(exampleBetaIndices);
exampleOutputInterval = 1;

exampleSelectionReport = table(thetaValues(:), exampleBetaValues(:), ...
    'VariableNames', {'theta', 'beta'});
disp('Selected parameter pairs for example abundance trajectories:')
disp(exampleSelectionReport)

speciesColumnNames = cellstr(compose("species_%02d", 1:model.S));
for g = 1:nTheta
    p = exampleBetaIndices(g);
    beta = betaValues(p);
    theta = thetaValues(g);

    % Record only this selected trajectory, not every trajectory in the sweep.
    N = model.N0 * ones(1, model.S);
    R = model.R0;
    T = model.T0 * ones(1, model.S);
    [burnTimes, burnAbundances, N, R, T] = simulateRecordedPhase( ...
        N, R, T, beta, theta, burnTime, exampleOutputInterval, ...
        maximumDtValues(p), explicitSafetyFactor, 0, model, []);

    exampleStream = RandStream('Threefry', 'Seed', masterSeed);
    exampleStream.Substream = nParameterPairs + 2 + g;
    [observationTimes, observationAbundances] = simulateRecordedPhase( ...
        N, R, T, beta, theta, observationTime, exampleOutputInterval, ...
        maximumDtValues(p), explicitSafetyFactor, sigma, model, ...
        exampleStream);

    % Join the phases once at the shared burn-in endpoint.
    exampleTimes = [burnTimes; burnTime + observationTimes(2:end)];
    exampleAbundances = [burnAbundances; observationAbundances(2:end, :)];
    examplePhases = [repmat("burn_in", numel(burnTimes), 1); ...
        repmat("observation", numel(observationTimes) - 1, 1)];
    exampleTable = table(exampleTimes, examplePhases, ...
        'VariableNames', {'time', 'phase'});
    exampleTable = [exampleTable, array2table(exampleAbundances, ...
        'VariableNames', speciesColumnNames)];

    exampleSheet = sprintf( ...
        'theta=%.2g, beta=%.3g', theta, beta);
    if g == 1
        sheetWriteMode = 'replacefile';
    else
        sheetWriteMode = 'overwritesheet';
    end
    writetable(exampleTable, exampleAbundanceFile, ...
        'Sheet', exampleSheet, 'WriteMode', sheetWriteMode)
end

% Register the workbooks with a batch job so nonshared clusters return them.
excelFileStoreKeys = ["surviving_fraction_trajectories.xlsx"; ...
    "richness_change_trajectories.xlsx"; ...
    "example_abundance_trajectories.xlsx"];
excelFiles = {survivalFile; richnessChangeFile; exampleAbundanceFile};
jobFileStore = getCurrentFileStore;
if ~isempty(jobFileStore)
    for fileIndex = 1:numel(excelFiles)
        copyFileToStore(jobFileStore, excelFiles{fileIndex}, ...
            excelFileStoreKeys(fileIndex))
    end
    fprintf('Excel files copied to the batch job FileStore.\n')
end
fprintf('Excel files written to: %s\n', exportDirectory)

% Export burn-in surviving fraction

burnDataExportDirectory = fullfile(pwd, 'optimal_beta_final_outputs');
if ~isfolder(burnDataExportDirectory)
    mkdir(burnDataExportDirectory)
end
burnSurvivalFile = fullfile( ...
    burnDataExportDirectory, 'burn_surviving_fraction.xlsx');
burnSurvivingFractionForExport = reshape( ...
    burnRichness / model.S, [nBeta, nTheta]);

for g = 1:nTheta
    thetaSheet = sprintf('theta=%.2g', thetaValues(g));
    burnSurvivalTable = table(betaValues(:), ...
        burnSurvivingFractionForExport(:, g), ...
        'VariableNames', {'beta', 'S_burn_over_S0'});
    if g == 1
        sheetWriteMode = 'replacefile';
    else
        sheetWriteMode = 'overwritesheet';
    end
    writetable(burnSurvivalTable, burnSurvivalFile, ...
        'Sheet', thetaSheet, 'WriteMode', sheetWriteMode)
end

% Make the new workbook retrievable when this section runs in a batch job.
burnFileStore = getCurrentFileStore;
if ~isempty(burnFileStore)
    copyFileToStore(burnFileStore, burnSurvivalFile, ...
        "burn_surviving_fraction.xlsx")
end
fprintf('Figure 2 data written to: %s\n', burnSurvivalFile)

%% Plot
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
colors = lines(nTheta);
legendLabels = compose('$\\theta = %.2g$', thetaValues);
richnessMarkerSize = 15;
richnessLineWidth = 3;
richnessFontSize = 60;
leftPanelMinimum = min( ...
    meanSurvivingFractionS0 - stdSurvivingFractionS0, ...
    [], 'all', 'omitnan');
rightPanelMinimum = min( ...
    meanRelativeRichnessChange - stdRelativeRichnessChange, ...
    [], 'all', 'omitnan');
if ~isfinite(leftPanelMinimum) || leftPanelMinimum >= 1
    leftPanelMinimum = 0.9;
end
if ~isfinite(rightPanelMinimum) || rightPanelMinimum >= 0
    rightPanelMinimum = -0.1;
end

% Figure 1: richness metrics shown as two square panels.
figure(1)
clf
set(gcf, 'Color', 'w', 'Position', [50 80 1900 900])
tiledlayout(1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: fraction surviving relative to the original S0 = 50 species.
nexttile
hold on
richnessHandles = gobjects(nTheta, 1);
for g = 1:nTheta
    richnessHandles(g) = errorbar( ...
        betaValues, meanSurvivingFractionS0(:, g), ...
        stdSurvivingFractionS0(:, g), ...
        'Color', colors(g, :), 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerSize', richnessMarkerSize, ...
        'MarkerFaceColor', colors(g, :), ...
        'MarkerEdgeColor', colors(g, :), ...
        'LineWidth', richnessLineWidth, 'CapSize', 10);
end
formatRichnessPanel( ...
    'Fraction of surviving species', richnessFontSize, betaValues, ...
    [leftPanelMinimum 1])
legend(richnessHandles, legendLabels, ...
    'Location', 'best', 'Orientation', 'vertical', 'FontSize', 42)

% Panel 2: fractional richness change relative to the burn-in state.
nexttile
hold on
for g = 1:nTheta
    errorbar(betaValues, meanRelativeRichnessChange(:, g), ...
        stdRelativeRichnessChange(:, g), ...
        'Color', colors(g, :), 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerSize', richnessMarkerSize, ...
        'MarkerFaceColor', colors(g, :), ...
        'MarkerEdgeColor', colors(g, :), ...
        'LineWidth', richnessLineWidth, 'CapSize', 10);
end
formatRichnessPanel( ...
    'Fractional change of species from steady state', ...
    richnessFontSize, betaValues, [rightPanelMinimum 0])

% Figure 2: surviving fraction at the end of deterministic burn-in.

burnSurvivingFraction = reshape( ...
    burnRichness / model.S, [nBeta, nTheta]);
burnPanelMinimum = min( ...
    burnSurvivingFraction, [], 'all', 'omitnan');
if ~isfinite(burnPanelMinimum) || burnPanelMinimum >= 1
    burnPanelMinimum = 0.9;
end

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
burnColors = lines(nTheta);
burnLegendLabels = compose('$\\theta = %.2g$', thetaValues);

figure(2)
clf
set(gcf, 'Color', 'w', 'Position', [180 100 900 850])
hold on
burnFigureHandles = gobjects(nTheta, 1);
for g = 1:nTheta
    burnFigureHandles(g) = plot( ...
        betaValues, burnSurvivingFraction(:, g), ...
        'Color', burnColors(g, :), 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerSize', 15, ...
        'MarkerFaceColor', burnColors(g, :), ...
        'MarkerEdgeColor', burnColors(g, :), 'LineWidth', 3);
end
set(gca, 'XScale', 'log', 'FontSize', 60, 'LineWidth', 1.5)
xlabel('$\beta$')
ylabel('Fraction of surviving species (before noise)')
xlim([min(betaValues) max(betaValues)])
ylim([burnPanelMinimum 1])
xticks([1e-2 1 1e2])
pbaspect([1 1 1])
box on
legend(burnFigureHandles, burnLegendLabels, ...
    'Location', 'best', 'Orientation', 'vertical', 'FontSize', 42)

function [N, R, T, firstNonfiniteTime, ...
    firstNegativeResourceTime] = simulateBurnIn( ...
    beta, theta, maximumDt, explicitSafetyFactor, burnTime, model)

N = model.N0 * ones(1, model.S);
R = model.R0;
T = model.T0 * ones(1, model.S);
firstNonfiniteTime = NaN;
firstNegativeResourceTime = NaN;
currentTime = 0;

while currentTime < burnTime
    [stepSize, invalidRate] = explicitStepSize(N, beta, theta, maximumDt, ...
        explicitSafetyFactor, model);
    if invalidRate
        firstNonfiniteTime = currentTime;
        break
    end
    stepSize = min(stepSize, burnTime - currentTime);

    [N, R, T] = deterministicStep(N, R, T, beta, theta, stepSize, model);
    currentTime = currentTime + stepSize;
    if any(~isfinite(N), 'all') || ~isfinite(R) || ...
            any(~isfinite(T), 'all')
        firstNonfiniteTime = currentTime;
        break
    end
    if R < 0 || any(T < 0, 'all')
        firstNegativeResourceTime = currentTime;
        break
    end

end
end

function [finalRichness, firstNonfiniteTime, ...
    firstNegativeResourceTime] = ...
    simulateObservationEnsemble(N, R, T, beta, theta, ...
    sigma, observationTime, maximumDt, explicitSafetyFactor, ...
    nTrajectories, model, stream)

% Each row is an independent stochastic trajectory.
N = repmat(N, nTrajectories, 1);
R = repmat(R, nTrajectories, 1);
T = repmat(T, nTrajectories, 1);
nonfiniteFlag = false(nTrajectories, 1);
firstNonfiniteTime = nan(nTrajectories, 1);
negativeResourceFlag = false(nTrajectories, 1);
firstNegativeResourceTime = nan(nTrajectories, 1);
currentTime = zeros(nTrajectories, 1);

while true
    invalidFlag = nonfiniteFlag | negativeResourceFlag;
    activeRows = ~invalidFlag & any(N > 0, 2) & ...
        currentTime < observationTime;
    if ~any(activeRows)
        break
    end

    [stepSizes, invalidActiveRates] = explicitStepSize( ...
        N(activeRows, :), beta, theta, ...
        maximumDt, explicitSafetyFactor, model);
    if any(invalidActiveRates)
        activeIndices = find(activeRows);
        invalidRateRows = activeIndices(invalidActiveRates);
        nonfiniteFlag(invalidRateRows) = true;
        firstNonfiniteTime(invalidRateRows) = currentTime(invalidRateRows);
        N(invalidRateRows, :) = 0;
        R(invalidRateRows) = 0;
        T(invalidRateRows, :) = 0;
        continue
    end
    stepSizes = min(stepSizes, observationTime - currentTime(activeRows));
    brownianIncrement = sqrt(stepSizes) .* ...
        randn(stream, nnz(activeRows), model.S);
    [activeN, activeR, activeT] = observationStep( ...
        N(activeRows, :), R(activeRows), T(activeRows, :), ...
        beta, theta, sigma, stepSizes, brownianIncrement, model);
    N(activeRows, :) = activeN;
    R(activeRows) = activeR;
    T(activeRows, :) = activeT;
    currentTime(activeRows) = currentTime(activeRows) + stepSizes;

    nonfiniteNow = false(nTrajectories, 1);
    nonfiniteNow(activeRows) = ~all(isfinite(activeN), 2) | ...
        ~isfinite(activeR) | ~all(isfinite(activeT), 2);
    newNonfinite = nonfiniteNow & ~nonfiniteFlag;
    firstNonfiniteTime(newNonfinite) = currentTime(newNonfinite);
    nonfiniteFlag = nonfiniteFlag | nonfiniteNow;

    negativeResourceNow = false(nTrajectories, 1);
    negativeResourceNow(activeRows) = activeR < 0 | any(activeT < 0, 2);
    newNegativeResource = negativeResourceNow & ~negativeResourceFlag;
    firstNegativeResourceTime(newNegativeResource) = ...
        currentTime(newNegativeResource);
    negativeResourceFlag = negativeResourceFlag | negativeResourceNow;

    % Isolate invalid rows so they cannot contaminate valid trajectories.
    invalidNow = nonfiniteNow | negativeResourceNow;
    N(invalidNow, :) = 0;
    R(invalidNow) = 0;
    T(invalidNow, :) = 0;

end

finalRichness = sum(N > 0, 2);
invalidFlag = nonfiniteFlag | negativeResourceFlag;
finalRichness(invalidFlag) = NaN;
end

function [N, R, T] = observationStep( ...
    N, R, T, beta, theta, sigma, dt, brownianIncrement, model)

[dN, dR, dT] = deterministicRates(N, R, T, beta, theta, model);
N = N + dN .* dt + sigma .* N .* brownianIncrement;
R = R + dR .* dt;
T = T + dT .* dt;
N(N < model.extinctionThreshold) = 0;
end

function [N, R, T] = deterministicStep(N, R, T, beta, theta, dt, model)

[dN, dR, dT] = deterministicRates(N, R, T, beta, theta, model);
N = N + dN .* dt;
R = R + dR .* dt;
T = T + dT .* dt;
N(N < model.extinctionThreshold) = 0;
end

function [dt, invalidRateRows] = explicitStepSize( ...
    N, beta, theta, maximumDt, explicitSafetyFactor, model)

% Bound every explicit resource-removal coefficient away from one.
localRemovalRate = beta + ...
    model.alpha .* model.gamma * theta .* N;
sharedRemovalRate = 1 + model.S * beta + ...
    model.gamma * (1 - theta) * sum(model.alpha .* N, 2);
invalidRateRows = ~all(isfinite(localRemovalRate), 2) | ...
    ~isfinite(sharedRemovalRate);
validLocalRates = localRemovalRate(~invalidRateRows, :);
validSharedRates = sharedRemovalRate(~invalidRateRows);
maximumRemovalRate = max(max(validLocalRates, [], 2), validSharedRates);
dt = min(maximumDt, explicitSafetyFactor ./ maximumRemovalRate);
end

function [dN, dR, dT] = deterministicRates(N, R, T, beta, theta, model)

% Continuous drift shared by deterministic and stochastic explicit steps.
dN = (model.alpha .* (theta * T + (1 - theta) * R) - model.k) .* N;
dT = beta * (R - T) - ...
    model.alpha .* model.gamma * theta .* T .* N;
dR = model.Rstar - R - beta * (model.S * R - sum(T, 2)) - ...
    model.gamma * R * (1 - theta) .* sum(model.alpha .* N, 2);
end

function [sampleTimes, abundanceHistory, N, R, T] = ...
    simulateRecordedPhase(N, R, T, beta, theta, duration, ...
    outputInterval, maximumDt, explicitSafetyFactor, sigma, model, stream)

% Record one trajectory at fixed physical times without storing solver steps.
sampleTimes = (0:outputInterval:duration).';
if sampleTimes(end) < duration
    sampleTimes(end + 1, 1) = duration;
end
abundanceHistory = zeros(numel(sampleTimes), model.S);
abundanceHistory(1, :) = N;
sampleIndex = 2;
currentTime = 0;
timeTolerance = 10 * eps(max(1, duration));
noiseActive = ~isempty(stream) && sigma > 0;

while currentTime < duration
    % Record a requested time reached to within floating-point precision.
    if sampleIndex <= numel(sampleTimes) && ...
            currentTime >= sampleTimes(sampleIndex) - timeTolerance
        abundanceHistory(sampleIndex, :) = N;
        sampleIndex = sampleIndex + 1;
        continue
    end

    % Once every species is extinct, all later abundances remain zero.
    if ~any(N > 0)
        abundanceHistory(sampleIndex:end, :) = repmat( ...
            N, numel(sampleTimes) - sampleIndex + 1, 1);
        break
    end

    [stepSize, invalidRate] = explicitStepSize( ...
        N, beta, theta, maximumDt, explicitSafetyFactor, model);
    if invalidRate
        error('A nonfinite adaptive rate occurred at time %.6g.', currentTime)
    end
    stepSize = min(stepSize, duration - currentTime);
    if sampleIndex <= numel(sampleTimes)
        stepSize = min( ...
            stepSize, sampleTimes(sampleIndex) - currentTime);
    end

    if noiseActive
        brownianIncrement = sqrt(stepSize) * ...
            randn(stream, 1, model.S);
        [N, R, T] = observationStep( ...
            N, R, T, beta, theta, sigma, stepSize, ...
            brownianIncrement, model);
    else
        [N, R, T] = deterministicStep( ...
            N, R, T, beta, theta, stepSize, model);
    end
    currentTime = currentTime + stepSize;

    if any(~isfinite(N), 'all') || ~isfinite(R) || ...
            any(~isfinite(T), 'all')
        error('A nonfinite recorded state occurred at time %.6g.', currentTime)
    end
    if R < 0 || any(T < 0, 'all')
        error('A negative recorded resource occurred at time %.6g.', currentTime)
    end

    if sampleIndex <= numel(sampleTimes) && ...
            currentTime >= sampleTimes(sampleIndex) - timeTolerance
        abundanceHistory(sampleIndex, :) = N;
        sampleIndex = sampleIndex + 1;
    end
end
end

function formatRichnessPanel(yLabelText, fontSize, betaValues, yLimits)

set(gca, 'XScale', 'log', 'FontSize', fontSize, 'LineWidth', 1.5)
xlabel('$\beta$')
ylabel(yLabelText)
xlim([min(betaValues) max(betaValues)])
ylim(yLimits)
xticks([1e-2 1 1e2])
pbaspect([1 1 1])
box on
end
