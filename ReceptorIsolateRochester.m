% ReceptorIsolateRochester
%
% Find modulations that isolates various photopigments, for various device
% models. The engine that drives this demo is ReceptorIsolate.
%
% This version is set up to aim at the ganglion cell AO experiments in
% Rochester.
%
% 8/8/16  dhb, jem  Wrote from ReceptorIsolateDemo in SilentSubstitutionToolbox

%% Clear and close
clear; close all;

%% Run in the directory that contains this file.
cd(fileparts(mfilename('fullpath')));

%% Read in receptor sensitivities
%
% We'll do this in the Rochester way
Photoreceptor.wavelengths = (350 : 1 : 700)';
Photoreceptor.macaque_scone      = StandardTemplate(430, Photoreceptor.wavelengths);
Photoreceptor.macaque_mcone      = StandardTemplate(530, Photoreceptor.wavelengths);
Photoreceptor.macaque_lcone      = StandardTemplate(561, Photoreceptor.wavelengths);
Photoreceptor.macaque_rod        = StandardTemplate(491, Photoreceptor.wavelengths);
Photoreceptor.macaque_melanopsin = StandardTemplate(480, Photoreceptor.wavelengths);

% Build up our T_receptors matrx for L, M, and S cones
S = WlsToS(Photoreceptor.wavelengths);
T_receptors = [Photoreceptor.macaque_lcone' ; Photoreceptor.macaque_mcone' ; Photoreceptor.macaque_scone'];

% Plot
theFig1 = figure; clf; hold on
plot(SToWls(S),T_receptors,'LineWidth',2);
xlabel('Wavelength (nm)')
ylabel('Sensitivity');
title('Normalized photoreceptor sensitivities');

%% Get LED primary spectra
%
% These get read out of a list of text files.
theLEDDir = 'LEDSpectra';
theLEDFiles = {'M405L2.txt', 'M455L3.txt', 'M530L3.txt' 'M590L3.txt'};
B_primary = zeros(length(Photoreceptor.wavelengths),length(theLEDFiles));
theFigLEDs = figure; clf; hold on
for ii = 1:length(theLEDFiles)
    theLEDFile = fullfile(theLEDDir,theLEDFiles{ii});
    theLEDSpectraRaw{ii} = dlmread(theLEDFile,'\t',1,0);
    theLEDSpecraSplined{ii} = interp1(theLEDSpectraRaw{ii}(:,1),theLEDSpectraRaw{ii}(:,2),Photoreceptor.wavelengths,'linear',0);
    B_primary(:,ii) = theLEDSpecraSplined{ii}*S(2);
    plot(theLEDSpectraRaw{ii}(:,1),theLEDSpectraRaw{ii}(:,2),'r','LineWidth',3);
    plot(Photoreceptor.wavelengths,theLEDSpecraSplined{ii},'k','LineWidth',1);
end
xlabel('Wavelength');
ylabel('LED Power');

%% Specify which receptor to target and set up background
%
%   1 = L, 2 = M, 3 = S
whichReceptorsToTarget = [2];
whichReceptorsToIgnore = [];
whichReceptorsToMinimize = [];
whichPrimariesToPin = [];
backgroundPrimary = 0.5*ones(length(theLEDFiles),1);
primaryHeadRoom = 0.02;
maxPowerDiff = 10000;

%% Specify ambient 
%
% If you have an actual measurement you can read it in here.
ambientSpd = zeros(size(Photoreceptor.wavelengths));

% User chooses whether to maximize contrast in targeted receptor classes or
% or get it as close to a specified value as possible.
%
% If we target, here we specify the same contrast for all targeted classes.
% This is not necessary, they can differ.  It just makes the demo code a
% bit simpler to yoke them since we only have to prompt for one number.
maximizeTargetContrast = GetWithDefault('\tMaximize contrast? [1 = yes, 0 = no]', 1);
if maximizeTargetContrast
    desiredContrast = [];
else
    desiredContrast = GetWithDefault('\tDesired contrast (applies to all targeted classes)?', 0.45)*ones(size(whichReceptorsToTarget));
end
fprintf('\n');

%% Call the optimization routine.
%
% Careful examaination of the arguments will reveal that the initialGuess
% for the primaries is set to the background value for the primaries, so
% that the constraints are all met at the start of the search.  The
% optimization routine is much happier when it is done this way -- bad
% things happen if you start with a guess that violates constraints.
modulationPrimary = ReceptorIsolate(T_receptors,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
    B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
    primaryHeadRoom, maxPowerDiff, desiredContrast, ambientSpd);

%% Compute the contrasts that we got.
backgroundReceptors = T_receptors*(B_primary*backgroundPrimary + ambientSpd);
modulationReceptors = T_receptors*B_primary*(modulationPrimary - backgroundPrimary);
contrastReceptors = modulationReceptors ./ backgroundReceptors;
fprintf('Obtained contrasts\n');
for j = 1:size(T_receptors,1)
    fprintf('\tClass %d: contrast = %0.4f\n',j,contrastReceptors(j));
end

% Modulation spectra
theFig2 = figure; hold on
plot(SToWls(S),B_primary*modulationPrimary,'r','LineWidth',2);
plot(SToWls(S),B_primary*backgroundPrimary,'k','LineWidth',2);
title('Modulation spectra');
xlim([380 780]);
xlabel('Wavelength');
ylabel('Power');
pbaspect([1 1 1]);

% Plot
theFig3 = figure; hold on
plot(modulationPrimary,'r','LineWidth',2);
plot(backgroundPrimary,'k','LineWidth',2);
title('Primary settings');
xlim([0 length(backgroundPrimary)]);
ylim([0 1]);
xlabel('Primary Number (nominal)');
ylabel('Setting');
