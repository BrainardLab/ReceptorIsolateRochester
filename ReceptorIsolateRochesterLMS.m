% ReceptorIsolateRochesterLMS
%
% Find modulations that do various things in LMS space, based on
% three primaries. The engine that drives this demo is ReceptorIsolate.
%
% This version is set up to aim at the ganglion cell AO experiments in
% Rochester.
%
% 08/08/16  dhb, jem  Wrote from ReceptorIsolateDemo in SilentSubstitutionToolbox
% 02/02/19  dhb       Modified version more matched to what Rochester currently wants to do.

%% Clear and close
clear; close all;

%% Run in the directory that contains this file.
cd(fileparts(mfilename('fullpath')));

%% Read in receptor sensitivities
%
% Here we use the CIE 2 degree fundamentals and
% spline to the wavelength spacing we feel like
% using.
wavelengths = (380:1:750)';
S = WlsToS(wavelengths);
load T_cones_ss2
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S);

% Specify luminance as weighted sum of L and M cones.  We want this to be
% a weighted sum of the L and M cones, as that matches how we typically
% think about post-receptoral channels.  To enforce this, we start with a
% luminance spectral sensitivity and fit it with a weigthed sum of the L
% and M cones that we're using.
load T_CIE_Y2
T_lum_tabulated = SplineCmf(S_CIE_Y2,T_CIE_Y2,S);
lumFactorsLM = (T_cones(1:2,:)'\T_lum_tabulated');
T_lum = (T_cones(1:2,:)'*lumFactorsLM)';

% Build up our T_receptors matrx for L, M, and S cones plus luminance as
% fourth row.
T_receptors = [T_cones ; T_lum];

% Plot. In this case, the plot also shows that the fit of the tabulated
% luminance as weighted sum of L and M cones is very good.  That's because
% in the CIE system, the tabulated luminance was constructed as a weighted
% sum of the L and M cones.
theFig1 = figure; clf; hold on
plot(SToWls(S),T_receptors(1,:)','r','LineWidth',3);
plot(SToWls(S),T_receptors(2,:)','g','LineWidth',3);
plot(SToWls(S),T_receptors(3,:)','b','LineWidth',3);
plot(SToWls(S),T_receptors(4,:)','k','LineWidth',3);
plot(SToWls(S),T_lum_tabulated','y-','LineWidth',2);
legend({'L cones', 'M_cones', 'S_cones', 'Luminance (fit)', 'Luminance (tabulated)'});
xlabel('Wavelength (nm)')
ylabel('Sensitivity');
title('Photoreceptor sensitivities');

%% Get LED primary spectra
%
% These get read out of a list of text files that provide the relative
% spectra, and the relative spectra are then scaled to be consistent with
% the specified maximum power in microwatts.
theLEDDir = 'LEDSpectra';
theLEDFiles = {'M405L2.txt','M530L3.txt' 'M590L3.txt'};
theLEDMaxPowerUW = [10 22 16];
theLEDPlotColors = ['b' 'g' 'r'];
B_primary = zeros(length(wavelengths),length(theLEDFiles));
theFigLEDs = figure; clf;
for ii = 1:length(theLEDFiles)
    % Read in relative spectra and then scale for max power.
    theLEDFile = fullfile(theLEDDir,theLEDFiles{ii});
    theRelativeLEDSpectraRaw{ii} = dlmread(theLEDFile,'\t',1,0);
    theRelativeLEDSpecraSplined{ii} = interp1(theRelativeLEDSpectraRaw{ii}(:,1),theRelativeLEDSpectraRaw{ii}(:,2),wavelengths,'linear',0);

    % Get constant by integrating relative spectra over wavelength.
    % Psychtoolbox, which underlies this code, likes to "think" of
    % power in units of power/wl-band rather than power/nm and we
    % follow that convention here.  If we are operating on 1 nm
    % wavelength sampling, as we were when this code was first written,
    % the two conventions collapse to the same thing.  Be a little careful
    % about this if you change the way you bring in your primary spectra.
    relativePower = sum(theRelativeLEDSpecraSplined{ii});
    theLEDSpectraScaled{ii} = theLEDMaxPowerUW(ii)*theRelativeLEDSpecraSplined{ii}/relativePower;
    
    % Build up matrix with primaries in the columns
    B_primary(:,ii) = theLEDSpectraScaled{ii};
    
    % Plot.
    %
    % Relative spectra, and splined version.
    subplot(1,2,1); hold on;
    plot(theRelativeLEDSpectraRaw{ii}(:,1),theRelativeLEDSpectraRaw{ii}(:,2),theLEDPlotColors(ii),'LineWidth',3);
    plot(wavelengths,theRelativeLEDSpecraSplined{ii},'k-','LineWidth',1);
    
    % Scaled version
    subplot(1,2,2); hold on;
    plot(wavelengths,theLEDSpectraScaled{ii},theLEDPlotColors(ii),'LineWidth',3);

end
subplot(1,2,1); hold on;
xlabel('Wavelength');
ylabel('Relative LED Power');
subplot(1,2,2); hold on;
xlabel('Wavelength');
ylabel('LED Power');

%% Specify ambient light
%
% The ambient light is not the background for the modulations. Rather it is
% the light delivered to the eyes when the settings on the primaries are all
% zero.  With LED's, this may well be zero.  With typical computer displays there
% is often sone light even when all the DAC values are 0, which is why we allow
% specification of something here. 
%
% Currently we assume it is dark when all primaries are set to 0.
ambientSpd = zeros(size(wavelengths));

%% Specify background in terms of primaries
backgroundPrimary = [0.45 0.55 0.50]';

%% Specify which receptor to target and set up background
%
%   1 = L, 2 = M, 3 = S, 4 = lum
whichReceptorsToTarget = [1 2 3 4];
whichReceptorsToIgnore = [];
whichReceptorsToMinimize = [];
whichPrimariesToPin = [];
desiredContrast = [0.5 0.5 0.5 0.5];
primaryHeadRoom = 0.02;
maxPowerDiff = 10000;

% % User chooses whether to maximize contrast in targeted receptor classes or
% % or get it as close to a specified value as possible.
% %
% % If we target, here we specify the same contrast for all targeted classes.
% % This is not necessary, they can differ.  It just makes the demo code a
% % bit simpler to yoke them since we only have to prompt for one number.
% maximizeTargetContrast = GetWithDefault('\tMaximize contrast? [1 = yes, 0 = no]', 1);
% if maximizeTargetContrast
%     desiredContrast = [];
% else
%     desiredContrast = GetWithDefault('\tDesired contrast (applies to all targeted classes)?', 0.45)*ones(size(whichReceptorsToTarget));
% end
% fprintf('\n');

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
