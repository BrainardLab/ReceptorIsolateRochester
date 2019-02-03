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

% Get CIE 1931 XYZ, just so we can compute xy chromaticity of the
% background.
load T_xyz1931;
T_xyz = SplineCmf(S_xyz1931,T_xyz1931,S);

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
% the light delivered to the eyes when the settings on the primaries are
% all zero.  With LED's, this may well be zero.  With typical computer
% displays there is often sone light even when all the DAC values are 0,
% which is why we allow specification of something here.
%
% This is specified as a spectrum, as in our lab we typically measure it
% directly and sometimes it is not a combination of the primarie (for
% example, if there is ambient in the room that that comes from some other
% source).
%
% Currently we assume it is dark when all primaries are set to 0, and this
% is probably fine as an approximation.
ambientSpd = zeros(size(wavelengths));

%% Specify background in terms of primaries
backgroundPrimary = [0.45 0.55 0.50]';

%% Specify which receptor to target and set up background
%
%   1 = L, 2 = M, 3 = S, 4 = lum
%
% Examples: Uncomment the two lines for each example and comment the same
% two lines for the other examples to see each in action. 
%
% 1) Produce specified isochromtic contrast (akd luminance contrast).
% Contrast has same value, here 50% on L, M, and S cones. Also 50% 
% luminance contrast. This approach will work as long as the specified
% contrasts are within gamut of the device. How much contrast you can
% get depends on how the background is set, on the ambient spectrum
% specified above, and on how much primary headroom we specify below.
% To check whether what was specified is in gamut, compare the obtained
% contrasts printed out when the program runs to what you requested. If
% they agree, then you're good. If they differ, then the requested
% modulation was not in gamut.
%   whichReceptorsToTarget = [1 2 3 4];
%   desiredContrast = [0.5 0.5 0.5 0.5];
%
% 2) Produce maximum luminance constrast. You can first run with all four
% receptors targeted and desired contrast as the empty matrix. This will
% give you a sense of the max contrast, but it doesn't enforce the
% isochromatic constraint. But you'll see that the answer is about 80%.
% Then you experiment.  Set desired contrast with all entries at 0.8 and
% see what comes back. You don't quite get 80% on all four, but it's close.
% then try 0.78. That works. Then 0.79. Not quite. So the answer to within
% a percent is 0.78. I agree, this is a little ugly but it would take me a
% while to make finding max contrast in a particular direction automatic.
% And, it only took me about a minute to find the 78% number, so it is
% perfectly practical.
%  whichReceptorsToTarget = [1 2 3 4];
%  desiredContrast = [0.78 0.78 0.78 0.78];
%
% 3) Isoluminant modulation. This works because it enforces that S cone and
% luminance contrast are zero, becuase neither 3 or 4 are included in the
% whichReceptorsToTarget list.  The contrast on these two is maximized
% given the contraints, and that comes out to producing an -L+M modulation.
% Note that luminance is defined by a weighted sum of L and M cone
% excitations, not L and M contrasts, so the ratio here is not exactly the 
% the same as that in the variable lumLMFactors. I should say that actual
% isoluminance probably doesn't work quite the way it is specified, because
% it probably does depend on adaptation to the background. None-the-less,
% this is pretty standard. And probably close enough as long as the
% background is somewhere close to neutral.
%  whichReceptorsToTarget = [1 2]T
%  desiredContrast = [];
%
% 4) Receptor isolating.  Specify which receptor you want to isolate and to
% ignore luminance. Here what works best is to get rid of the luminance
% row in the T_receptors matrix.
T_receptors = T_receptors(1:3,:);
whichReceptorsToTarget = [3];
desiredContrast = [];

%% Specify headroom
%
% This parameter says how close to the edge of the gamut we allow a
% modulation to get.  So if we specify 0.02, the primary values can range
% between 0.02 and 0.98. I like to leave a little headroom, becuase devices
% can get a little wonky right at their extrema, and sometimes the gamut
% change over time and I like to avoid having to change stimuli if the
% gamut gets a little smaller.
primaryHeadRoom = 0.02;

%% Other parameters
%
% This really doesn't do anything here, leave it as a big number.
maxPowerDiff = 10000;

% You don't need to much with these.
whichReceptorsToIgnore = [];
whichReceptorsToMinimize = [];
whichPrimariesToPin = [];

%% Call the optimization routine.
%
% Careful examaination of the arguments will reveal that the initialGuess
% for the primaries is set to the background value for the primaries, so
% that the constraints are all met at the start of the search.  The
% optimization routine is much happier when it is done this way -- bad
% things happen if you start with a guess that violates constraints.
positiveModulationPrimary = ReceptorIsolate(T_receptors,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
    B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
    primaryHeadRoom, maxPowerDiff, desiredContrast, ambientSpd);

%% Compute change in primaries.
% For a sinusoidal modulation, you'd multipy deltaPrimary by the sinusoidal
% value (ranging between -1 and 1) and add this to the background primary
deltaPrimary = positiveModulationPrimary-backgroundPrimary;
negativeModulationPrimary = backgroundPrimary-deltaPrimary;

%% Compute the contrasts that we got.
% 
% We put back luminance into T_receptors here, which we take out in some of
% the examples above.  That just lets us see luminance contrast even if we
% didn't specify it when we computed.  Note that L and M cone isolating
% modulations also produce luminance contrast.
T_receptors = [T_cones ; T_lum];
backgroundReceptors = T_receptors*(B_primary*backgroundPrimary + ambientSpd);
modulationReceptors = T_receptors*B_primary*(positiveModulationPrimary - backgroundPrimary);
contrastReceptors = modulationReceptors ./ backgroundReceptors;
fprintf('Obtained contrasts\n');
for j = 1:size(T_receptors,1)
    fprintf('\tClass %d: contrast = %0.4f\n',j,contrastReceptors(j));
end

%% Compute and report CIE xy chromaticity of background
% 
% Note that luminance (Y) is in arbitrary units because we
% have not scaled XYZ to match units of power used to in
% specifying LEDs.
%
% Equal energy white xy chromaticity is 0.33, 0.33.
backgroundXYZ = T_xyz*B_primary*backgroundPrimary;
backgroundxyY = XYZToxyY(backgroundXYZ);
fprintf('Background 1931 xy chromaticity: %0.3f, %0.3f\n',backgroundxyY(1),backgroundxyY(2));

% Plot some spectra
theFig2 = figure; hold on
plot(SToWls(S),B_primary*positiveModulationPrimary,'r','LineWidth',2);
plot(SToWls(S),B_primary*backgroundPrimary,'k','LineWidth',2);
title('Modulation spectra');
xlim([380 780]);
xlabel('Wavelength');
ylabel('Power');
pbaspect([1 1 1]);

% Plot primaries for background as well as high and low end of modulation
theFig3 = figure; hold on
plot(positiveModulationPrimary,'ro','MarkerFaceColor','r','MarkerSize',12);
plot(negativeModulationPrimary,'go','MarkerFaceColor','g','MarkerSize',12);
plot(backgroundPrimary,'ko','MarkerFaceColor','k','MarkerSize',12);
title('Primary values');
xlim([0 length(backgroundPrimary)+ 1]);
ylim([0 1]);
xlabel('Primary Number (nominal)');
ylabel('Primary Value');
