function ActionSpectra
%{
    Function to plot action spectra of photoreceptor sensitivites, light sources, flurophores, 
    filters and others. 

    Defined spectra:
    Intrinsic Opsins:
        - macaque photoreceptor action spectra - Schnapf etal 1988 Vis Neurosci
        - mouse photoreceptor action spectra - Wang etal 2011 J Neurosci
        - marmoset photoreceptor action spectra - In the marmoset three alleles code opsins that are 
            most sensitive to 543, 556, or 563 nm (Travis et al., 1988; Tovée et al., 1992; 
            Williams et al., 1992; Hunt et al., 1993; Shyue et al., 1995)

    Thorlabs LED:
        - M365L2 - ThorLabs
        - M405L2 - ThorLabs
        - M420L3 - ThorLabs
        - M455L3 - ThorLabs
        - M470L3 - ThorLabs
        - M490L3 - ThorLabs
        - M505L3 - ThorLabs
        - M530L3 - ThorLabs
        - M565L3 - ThorLabs
        - M590L3 - ThorLabs
        - M595L3 - ThorLabs
        - M617L3 - ThorLabs
        - M625L3 - ThorLabs
        - MCWHL5 - ThorLabs (cool white)
        - MWWHL4 - ThorLabs (warm white)

    Flurophore:
        - Alexa488      - Life Technologies
        - Alexa568      - Life Technologies
        - YFP           - Life Technologies
        - CFP           - Life Technologies
        - GFP           - 
        - Citrine       -
        - Rhodamine Red -
        - Texas Red     -
        - DAPI          - Life Technologies
        - Fluorescein   - Life Technologies
        - mCherry       -
        - tdTomato      - chroma
        - Cy3           - chroma
        - GCaMP5g       - Akerboom etal 2012 J Neruosci
        - RCaMP1e       - Akerboom etal 2013 Front Molec Neurosci, Abs: Fig 2; Ems: Fig 8D
        - RCaMP2        - Inoue etal 2015 Nat Methods
        - R-GECO1       - Zhao etal 2011 Science (1P spectra), Dana etal 2014 (SFN poster, 2P spectrum)
        - peredox       - Hung etal 2011 Cell Metabol. NADH-NAD+ FLIM sensor
        - DsRed         - (from Semrock SearchLight)

    Filters:
        - BP545/25       - Zeiss/Delta
        - BP605/70       - Zeiss/Delta
        - EO65141_455_10 - Edmund Optics
        - EO65160_568_10 - Edmund Optics
        - ET535_70       - Chroma
        - FF01_448_20    - Semrock
        - FF01_460_14    - Semrock
        - FF01_469_35    - Semrock
        - FF01_470_28    - Semrock
        - FF01_475_42    - Semrock
        - FF01_496_LP    - Semrock
        - FF01_500_10    - Semrock
        - FF01_500_15    - Semrock
        - FF01_504_12    - Semrock
        - FF01_520_35    - Semrock
        - FF01_520_70    - Semrock
        - FF01_525_45    - Semrock
        - FF01_530_55    - Semrock
        - FF01_550_49    - Semrock
        - FF01_565_24    - Semrock
        - FF01_571_72    - Semrock
        - FF01_580_60    - Semrock
        - FF01_582_75    - Semrock
        - FF01_630_92    - Semrock

        Zeiss Filter sets:
        - Filter Set 43 - BP545/25 & BP605/70 - cy3, tdTomato (https://www.micro-shop.zeiss.com/?s=227735774c82f69&l=en&p=us&f=f&a=v&b=f&id=000000-1114-101&o=)

    Biology:
        - channelrhodopsin 2 - Nagel etal 2003 PNAS
        - mouseWholeEyeTrasmittance - Henricksson etal 2010 Exp Eye Res
        - ReaChR - red-shifted ChR - Lin etal 2013 Nat Neurosci
        - C1V1 - red-shifted ChR - Yhizar etal 2011 Nature
        - ChrimsonR - red-shifted channelrhodopsin - Klapoetke etal 2014 Nat Methods
        - chronos - blue light activated channelrhodopsin - Klapoetke etal 2009 Nat Methods
        - jaws - red-shifted halorhodopsin - Chuong etal 2014 Nat Neurosci
        
    Other:
        - Solar irradiance at sea level (http://www.newport.com/Introduction-to-Solar-Radiation/411919/1033/content.aspx)
            - irradiance units (Wm^-2nm^-1)

    TODO:


    Changelog:
        - tidied code

    Written By: Kenny Cheong
    Date: 1 Janurary 2015

%}
% ---- Construction of spectra curves ----
    % construct photoreceptor sensitivity spectra
        wavelengths = (350 : 1 : 700)';
        Photoreceptor.WAVELENGTH = wavelengths; 
        
        % mouse sensitivities (Wang etal 2011 J Neurosci)
        Photoreceptor.mouse_scone      = StandardTemplate(360, wavelengths); 
        Photoreceptor.mouse_mcone      = StandardTemplate(508, wavelengths); 
        Photoreceptor.mouse_rods       = StandardTemplate(498, wavelengths);
        Photoreceptor.mouse_melanopsin = StandardTemplate(480, wavelengths); % lmax from Fig 1 Lall etal 2010 Neuron
        
        % macaque sensitivities (Schnapf etal 1988 Vis Neurosci)
        Photoreceptor.macaque_scone      = StandardTemplate(430, wavelengths);
        Photoreceptor.macaque_mcone      = StandardTemplate(530, wavelengths);
        Photoreceptor.macaque_lcone      = StandardTemplate(561, wavelengths);
        Photoreceptor.macaque_rod        = StandardTemplate(491, wavelengths);
        Photoreceptor.macaque_melanopsin = StandardTemplate(480, wavelengths);
        
        % marmoset sensitivities ()
        Photoreceptor.marmoset_scone      = StandardTemplate(423, wavelengths);
        Photoreceptor.marmoset_543        = StandardTemplate(543, wavelengths);
        Photoreceptor.marmoset_556        = StandardTemplate(556, wavelengths);
        Photoreceptor.marmoset_563        = StandardTemplate(563, wavelengths);
        Photoreceptor.marmoset_rod        = StandardTemplate(491, wavelengths); % ?
        Photoreceptor.marmoset_melanopsin = StandardTemplate(480, wavelengths); % ?       
        
        % optogenetic molecules
        Optogenetics.WAVELENGTH = wavelengths; 
        Optogenetics.ChR2  = StandardTemplate(460, wavelengths); % lmax from Nagel etal 2003 PNAS
        Optogenetics.CatCh = StandardTemplate(474, wavelengths); % lmax from Kleinlogel etal 2011 Nat Neurosci
            % lmax of 463 results in a better with with the standard template than the reported 474
            % according to supp. fig. 2, few data points were measured posibly resulting in poor lmax estimation

    % load in other biological data
        Biology = loadTablesFromFolder(fullfile(pwd, 'Biology'));
        
    % load in led spectra
        % automatically load led spectra from txt files in .\LED folder
        % see existing file for example format
        LED = loadTablesFromFolder(fullfile(pwd, 'LED'));
        
%         % special case normalisation
%         LED.M565L3_real.INTENSITY      = LED.M565L3_real.INTENSITY - min(LED.M565L3_real.INTENSITY);
%         LED.M565L3_real.NORM_INTENSITY = LED.M565L3_real.INTENSITY ./ max(LED.M565L3_real.INTENSITY);
        
    % load in dichroic 
        Optics = loadTablesFromFolder(fullfile(pwd, 'Optics'));
        
    % load in filters
        Filter = loadTablesFromFolder(fullfile(pwd, 'Filter'));
        
    % load in flurophores
        Flurophore = loadTablesFromFolder(fullfile(pwd, 'Flurophore'));
        
    % load in other
        Other = loadTablesFromFolder(fullfile(pwd, 'Other'));        
    
    % load in measured system spectra
        % measurements were taked at the animal pupil through system / stimulus aparatus and pellicle
        System = loadTablesFromFolder(fullfile(pwd, 'System'));

        % normalise curves
        System_fieldnames = fieldnames(System);
        for qq = 1:length(System_fieldnames)
            System.(System_fieldnames{qq}).INTENSITY = System.(System_fieldnames{qq}).INTENSITY - min(System.(System_fieldnames{qq}).INTENSITY);
            System.(System_fieldnames{qq}).NORM_INTENSITY = System.(System_fieldnames{qq}).INTENSITY ./ max(System.(System_fieldnames{qq}).INTENSITY);
        end
        
    % load in measured power calibration curves
        PowerCurves = loadTablesFromFolder(fullfile(pwd, 'PowerCurves'));
        
    % references
    References.macaque = 'Schnapf etal 1988 Vis Neurosci';
    References.mouse   = 'Wang etal 2011 J Neurosci';
    References.marmoset  = 'Tovee etal 1992 VisRes';
    References.ChR2      = 'Nagel etal 2003 PNAS';
    References.ChrimsonR = 'Klapoetke etal 2014 Nat Methods';
    References.Chronos   = 'Klapoetke etal 2014 Nat Methods';
    References.ReaChR = 'Lin etal 2013 Nat Neurosci';
    References.C1V1 = 'Yhizar etal 2011 Nature';
    References.jaws = 'Chuong etal 2014 Nat Neurosci';
    References.GCaMP5g = 'Akerboom etal 2012 J Neruosci';
    References.RCaMP1e = 'Akerboom etal 2013 Front Molec Neurosci';
    References.RCaMP2  = 'Inoue etal 2015 Nat Methods';
    References.RGECO1  = '1P: Zhao etal 2011 Science; 2P: Dana etal 2014';
    References.peredox = 'Hung etal 2011 Cell Metabol';
    
    References.mouseWholeEyeTrasmittance = 'Henricksson etal 2010 Exp Eye Res';

        
'pause'
% ==== Plots ====
% manually code figures as there are too many possilbe combinations of plots to justify a standard format or GUI

%% -- plot --
% histology
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % Fluorephores
        plot(Flurophore.DAPI.WAVELENGTH, Flurophore.DAPI.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-'); hold on
        plot(Flurophore.DAPI.WAVELENGTH, Flurophore.DAPI.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on        
        
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-'); hold on
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        plot(Flurophore.Fluorescein.WAVELENGTH, Flurophore.Fluorescein.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-'); hold on
        plot(Flurophore.Fluorescein.WAVELENGTH, Flurophore.Fluorescein.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on        

        plot(Flurophore.tdTomato.WAVELENGTH, Flurophore.tdTomato.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-'); hold on
        plot(Flurophore.tdTomato.WAVELENGTH, Flurophore.tdTomato.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        plot(Flurophore.Cy3.WAVELENGTH, Flurophore.Cy3.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-'); hold on
        plot(Flurophore.Cy3.WAVELENGTH, Flurophore.Cy3.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on        
        
        plot(Flurophore.Alexa647.WAVELENGTH, Flurophore.Alexa647.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-'); hold on
        plot(Flurophore.Alexa647.WAVELENGTH, Flurophore.Alexa647.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        plot(Flurophore.Cy5.WAVELENGTH, Flurophore.Cy5.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-'); hold on
        plot(Flurophore.Cy5.WAVELENGTH, Flurophore.Cy5.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Histology spectra')
        xlim([350 700])
        ylim([0 1.01])

        l = legend('DAPI - excitation', 'DAPI - emission',...
            'GCaMP5g - excitation', 'GCaMP5g - emission', 'Fluorescein - excitation', 'Fluorescein - emission',...
            'tdTomato - excitation', 'tdTomato - emission', 'Cy3 - excitation', 'Cy3 - emission',...
            'Alexa647 - excitation', 'Alexa647 - emission', 'Cy5 - excitation', 'Cy5 - emission');
               
        set(l, 'Fontsize', 12, 'Box', 'off')     
        
        save_name = 'histology.pdf';
        save_dir  = fullfile('.', 'Figures');
        save_path = fullfile(save_dir, save_name);
        exportpdf(f, save_path)


%% -- Plot --
% compute stimulus equivalency for mcone for 455 vs 565
    % mcone and 455
    mouse_mcone_455 = StandardTemplate(508, LED.M455L3.WAVELENGTH);
    mouse_mcone_455_product  = mouse_mcone_455 .* LED.M455L3.NORM_INTENSITY;
    mouse_mcone_455_integral = trapz(LED.M455L3.WAVELENGTH, mouse_mcone_455_product);
    
    mouse_mcone_565 = StandardTemplate(508, LED.M565L3.WAVELENGTH);
    mouse_mcone_565_product  = mouse_mcone_565 .* LED.M565L3.NORM_INTENSITY;
    mouse_mcone_565_integral = trapz(LED.M565L3.WAVELENGTH, mouse_mcone_565_product);
     
    ratio_455_565 = mouse_mcone_565_integral / mouse_mcone_455_integral; % 565 as a multiple of 455
    ratio_455_565_inv = 1/ratio_455_565; % ratio of 565 for 1 unit intensity of 455
    
    % approximate power using peak led wavelength
    % use photon_calculator
    
    f = figure('position', [541         152        1386         719]);
        plot(LED.M455L3.WAVELENGTH, mouse_mcone_455); hold on
        plot(LED.M455L3.WAVELENGTH, LED.M455L3.NORM_INTENSITY); hold on
        plot(LED.M455L3.WAVELENGTH, mouse_mcone_455_product); hold on

        plot(LED.M565L3.WAVELENGTH, LED.M565L3.NORM_INTENSITY); hold on
        plot(LED.M565L3.WAVELENGTH, mouse_mcone_565_product); hold on    

        xlim([350 700])

%% -- Plot -- 
    % compare ChR2 activation from 455 vs 565
    % conclusion: relative activation of ChR2 from 455 : 565 => 1 : 0.3850
    
    % ChR2 and 455
    chr2_455 = StandardTemplate(460, LED.M455L3.WAVELENGTH);
    chr2_455_product  = chr2_455 .* LED.M455L3.NORM_INTENSITY;
    chr2_455_integral = trapz(LED.M455L3.WAVELENGTH, chr2_455_product);
    
    % ChR2 and 565
    chr2_565 = StandardTemplate(460, LED.M565L3.WAVELENGTH);
    chr2_565_product  = chr2_565 .* LED.M565L3.NORM_INTENSITY;
    chr2_565_integral = trapz(LED.M565L3.WAVELENGTH, chr2_565_product);
     
    ratio_455_565     = chr2_565_integral / chr2_455_integral; % 565 as a multiple of 455
    ratio_455_565_inv = 1/ratio_455_565; % ratio of 565 for 1 unit intensity of 455
    
    % approximate power using peak led wavelength
    % use photon_calculator
    
    f = figure('position', [541         152        1386         719]);
        plot(LED.M455L3.WAVELENGTH, chr2_455); hold on
        plot(LED.M455L3.WAVELENGTH, LED.M455L3.NORM_INTENSITY); hold on
        plot(LED.M455L3.WAVELENGTH, chr2_455_product); hold on

        plot(LED.M565L3.WAVELENGTH, LED.M565L3.NORM_INTENSITY); hold on
        plot(LED.M565L3.WAVELENGTH, chr2_565_product); hold on    

        xlim([350 700])

%% -- Plot --
    % plot ChR2 with multiple led spectra to design a combination led stimulus that maximally activates ChR2
    f = figure('position', [541         152        1386         719]);
        plot(Optogenetics.WAVELENGTH, Optogenetics.ChR2, '--'), hold on
    
        plot(LED.M405L2.WAVELENGTH, LED.M405L2.NORM_INTENSITY); hold on
        plot(LED.M420L3.WAVELENGTH, LED.M420L3.NORM_INTENSITY); hold on
        plot(LED.M455L3.WAVELENGTH, LED.M455L3.NORM_INTENSITY); hold on
        plot(LED.M470L3.WAVELENGTH, LED.M470L3.NORM_INTENSITY); hold on
        plot(LED.M490L3.WAVELENGTH, LED.M490L3.NORM_INTENSITY); hold on
        plot(LED.M505L3.WAVELENGTH, LED.M505L3.NORM_INTENSITY); hold on
    
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350 700])
        ylim([0 1.01])

        l = legend('ChR2', '405', '420', '455', '470', '490', '505');
        set(l, 'Fontsize', 12, 'Box', 'off')
        
%% -- Plot --
    % Compare 565 datasheet and measured spectra
    f = figure('position', [541         152        1386         719]);
        plot(LED.M565L3.WAVELENGTH, LED.M565L3.NORM_INTENSITY); hold on
        plot(System.M565L3_system.WAVELENGTH, System.M565L3_system.NORM_INTENSITY); hold on
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350 700])
        ylim([0 1.01])
        
        l = legend('M565L3 - datasheet', 'M565L3 - measured');
        set(l, 'Fontsize', 12, 'Box', 'off')
        
        
%% -- plot --
    % LED stimulus aparatus
    f = figure('position', [541         152        1386         719]);
        area(LED.M565L3.WAVELENGTH, LED.M565L3.NORM_INTENSITY, 'FaceColor', 'g', 'EdgeColor', 'none'); hold on
        area(LED.M455L3.WAVELENGTH, LED.M455L3.NORM_INTENSITY, 'FaceColor', 'b', 'EdgeColor', 'none'); hold on
    
        plot(Optics.T510lpxrxt.WAVELENGTH, Optics.T510lpxrxt.NORM_TRANSMISSION, 'r') % dichroic
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('LED stimulus')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350 700])
        ylim([0 1.01])
        
        l = legend('565 - measured', '455 - measured', 'T510lpxrxt');
        set(l, 'Fontsize', 12, 'Box', 'off')
        
        
%% -- Plot --
    % compare 455 stimulus led spectra without filter and with filter 460/14
    f = figure('position', [541         152        1386         719]);
    
        plot(System.M455L3_noFilter_system.WAVELENGTH, System.M455L3_noFilter_system.NORM_INTENSITY); hold on
        plot(System.M455L3_460_14_system.WAVELENGTH, System.M455L3_460_14_system.NORM_INTENSITY); hold on
    
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_rods ), hold on
        plot(Optogenetics.WAVELENGTH, Optogenetics.ChR2), hold on
        
        plot(Filter.FF01_460_14.WAVELENGTH, Filter.FF01_460_14.TRANSMISSION, 'color', colorHex2RGB('#a200ff'), 'LineStyle', '--', 'LineWidth', 0.5); hold on

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('455 LED with and without filter')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350 700])
        ylim([0 1.01])
        
        l = legend('455 - no filter', '455 + FF01-460/14', 'mcone', 'scone', 'rods', 'ChR2', 'FF01-460/14');
        set(l, 'Fontsize', 12, 'Box', 'off')
        
%% -- Plot --
    % evaluate 455LED filter
    f = figure('position', [541         152        1386         719]);
        plot(System.M455L3_noFilter_system.WAVELENGTH, System.M455L3_noFilter_system.NORM_INTENSITY); hold on
        
        plot(wavelengths, mouse_scone), hold on
        plot(wavelengths, chr2), hold on
        
        plot(Filter.FF01_460_14.WAVELENGTH, Filter.FF01_460_14.TRANSMISSION); hold on
        % plot(Filter.FF01_469_35.WAVELENGTH, Filter.FF01_469_35.TRANSMISSION); hold on
        % plot(Filter.FF01_470_28.WAVELENGTH, Filter.FF01_470_28.TRANSMISSION); hold on
        plot(Filter.FF01_475_42.WAVELENGTH, Filter.FF01_475_42.TRANSMISSION); hold on

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('Evaluate 455LED filters')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350 700])
        ylim([0 1.01])
        
        l = legend('455 - no filter', 'scone', 'ChR2', 'FF01-460/14', 'FF01-475/42');
        set(l, 'Fontsize', 12, 'Box', 'off')

        % calculaiton of 455LED overlap with scone curve with and without filters
        % no filter
            scone_curve = StandardTemplate(360, System.M455L3_noFilter_system.WAVELENGTH);
            M455L3_scone_nofilter = scone_curve .* System.M455L3_noFilter_system.NORM_INTENSITY;
            M455L3_scone_nofilter_int = trapz(M455L3_scone_nofilter); % integral under curve
        
        % Filter.FF01_460_14
            FF01_460_14_transmission = interp1(Filter.FF01_460_14.WAVELENGTH, Filter.FF01_460_14.TRANSMISSION, System.M455L3_noFilter_system.WAVELENGTH);
            M455L3_scone_FF01_460_14 = scone_curve .* FF01_460_14_transmission .* System.M455L3_noFilter_system.NORM_INTENSITY;
            M455L3_scone_FF01_460_14(isnan(M455L3_scone_FF01_460_14)) = 0;
            M455L3_scone_FF01_460_14_int = trapz(M455L3_scone_FF01_460_14); % integral under curve
            
            M455L3_scone_FF01_460_14_nofilter_ratio = M455L3_scone_FF01_460_14_int / M455L3_scone_nofilter_int;
           
        % Filter.FF01_460_14
            FF01_475_42_transmission = interp1(Filter.FF01_475_42.WAVELENGTH, Filter.FF01_475_42.TRANSMISSION, System.M455L3_noFilter_system.WAVELENGTH);
            M455L3_scone_FF01_475_42 = scone_curve .* FF01_475_42_transmission .* System.M455L3_noFilter_system.NORM_INTENSITY;
            M455L3_scone_FF01_475_42_int = trapz(M455L3_scone_FF01_475_42); % integral under curve     
            
            M455L3_scone_FF01_475_42_nofilter_ratio = M455L3_scone_FF01_475_42_int / M455L3_scone_nofilter_int;
        
            % values:
            fprintf('integral - NOFILTER = %g\n', M455L3_scone_nofilter_int)
            fprintf('integral - FF01-460/14 = %g; relative to NOFILTER = %0.2g%%\n', M455L3_scone_FF01_460_14_int, M455L3_scone_FF01_460_14_nofilter_ratio * 100)
            fprintf('integral - FF01-475/42 = %g; relative to NOFILTER = %0.2g%%\n', M455L3_scone_FF01_475_42_int, M455L3_scone_FF01_475_42_nofilter_ratio * 100)
            
%% -- Plot --
% Compute equvalence of 455 - 565
    % plot spectra
    f = figure('position', [541         152        1386         719]);
        % 365 LED no filter
        area(System.M365L2_system.WAVELENGTH, System.M365L2_system.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#a200ff'), 'EdgeColor', 'none'); hold on
    
        % 455 LED with 475/42 filter
        area(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#728ef4'), 'EdgeColor', 'none'); hold on
        plot(System.M455L3_475_42_455_10.WAVELENGTH, System.M455L3_475_42_455_10.NORM_INTENSITY, 'LineWidth', 1.5), hold on
        
        % 565 LED with 571/72 filter
        area(System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#f9d59c'), 'EdgeColor', 'none'); hold on
        plot(System.M565L3_571_72_568_10.WAVELENGTH, System.M565L3_571_72_568_10.NORM_INTENSITY, 'LineWidth', 1.5), hold on
        
        % photoreceptors
        plot(Photoreceptor.WAVELENGTH, mouse_mcone, 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, mouse_scone, 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, mouse_rods , 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, mouse_melanopsin , 'LineWidth', 1.5), hold on
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Calculating equvalence 455 - 565')
        xlim([350 700])
        ylim([0  1.01])

        l = legend('LED 365 - measured', 'LED 455 - 475/42 - measured', 'LED 455 - 475/42 - 455/10 - measured', 'LED 565 - 571/72 - measured', 'LED 565 - 571/72 - 568/10', 'mouse m-cone', 'mouse s-cone', 'mouse rods', 'mouse melanosin');
               
        set(l, 'Fontsize', 12, 'Box', 'off')     
        
    % plot power conversion curves
    f = figure('position', [1229          24 1260         937]);
        subplot(2, 2, 1)
            plot(PowerCurves.Power_400_20150615.VOLTAGE, PowerCurves.Power_400_20150615.POWER, 'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'o'); hold on
            plot(PowerCurves.Power_455_475_42_20150615.VOLTAGE, PowerCurves.Power_455_475_42_20150615.POWER, 'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'o'); hold on
            plot(PowerCurves.Power_565_571_72_20150615.VOLTAGE, PowerCurves.Power_565_571_72_20150615.POWER, 'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'o'); hold on
            
            box off
            set(gca, 'tickdir', 'out', 'FontSize', 12)
            xlabel('Voltage (V)')
            ylabel('Power (uW)')
            xlim([0 5])
            
            l = legend('400', '455 - 475/42', '565 - 571/72', 'location', 'northwest');
            
            set(l, 'Fontsize', 12, 'Box', 'off')
            
        subplot(2, 2, 2)
            plot(PowerCurves.Power_455_475_42_20150615.POWER, PowerCurves.Power_455_475_42_455_10_20150615.POWER, 'LineWidth', 1.5, 'LineStyle', 'none', 'Marker', 'o'); hold on
            
            box off
            set(gca, 'tickdir', 'out', 'FontSize', 12)
            xlabel('Power (uW) - 455 - 475/42')
            ylabel('Power (uW) - 455 - 475/42 - 455/10')
            
            l = legend('455 - 475/42 vs 455 - 475/42 - 455/10', 'location', 'northwest');
            
            set(l, 'Fontsize', 12, 'Box', 'off')            
        
        subplot(2, 2, 3)
            plot(PowerCurves.Power_565_571_72_20150615.POWER, PowerCurves.Power_568_571_72_568_10_20150615.POWER, 'LineWidth', 1.5, 'LineStyle', 'none', 'Marker', 'o'); hold on
            
            box off
            set(gca, 'tickdir', 'out', 'FontSize', 12)
            xlabel('Power (uW) - 565 - 571/72')
            ylabel('Power (uW) - 565 - 571/72 - 568/10')
            
            l = legend('565 - 571/72 vs 455 - 571/72 - 568/10', 'location', 'northwest');
            
            set(l, 'Fontsize', 12, 'Box', 'off')  
            
            
    % compute light equivalence
        % compute scale factor for measured part of spectrum to complete spectrum
            integral.M455L3_475_42 = trapz(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY);
            integral.M455L3_475_42_455_10 = trapz(System.M455L3_475_42_455_10.WAVELENGTH, System.M455L3_475_42_455_10.NORM_INTENSITY);
            
            M455_measured_ratio = integral.M455L3_475_42_455_10 / integral.M455L3_475_42; % proportion of LED curve that is measured
            M455_scale = 1/M455_measured_ratio;

            integral.M565L3_571_72 = trapz(System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY);
            integral.System.M565L3_571_72_568_10 = trapz(System.M565L3_571_72_568_10.WAVELENGTH, System.M565L3_571_72_568_10.NORM_INTENSITY);
            
            M565_measured_ratio = integral.System.M565L3_571_72_568_10 / integral.M565L3_571_72; % proportion of LED curve that is measured
            M565_scale = 1/M565_measured_ratio;
        
            
        % product with mcone curve
            M455_mcone = crossProduct(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY, Photoreceptor.WAVELENGTH, mouse_mcone);
            M565_mcone = crossProduct(System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY, Photoreceptor.WAVELENGTH, mouse_mcone);
            
            mcone_M455_M565 = M455_mcone.INTEGRAL / M565_mcone.INTEGRAL; % ratio of M cone stimulation from 455 vs 565
            
            
        % convert equivalence for 565 to 455 for mcone, with spectrum correction
            M455_power = 0:100; % uw (power of 455 that i want to calculate 565 equivalent)
            
            for qq = 1:length(M455_power)
                M455_455_power = interp1(PowerCurves.Power_455_475_42_20150615.POWER, PowerCurves.Power_455_475_42_455_10_20150615.POWER, M455_power(qq));
                M455_photons = power2density(M455_455_power, 455); % <== ENTER POWER MEAUSREMENT HERE
                M455_photons_scaled = M455_photons * M455_scale; % quanta from 

                M565_photons_scaled = M455_photons_scaled * mcone_M455_M565; % convert to equiv 565 quanta

                M565_photons = M565_photons_scaled / M565_scale;
                M565_568_power = density2power(M565_photons, 568); % <== ENTER POWER MEAUSREMENT HERE
            
                M565_power(qq) = interp1(PowerCurves.Power_568_571_72_568_10_20150615.POWER, PowerCurves.Power_565_571_72_20150615.POWER, M565_568_power);
            end
            
            mcone_M455_M565_corrected = nanmean(M565_power ./ M455_power);
            
        subplot(2, 2, 4)
            plot(M455_power, M565_power, 'LineWidth', 1.5, 'LineStyle', 'none', 'Marker', 'o'); hold on
            
            text(0.01, 0.8, sprintf('spectra corrected 455:565 ratio = 1:%0.2f', mcone_M455_M565_corrected), 'units', 'normalized', 'fontsize', 12)
            
            box off
            set(gca, 'tickdir', 'out', 'FontSize', 12)
            xlabel('Power (uW) - 455 - 475/42')
            ylabel('Power (uW) - 565 - 571/72')
            
            l = legend('455 - 475/42 vs 565 - 571/72', 'location', 'northwest');
            
            set(l, 'Fontsize', 12, 'Box', 'off')  
        
%% -- plot --
    % GCaMP5g vs R-GECO1
    f = figure('position', [541         152        1386         719]);
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'color', colorHex2RGB('#66cccb')); hold on
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, 'color', colorHex2RGB('#66cccb'), 'Linestyle', '--'); hold on
        
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EXCITATION, 'Linestyle', ':', 'lineWidth', 2)
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EMISSION, 'Linestyle', '--', 'lineWidth', 2)   

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('GCaMP RCaMP')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])

        l = legend('GCaMP5g excitation', 'GCaMP5g emission', 'R-GECO1 excitation', 'R-GECO1 emission', 'location', 'northwest');
        set(l, 'Fontsize', 12, 'Box', 'off')     
    
%% -- Plot --
% compare RCaMP with Alexa 568
    f = figure('position', [541         152        1386         719]);
        plot(Flurophore.RCaMP1e.WAVELENGTH,  Flurophore.RCaMP1e.NORM_EMISSION); hold on
        plot(Flurophore.RCaMP1e.WAVELENGTH,  Flurophore.RCaMP1e.NORM_EXCITATION); hold on
        
        plot(Flurophore.Alexa568.WAVELENGTH, Flurophore.Alexa568.NORM_EXCITATION, '--');
        plot(Flurophore.Alexa568.WAVELENGTH, Flurophore.Alexa568.NORM_EMISSION,   '--');
        
        plot([561 561], [0 1], 'color', colorHex2RGB('#4af406'), 'LineStyle', '--', 'LineWidth', 1.5)

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('Comparing RCaMP1e with Alexa568')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('RCaMP1e emission', 'RCaMP1e excitation', 'Alexa568 excitation', 'Alexa568 emission');
        set(l, 'Fontsize', 12, 'Box', 'off')
   
%% -- plot --
% RCaMP mouse experiment
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % LEDs
        area(System.M365L2_system.WAVELENGTH, System.M365L2_system.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#a200ff'), 'EdgeColor', 'none'); hold on
        area(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#728ef4'), 'EdgeColor', 'none'); hold on
    
        % calcium indicator spectra
        plot(Flurophore.RCaMP1e.WAVELENGTH, Flurophore.RCaMP1e.NORM_EMISSION, '--'); hold on
        plot(Flurophore.RCaMP1e.WAVELENGTH, Flurophore.RCaMP1e.NORM_EXCITATION, '--'); hold on        
    
        % photoreceptor sensitivites
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone, 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone, 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_rods , 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin, 'LineWidth', 1.5); hold on
        
        % channelrhodopsins
        plot(Optogenetics.WAVELENGTH,  Optogenetics.ChR2, 'LineWidth', 1.5, 'Linestyle', '-.'), hold on
        
        % barrier filter
        plot(Filter.FF01_630_92.WAVELENGTH, Filter.FF01_630_92.TRANSMISSION), hold on
        
        % fluorescence excitation laser
        plot(System.L561_system.WAVELENGTH, System.L561_system.NORM_INTENSITY, 'color', colorHex2RGB('#4af406'), 'LineStyle', '--', 'LineWidth', 1);
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('RCaMP mouse experiment')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('365 - measured', '455(475/42) - measured', 'RCaMP1e emssion', 'RCaMP1e excitation', 'mouse M-cone', 'mouse S-cone', 'mouse rod', 'mouse melanopsin', 'ChR2', ...
            'FF01-630/92 - emission filter', '561 nm imaging laser', 'location', 'northeast');
        set(l, 'Fontsize', 12, 'Box', 'off') 
        
    save_name = 'RCaMP mouse experiment.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)
    
%% -- plot --
% notes: 
% - GFP emission overlap with barrir filter is negligible 
% - M490L3 LED is more ideal to activate Chronos, but has much more overlap with jRGECO1a excitation and may cause response artefact.

% JRGECO1a mouse experiment
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % LEDs
        area(System.M365L2_system.WAVELENGTH, System.M365L2_system.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#a200ff'), 'EdgeColor', 'none'); hold on
        area(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#728ef4'), 'EdgeColor', 'none'); hold on
        area(LED.M490L3.WAVELENGTH, LED.M490L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#81def2'), 'EdgeColor', 'none'); hold on
        
        % calcium indicator spectra
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EMISSION, '--'); hold on
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EXCITATION, '--'); hold on        
    
        % photoreceptor sensitivites
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone, 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone, 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_rods , 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin, 'LineWidth', 1.5); hold on
        
        % channelrhodopsins
%         plot(Optogenetics.WAVELENGTH,  Optogenetics.ChR2, 'LineWidth', 1.5, 'Linestyle', '-.'), hold on
        plot(Biology.chronos.WAVELENGTH, Biology.chronos.NORM_RESPONSE, '-.'); hold on
%         plot(Flurophore.GFP.WAVELENGTH, Flurophore.GFP.NORM_EXCITATION, '-.'); hold on
%         plot(Flurophore.GFP.WAVELENGTH, Flurophore.GFP.NORM_EMISSION, '-.'); hold on
        
        % barrier filter
        plot(Filter.FF01_630_92.WAVELENGTH, Filter.FF01_630_92.TRANSMISSION), hold on
        
        % fluorescence excitation laser
        plot(System.L561_system.WAVELENGTH, System.L561_system.NORM_INTENSITY, 'color', colorHex2RGB('#4af406'), 'LineStyle', '--', 'LineWidth', 1);
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('JRGECO1a mouse experiment')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('365 - measured', '455(475/42) - measured', '490 LED',...
            'R-GECO1 emssion', 'R-GECO1 excitation',...
            'mouse M-cone', 'mouse S-cone', 'mouse rod', 'mouse melanopsin', 'chronos',...
            'FF01-630/92 - emission filter', '561 nm imaging laser', 'location', 'northeast');
        set(l, 'Fontsize', 12, 'Box', 'off') 
        
    save_name = 'JRGECO1a mouse experiment.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)        

%% -- plot --
% JRGECO1a-2P mouse experiment
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % LEDs
        area(System.M365L2_system.WAVELENGTH, System.M365L2_system.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#a200ff'), 'EdgeColor', 'none'); hold on
        area(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#728ef4'), 'EdgeColor', 'none'); hold on
    
        % calcium indicator spectra
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EXCITATION, '--'); hold on   
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EXCITATION_2P, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EMISSION, '--'); hold on
    
        % photoreceptor sensitivites
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone, 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone, 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_rods , 'LineWidth', 1.5), hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin, 'LineWidth', 1.5); hold on
        
        % channelrhodopsins
%         plot(Optogenetics.WAVELENGTH,  Optogenetics.ChR2, 'LineWidth', 1.5, 'Linestyle', '-.'), hold on
        plot(Biology.chronos.WAVELENGTH, Biology.chronos.NORM_RESPONSE, '-.'); hold on
        
        % barrier filter
        plot(Filter.FF01_630_92.WAVELENGTH, Filter.FF01_630_92.TRANSMISSION), hold on
        
        % fluorescence excitation laser
        plot(System.L561_system.WAVELENGTH, System.L561_system.NORM_INTENSITY, 'color', colorHex2RGB('#4af406'), 'LineStyle', '--', 'LineWidth', 1);
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('JRGECO1a 2P mouse experiment')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  1100])
        ylim([0   1.01])
        
        l = legend('365 - measured', '455(475/42) - measured', 'R-GECO1 excitation', 'R-GECO1 2P excitation', 'R-GECO1 emssion', ...
            'mouse M-cone', 'mouse S-cone', 'mouse rod', 'mouse melanopsin', 'chronos', ...
            'FF01-630/92 - emission filter', '561 nm imaging laser', 'location', 'northeast');
        set(l, 'Fontsize', 12, 'Box', 'off') 
        
    save_name = 'JRGECO1a 2P mouse experiment.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)       
    
%% -- Plot --
% GCaMP mouse experiment (with chrimson)
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % LEDs
        area(System.M365L2_system.WAVELENGTH, System.M365L2_system.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#a200ff'), 'EdgeColor', 'none'); hold on
        area(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#728ef4'), 'EdgeColor', 'none'); hold on
%         area(System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#f9d59c'), 'EdgeColor', 'none'); hold on
        area(System.M595L3_BLP_594R.WAVELENGTH, System.M595L3_BLP_594R.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#ff9999'), 'EdgeColor', 'none'); hold on


        % fluorescence excitation laser
        plot(System.L488_system.WAVELENGTH, System.L488_system.NORM_INTENSITY, 'color', colorHex2RGB('#00f7ff'), 'LineStyle', '-.', 'LineWidth', 1);
    
        % photoreceptor sensitivites
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_rods , 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin, 'LineWidth', 1.5); hold on
        
        % channelrhodopsins
        plot(Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE, '-.'); hold on
%         plot(Biology.C1V1.WAVELENGTH, Biology.C1V1.NORM_RESPONSE, '-.'); hold on
        
        % calcium indicator spectra
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION,   'color', colorHex2RGB('#00aba9'), 'LineStyle', ':'); hold on
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'color', colorHex2RGB('#66cccb'), 'LineStyle', ':'); hold on
        
        % barrier filter
        plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION, 'color', colorHex2RGB('#36ff00'), 'LineStyle', '-.'), hold on
        

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('GCaMP mouse experiment')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('365 - measured', '455 - 475/42 - measured', '565 - 571/72 - measured', '488 laser', ...
             'mcone', 'scone', 'rod', 'melanopsin', 'chrimson', 'GCaMP5g - excitation', 'GCaMP5g - emission', 'barrier - 520/35');
        
        set(l, 'Fontsize', 12, 'Box', 'off') 
        
    save_name = 'GCaMP (Chrimson) mouse experiment.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)
    
%% -- light power in daylight --
% put stimulation power in context of daylight spectral irradiance

    % convert LED power to irradiance 
    LED_FOV_r    = 8.6/2 * 33;        % convert deg to µm
    FOV_LED_um   = pi * LED_FOV_r^2;  % µm^2
    FOV_LED_cm   = FOV_LED_um / 10^8; % cm^2
    
    led_365_power = 20;  % < enter value here
    led_620_power = 100; % < enter value here
    
    led_365_I = led_365_power / FOV_LED_cm;
    led_620_I = led_620_power / FOV_LED_cm;
    
    solar_I_umcm = Other.SolarIrradianceSeaLevel.IRRADIANCE .* 10^6 ./ 10^4; % convert to µW/cm^2
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);

        plot(Other.SolarIrradianceSeaLevel.WAVELENGTH, solar_I_umcm, 'k-'); hold on
        plot(365, led_365_I, 'bo')
        plot(620, led_620_I, 'rs')

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12, 'yscale', 'log')
        title('Irrandiance for daylight and stimuli')
        xlabel('wavelength (nm)')
        ylabel('Irradiance (µW/cm^2/nm^1)')
        xlim([350  700])
%         ylim([0   1.01]) 
        legend('Normal incident daylight at sea level on a clear day', '365 nm LED', '620 nm LED')
        
        save_name = 'Irradiance of daylight and stimuli.pdf';
        save_dir  = fullfile('.', 'Figures');
        save_path = fullfile(save_dir, save_name);
        exportpdf(f, save_path)


    % plot solar irradiance in terms of LED power as in our setup
    solar_power = solar_I_umcm .* FOV_LED_cm;
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);

        plot(Other.SolarIrradianceSeaLevel.WAVELENGTH, solar_power, 'k-'); hold on

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('Irrandiance for daylight and stimuli')
        xlabel('wavelength (nm)')
        ylabel('power (µW/nm^1)')
        xlim([350  700])  
        
        save_name = 'Irradiance of daylight as LED power.pdf';
        save_dir  = fullfile('.', 'Figures');
        save_path = fullfile(save_dir, save_name);
        exportpdf(f, save_path)
    

%% -- Plot --
% LED stimuli with peak wavelengths
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % LEDs
        area(System.M365L2_system.WAVELENGTH, smooth(System.M365L2_system.NORM_INTENSITY, 20), 'FaceColor', colorHex2RGB('#a200ff'), 'EdgeColor', 'none'); hold on
        area(System.M455L3_475_42.WAVELENGTH, smooth(System.M455L3_475_42.NORM_INTENSITY, 20), 'FaceColor', colorHex2RGB('#728ef4'), 'EdgeColor', 'none'); hold on
        area(System.M565L3_571_72.WAVELENGTH, smooth(System.M565L3_571_72.NORM_INTENSITY, 20), 'FaceColor', colorHex2RGB('#f9d59c'), 'EdgeColor', 'none'); hold on
        area(System.M595L3_BLP_594R.WAVELENGTH, smooth(System.M595L3_BLP_594R.NORM_INTENSITY, 20), 'FaceColor', colorHex2RGB('#ff9999'), 'EdgeColor', 'none'); hold on
    
        [~, M365_max_ind] = max(smooth(System.M365L2_system.NORM_INTENSITY, 20));
        [~, M455_max_ind] = max(smooth(System.M455L3_475_42.NORM_INTENSITY, 20));
        [~, M565_max_ind] = max(smooth(System.M565L3_571_72.NORM_INTENSITY, 20));
        [~, M595_max_ind] = max(smooth(System.M595L3_BLP_594R.NORM_INTENSITY, 20));
        
        text(0.5, 0.95, sprintf('max wavelengths: %03.0f, %03.0f, %03.0f, %03.0f', System.M365L2_system.WAVELENGTH(M365_max_ind),...
            System.M455L3_475_42.WAVELENGTH(M455_max_ind), System.M565L3_571_72.WAVELENGTH(M565_max_ind), System.M595L3_BLP_594R.WAVELENGTH(M595_max_ind)),...
            'units', 'normalized', 'fontsize', 12, 'horizontalalignment', 'center')
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('LED stimuli 1P mouse')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('365 - measured', '455 - 475/42 - measured', '565 - 571/72 - measured', 'LED 595 - BLP01-594R - measured');
        
        set(l, 'Fontsize', 12, 'Box', 'off') 
        
        
        save_name = '1P mouse LED stimuli.pdf';
        save_dir  = fullfile('.', 'Figures');
        save_path = fullfile(save_dir, save_name);
        exportpdf(f, save_path)

%% -- Plot --
% mouse GCaMP.jaws experiment

f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % LEDs
        area(System.M595L3_BLP_594R.WAVELENGTH, System.M595L3_BLP_594R.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#ff9999'), 'EdgeColor', 'none'); hold on

        % fluorescence excitation laser
        plot(System.L488_system.WAVELENGTH, System.L488_system.NORM_INTENSITY, 'color', colorHex2RGB('#00f7ff'), 'LineStyle', '-.', 'LineWidth', 1);
    
        % photoreceptor sensitivites
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_rods , 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin, 'LineWidth', 1.5); hold on
        
        % channelrhodopsins
        plot(Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE, '-.'); hold on
        plot(Biology.jaws.WAVELENGTH, Biology.jaws.NORM_RESPONSE, '-.'); hold on
        
        % calcium indicator spectra
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION,   'color', colorHex2RGB('#00aba9'), 'LineStyle', ':'); hold on
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'color', colorHex2RGB('#66cccb'), 'LineStyle', ':'); hold on
        
        % barrier filter
        plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION, 'color', colorHex2RGB('#36ff00'), 'LineStyle', '-.'), hold on
        


        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('GCaMP-Chrimson mouse experiment')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('LED 595 - BLP01-594R - measured', '488 laser', ...
             'mcone', 'scone', 'rod', 'melanopsin',...
             'chrimson', 'jaws', 'GCaMP5g - excitation', 'GCaMP5g - emission',...
             'barrier - 520/35');
        
        set(l, 'Fontsize', 12, 'Box', 'off') 

    save_name = 'GCaMP-jaws mouse experiment.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)        
        
%% - plot sensitivity curve of jaws - 
% Fig 2e - Chuong etal 2014 Nat Neurosci

    data = [15.9819819819820	21.5116279069773
            16.5225225225225	105.232558139535
            17.0720720720721	134.883720930233
            18.0000000000000	197.674418604651];
        
    data(:,1) = 10.^data(:,1); % convert from log scale
    
    % reproduce figure from paper
    f = figure('color', [1 1 1]);
        semilogx(data(:,1), data(:,2), 'o-')
        
        ylim([0 250])
        xlim([10^15, 10^19])
        set(gca, 'box', 'off', 'tickdir', 'out', 'FontSize', 12)
        xlabel('600 nm light intensity (photons/cm^2/s)')
        ylabel('mean ganglion cell spiking rate (Hz)')
        
    % convert light power units
    power = density2power(data(:,1), 600); % power in µW
    
    f = figure('color', [1 1 1]);
        semilogx(power, data(:,2), 'o-')
        
        ylim([0 250])
%         xlim([10^15, 10^19])
        set(gca, 'box', 'off', 'tickdir', 'out', 'FontSize', 12)
        xlabel('600 nm light intensity (µW/cm^2)')
        ylabel('mean ganglion cell spiking rate (Hz)')    
    
        
    % convert to LED power with our setup
    LED_FOV_r    = 8.6/2 * 33;        % convert deg to µm
    FOV_LED_um   = pi * LED_FOV_r^2;  % µm^2
    FOV_LED_cm   = FOV_LED_um / 10^8; % cm^2
    
    raw_power = power * FOV_LED_cm;
    
    f = figure('color', [1 1 1]);
        semilogx(raw_power, data(:,2), 'o-')
        
        ylim([0 250])
%         xlim([10^15, 10^19])
        set(gca, 'box', 'off', 'tickdir', 'out', 'FontSize', 12)
        xlabel('600 nm light intensity (µW over 8.6 deg field)')
        ylabel('mean ganglion cell spiking rate (Hz)')    
        title('Sensitivity of jaws to 600 nm light')
        
        save_name = 'Sensitivity of jaws to 600 nm light.pdf';
        save_dir  = fullfile('.', 'Figures');
        save_path = fullfile(save_dir, save_name);
        exportpdf(f, save_path)
    
%% -- Plot --
% evaluate better LED for mouse GCaMP.Chrimson experiment
% Conclusion:
%   LED: 595 nm 
%   LED filter: BLP01_594R

f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % LEDs
%         area(System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#f9d59c'), 'EdgeColor', 'none'); hold on
        area(System.M595L3_BLP_594R.WAVELENGTH, System.M595L3_BLP_594R.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#ff9999'), 'EdgeColor', 'none'); hold on

%         % alternative LEDs with filter % theoretical
%         M595L3_BLP01_594R = crossProduct(LED.M595L3.WAVELENGTH, LED.M595L3.NORM_INTENSITY, Filter.BLP01_594R.WAVELENGTH, Filter.BLP01_594R.TRANSMISSION);
%         area(M595L3_BLP01_594R.WAVELENGTH, M595L3_BLP01_594R.SPECTRA, 'FaceColor', colorHex2RGB('#ff9999'), 'EdgeColor', 'none'); hold on

        % fluorescence excitation laser
        plot(System.L488_system.WAVELENGTH, System.L488_system.NORM_INTENSITY, 'color', colorHex2RGB('#00f7ff'), 'LineStyle', '-.', 'LineWidth', 1);
    
        % photoreceptor sensitivites
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_rods , 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin, 'LineWidth', 1.5); hold on
        
        % channelrhodopsins
        plot(Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE, '-.'); hold on
        
        % calcium indicator spectra
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION,   'color', colorHex2RGB('#00aba9'), 'LineStyle', ':'); hold on
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'color', colorHex2RGB('#66cccb'), 'LineStyle', ':'); hold on
        
        % barrier filter
        plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION, 'color', colorHex2RGB('#36ff00'), 'LineStyle', '-.'), hold on
        


        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('GCaMP-Chrimson mouse experiment')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('LED 595 - BLP01-594R - measured', '488 laser', ...
             'mcone', 'scone', 'rod', 'melanopsin',...
             'chrimson', 'GCaMP5g - excitation', 'GCaMP5g - emission',...
             'barrier - 520/35');
        
        set(l, 'Fontsize', 12, 'Box', 'off') 

    save_name = 'GCaMP-Chrimson mouse experiment - evaulation.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)
    
    % compute integrals for LEDS and Chrimson
    cutoff = 550; % cuttoff frequency for red led filter
    
    % first convolve LED with filter
    LED590_filter = crossProduct(LED.M590L3.WAVELENGTH, LED.M590L3.NORM_INTENSITY, Filter.BLP01_594R.WAVELENGTH, Filter.BLP01_594R.TRANSMISSION);
    LED595_filter = crossProduct(LED.M595L3.WAVELENGTH, LED.M595L3.NORM_INTENSITY, Filter.BLP01_594R.WAVELENGTH, Filter.BLP01_594R.TRANSMISSION);
    LED617_filter = crossProduct(LED.M617L3.WAVELENGTH, LED.M617L3.NORM_INTENSITY, Filter.BLP01_594R.WAVELENGTH, Filter.BLP01_594R.TRANSMISSION);
    LED625_filter = crossProduct(LED.M625L3.WAVELENGTH, LED.M625L3.NORM_INTENSITY, Filter.BLP01_594R.WAVELENGTH, Filter.BLP01_594R.TRANSMISSION);
    
    
%     LED590_filter = crossProduct(LED.M590L3.WAVELENGTH, LED.M590L3.NORM_INTENSITY, Filter.FF01_593_LP.WAVELENGTH, Filter.FF01_593_LP.TRANSMISSION);
%     LED595_filter = crossProduct(LED.M595L3.WAVELENGTH, LED.M595L3.NORM_INTENSITY, Filter.FF01_593_LP.WAVELENGTH, Filter.FF01_593_LP.TRANSMISSION);
%     LED617_filter = crossProduct(LED.M617L3.WAVELENGTH, LED.M617L3.NORM_INTENSITY, Filter.FF01_593_LP.WAVELENGTH, Filter.FF01_593_LP.TRANSMISSION);
%     LED625_filter = crossProduct(LED.M625L3.WAVELENGTH, LED.M625L3.NORM_INTENSITY, Filter.FF01_593_LP.WAVELENGTH, Filter.FF01_593_LP.TRANSMISSION);
    
    % LED and Chrimson
    chrimson_590 = crossProduct(LED590_filter.WAVELENGTH, LED590_filter.SPECTRA, Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE);
    chrimson_595 = crossProduct(LED595_filter.WAVELENGTH, LED595_filter.SPECTRA, Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE);
    chrimson_617 = crossProduct(LED617_filter.WAVELENGTH, LED617_filter.SPECTRA, Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE);
    chrimson_625 = crossProduct(LED625_filter.WAVELENGTH, LED625_filter.SPECTRA, Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE);
    
    % LED and M-opsin
    mopsin_590 = crossProduct(LED590_filter.WAVELENGTH, LED590_filter.SPECTRA, Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone);
    mopsin_595 = crossProduct(LED595_filter.WAVELENGTH, LED595_filter.SPECTRA, Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone);
    mopsin_617 = crossProduct(LED617_filter.WAVELENGTH, LED617_filter.SPECTRA, Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone);
    mopsin_625 = crossProduct(LED625_filter.WAVELENGTH, LED625_filter.SPECTRA, Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone);
    
    % LED and Melanopsin
    melanopsin_590 = crossProduct(LED590_filter.WAVELENGTH, LED590_filter.SPECTRA, Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin);
    melanopsin_595 = crossProduct(LED595_filter.WAVELENGTH, LED595_filter.SPECTRA, Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin);
    melanopsin_617 = crossProduct(LED617_filter.WAVELENGTH, LED617_filter.SPECTRA, Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin);
    melanopsin_625 = crossProduct(LED625_filter.WAVELENGTH, LED625_filter.SPECTRA, Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin);
    
%     figure
%         plot(chrimson_590.WAVELENGTH, chrimson_590.SPECTRA); hold on
%         plot(chrimson_595.WAVELENGTH, chrimson_595.SPECTRA); hold on
%         plot(chrimson_617.WAVELENGTH, chrimson_617.SPECTRA); hold on
%         plot(chrimson_625.WAVELENGTH, chrimson_625.SPECTRA); hold on
%         
%         box off
%         set(gca, 'tickdir', 'out', 'FontSize', 12)
%         xlabel('wavelength (nm)')
%         ylabel('LED-Chrimson spectra product')
%         xlim([cutoff 700])
%     
%         l = legend('LED 590', 'LED 595', 'LED 617', 'LED 625');
        
    % Compute vector of intergrals at "typical LED output power":
        % for reference the current 565 nm LED has a typical power rating of 979 mW
    chrimson_ints   = [chrimson_590.INTEGRAL * 170,   chrimson_595.INTEGRAL * 502,   chrimson_617.INTEGRAL * 650,   chrimson_625.INTEGRAL * 770];
    mopsin_ints     = [mopsin_590.INTEGRAL * 170,     mopsin_595.INTEGRAL * 502,     mopsin_617.INTEGRAL * 650,     mopsin_625.INTEGRAL * 770];
    melanopsin_ints = [melanopsin_590.INTEGRAL * 170, melanopsin_595.INTEGRAL * 502, melanopsin_617.INTEGRAL * 650, melanopsin_625.INTEGRAL * 770];
    
    figure
        bar([chrimson_ints; mopsin_ints; melanopsin_ints])
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12, 'xticklabel', {'Chrimson', 'M-opsin', 'melanopsin'})
%         title(sprintf('cutoff = %d nm', cutoff))
        xlabel('opsin')
        ylabel('LED-opsin integral at typical power rating')
        
        legend('LED 590', 'LED 595', 'LED 617', 'LED 625')
        
        
        
        
        
%% -- PLOT --
% compute effective stimulus contrast for GCaMP imaging

% TODO: make legend labels dynamic, so that i can change the power values and observe the plot chages

% need to account for area
% LED FOV on retina is 8.6 deg circular
% imaging FOV 5x6.7 deg rectangular field

    % Calcualte FOV of lights as area
        % 1 deg visual angle in mouse = 33 µm (Geng et al. 2012. Biomed Optics Exp)

        LED_FOV_r = 8.6/2 * 33;       % convert to µm
        FOV.LED   = pi * LED_FOV_r^2; % µm^2

        imaging_FOV_x = 6.7 * 33; % convert to µm
        imaging_FOV_y = 5 * 33;   % convert to µm
        FOV.imaging   = imaging_FOV_x * imaging_FOV_y; % µm^2

        
    % light powers:
    power_365 = 20;  % (uW)
    power_565 = 100; % (uW)
    power_488 = 103; % (uW)

    % S-cone
    scone_365 = crossProduct(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone, System.M365L2_system.WAVELENGTH, System.M365L2_system.NORM_INTENSITY);
    scone_565 = crossProduct(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone, System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY);
    scone_488 = crossProduct(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone, System.L488_system.WAVELENGTH, System.L488_system.NORM_INTENSITY);

    % M-cone
    mcone_365 = crossProduct(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone, System.M365L2_system.WAVELENGTH, System.M365L2_system.NORM_INTENSITY);
    mcone_565 = crossProduct(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone, System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY);
    mcone_488 = crossProduct(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone, System.L488_system.WAVELENGTH, System.L488_system.NORM_INTENSITY);

    % calculate relative activation, normalised for FOV
    scone_activation.l365 = scone_365.INTEGRAL .* power_365 ./ FOV.LED;
    scone_activation.l565 = scone_565.INTEGRAL .* power_565 ./ FOV.LED;
    scone_activation.l455 = scone_488.INTEGRAL .* power_488 ./ FOV.imaging;
    
    mcone_activation.l365 = mcone_365.INTEGRAL .* power_365 ./ FOV.LED;
    mcone_activation.l565 = mcone_565.INTEGRAL .* power_565 ./ FOV.LED;
    mcone_activation.l455 = mcone_488.INTEGRAL .* power_488 ./ FOV.imaging;    
    
%     scone_activation.l365 = scone_365.INTEGRAL .* power_365 ;
%     scone_activation.l565 = scone_565.INTEGRAL .* power_565 ;
%     scone_activation.l455 = scone_488.INTEGRAL .* power_488 ;
%     
%     mcone_activation.l365 = mcone_365.INTEGRAL .* power_365 ;
%     mcone_activation.l565 = mcone_565.INTEGRAL .* power_565 ;
%     mcone_activation.l455 = mcone_488.INTEGRAL .* power_488 ;    
    
    % plot relative activation
    f = figure;
        bar([scone_activation.l365, scone_activation.l455, scone_activation.l565; mcone_activation.l365, mcone_activation.l455, mcone_activation.l565])
    
        set(gca, 'TickDir', 'out', 'Box', 'off', 'XTickLabel', {'S-cone', 'M-cone'})
        legend('365 LED (20 µW)', '488 laser (103 µW)', '565 LED (100 µW)', 'location', 'northwest')
        legend BOXOFF
        title('Relative activation')
        
    % plot relative activation as a ratio of scone/mcone
    f = figure;
        bar([scone_activation.l365 / mcone_activation.l365, scone_activation.l565 / mcone_activation.l565; mcone_activation.l365 / scone_activation.l365, mcone_activation.l565/ scone_activation.l565]) 
    
        set(gca, 'TickDir', 'out', 'Box', 'off', 'XTickLabel', {'S-cone/M-cone', 'M-cone/S-cone'})
        legend('365 LED (20 µW)', '565 LED (100 µW)', 'location', 'northwest')
        legend BOXOFF
        title('Relative activation')   
    
        
        
    % compute relative contrast, stimuli vs imaging light
    scone_contrast.l365 = ((scone_activation.l365 + scone_activation.l455) - scone_activation.l455) ./ ((scone_activation.l365 + scone_activation.l455) + scone_activation.l455);
    scone_contrast.l565 = ((scone_activation.l565 + scone_activation.l455) - scone_activation.l455) ./ ((scone_activation.l565 + scone_activation.l455) + scone_activation.l455);
    
    mcone_contrast.l365 = ((mcone_activation.l365 + mcone_activation.l455) - mcone_activation.l455) ./ ((mcone_activation.l365 + mcone_activation.l455) + mcone_activation.l455);
    mcone_contrast.l565 = ((mcone_activation.l565 + mcone_activation.l455) - mcone_activation.l455) ./ ((mcone_activation.l565 + mcone_activation.l455) + mcone_activation.l455);    
    
    % plot
    f = figure;
        bar([scone_contrast.l365, scone_contrast.l565; mcone_contrast.l365, mcone_contrast.l565])
    
        set(gca, 'TickDir', 'out', 'Box', 'off', 'XTickLabel', {'S-cone', 'M-cone'})
        ylim([0 1])
        legend('365 LED (20 µW)', '565 LED (100 µW)', 'location', 'north')
        legend BOXOFF
        title('Relative contrast imaging with 488 nm (103 µW)')    
        ylabel('Relative cone contrast')
        
%% -- Calculator --
% Convert absolute power to power density (power per unit area)

    % LED FOV on retina is 8.6 deg circular
    % imaging FOV 5x6.7 deg rectangular field        
    LED_FOV_r  = 8.6/2 * 33;       % convert to µm
    FOV.LED    = pi * LED_FOV_r^2; % µm^2 
    FOV.LED_cm = pi * (LED_FOV_r * 10^-4)^2; % cm^2 
    
    LED_power = [0:5:100]; % µm
    LED_power = LED_power * 10^-3;
    
    LED_powerDensity    = LED_power ./ FOV.LED;
    LED_powerDensity_cm = LED_power ./ FOV.LED_cm;
    
    
    figure
        plot(LED_power, LED_powerDensity, 'o-')
        
        set(gca, 'tickdir', 'out', 'FontSize', 12, 'box', 'off')
        xlabel('LED power (µW)')
        ylabel('LED power density (µW/µm^2)')

    figure
        plot(LED_power, LED_powerDensity_cm, 'o-')
        
        set(gca, 'tickdir', 'out', 'FontSize', 12, 'box', 'off')
        xlabel('LED power (µW)')
        ylabel('LED power density (µW/cm^2)')    
        
%% -- Calculator --
% Convert light power of imaging light 
    % compute imaging FOV
    imaging_FOV_x = 6.7 * 33; % convert to µm
    imaging_FOV_y = 5 * 33;   % convert to µm
    FOV.imaging    = imaging_FOV_x * imaging_FOV_y; % µm^2        
    FOV.imaging_cm = FOV.imaging * 10^-8; % cm^2  
        
    % light power
    L488_power = 100; % µW
    L561_power = 100; % µW
        
    % compute power per unit area (intensity)
    L488_power_intensity = L488_power * 10^-6 / FOV.imaging_cm; % J/cm^2/sec
    L561_power_intensity = L561_power * 10^-6 / FOV.imaging_cm; % J/cm^2/sec
    
    % energy density for 1 trial (60 sec)
    L488_power_intensity * 60
    L561_power_intensity * 60
    
    
        
%% -- Plot --
% AO imaging system with YFP
    f = figure('position', [541         152        1386         719]);
        area(Flurophore.YFP.WAVELENGTH, Flurophore.YFP.NORM_EMISSION, 'FaceColor', colorHex2RGB('#faf966'), 'EdgeColor', 'none'); hold on
        plot(Flurophore.YFP.WAVELENGTH, Flurophore.YFP.NORM_EXCITATION, 'color', colorHex2RGB('#c8c751')); hold on
        
        plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION), hold on
        
        plot([488 488], [0 1], 'color', colorHex2RGB('#1502ea'), 'LineStyle', '--', 'LineWidth', 1.5)
        plot([515 515], [0 1], 'LineStyle', '--', 'LineWidth', 1.5)
        plot([561 561], [0 1], 'LineStyle', '--', 'LineWidth', 1.5)
        plot([640 640], [0 1], 'LineStyle', '--', 'LineWidth', 1.5)

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('Imaging YFP')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('YFP emssion', 'YFP excitation', 'FF01-520/35 - emission filter', ...
            '488 nm', '515 nm', '561 nm', '640 nm',  'location', 'northeast');
        set(l, 'Fontsize', 12, 'Box', 'off')          
        
%% -- plot --
% AO imaging system with GFP, GCaMP
    f = figure('position', [541         152        1386         719]);
        area(Flurophore.GFP.WAVELENGTH, Flurophore.GFP.NORM_EMISSION, 'FaceColor', colorHex2RGB('#faf966'), 'EdgeColor', 'none'); hold on
        plot(Flurophore.GFP.WAVELENGTH, Flurophore.GFP.NORM_EXCITATION, 'color', colorHex2RGB('#c8c751')); hold on
        
        area(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, 'FaceColor', colorHex2RGB('#00aba9'), 'EdgeColor', 'none'); hold on
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'color', colorHex2RGB('#66cccb')); hold on
        
        plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION), hold on
        
        plot([488 488], [0 1], 'color', colorHex2RGB('#1502ea'), 'LineStyle', '--', 'LineWidth', 1.5)
        plot([515 515], [0 1], 'LineStyle', '--', 'LineWidth', 1.5)
        plot([561 561], [0 1], 'LineStyle', '--', 'LineWidth', 1.5)
        plot([640 640], [0 1], 'LineStyle', '--', 'LineWidth', 1.5)

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('Imaging YFP')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('GFP emssion', 'GFP excitation', 'GCaMP5g emssion', 'GCaMP5g excitation', 'FF01-520/35 - emission filter', ...
            '488 nm', '515 nm', '561 nm', '640 nm',  'location', 'northwest');
        set(l, 'Fontsize', 12, 'Box', 'off')  
        
%% -- plot -- 
    % compute photon density delivered by light stimulus
    
    % parameters:
    led_fov     = 8.6; % degrees circular
    wavelength  = 565; % nm
    light_power = logspace(-1, 2.3802, 100); % uW
    
    % compute fov on retina in cm^2
    r_cm = led_fov / 2 * 0.0032; % 32 um per degree (mouse) (cm)
    r_mm = r_cm * 10; % (mm)
    a_cm = pi * r_cm^2; % area on retina
    
    % compute quanta
    for qq = 1:length(light_power)
        quanta(qq) = power2density(light_power(qq),  wavelength);   
    end
    
    % compute density
    density = quanta ./ a_cm;
    
    % plot
    figure
        semilogy(light_power, density, 'o')
    
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title(sprintf('Power to density (%g nm, fov = %g deg)', wavelength, led_fov))
        xlabel('power (uW)')
        ylabel('photon density (photons.sec^{-1}.cm^{-2})')
        
   
    % calcuate irradiance (mW / mm^2)
        led_power = 100; % µW
        a_mm = pi * r_mm^2; % area on retina (mm)
        irrad = (led_power / 1000) / a_mm % (mW / mm^2)

        % back calculate irradiance to power
        irrad = 0.015; % mW/mm^2
        led_power = irrad * a_mm * 1000
        
        
    
%% -- Plot -- 
    f = figure('position', [541         152        1386         719]);
%         area(LED.M565L3.WAVELENGTH, LED.M565L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#f9d59c'), 'EdgeColor', 'none'); hold on
%         area(LED.M455L3.WAVELENGTH, LED.M455L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#1544ee'), 'EdgeColor', 'none'); hold on
        area(System.M565L3_system.WAVELENGTH, System.M565L3_system.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#f9d59c'), 'EdgeColor', 'none'); hold on
        area(System.M455L3_system.WAVELENGTH, System.M455L3_system.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#1544ee'), 'EdgeColor', 'none'); hold on
        area(LED.M365L2.WAVELENGTH, LED.M365L2.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#a200ff'), 'EdgeColor', 'none'); hold on

        plot([488 488], [0 1], '-.', 'LineWidth', 2.5);
        
        plot(wavelengths, mouse_mcone, 'LineWidth', 1.5); hold on
        plot(wavelengths, mouse_scone, 'LineWidth', 1.5); hold on
        plot(wavelengths, mouse_rods , 'LineWidth', 1.5); hold on
        plot(wavelengths, mouse_melanopsin, 'LineWidth', 1.5); hold on
        
        plot(wavelengths, chr2, '--'), hold on
        
        plot(Biology.ReaChR_max.WAVELENGTH, Biology.ReaChR_max.NORM_RESPONSE, '--')
        plot(Biology.ReaChR_steady.WAVELENGTH, Biology.ReaChR_steady.NORM_RESPONSE, '--')
        
        plot(Biology.C1V1.WAVELENGTH, Biology.C1V1.NORM_RESPONSE, '--')
        
%         plot(Flurophore.RCaMP2.WAVELENGTH, Flurophore.RCaMP2.NORM_EXCITATION)
%         plot(Flurophore.RCaMP2.WAVELENGTH, Flurophore.RCaMP2.NORM_EMISSION)

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Opsin action spectra and stimulation LEDs')
        xlim([350 700])
        ylim([0 1.01])

        l = legend('565 - system', '455 - system', '365 - datasheet', '488 laser', 'mcone', ...
                   'scone', 'rod', 'melanopsin', 'chr2', 'ReaChR max', 'ReaChR steady', 'C1V1');

               
        set(l, 'Fontsize', 12, 'Box', 'off')

%% -- plot --
% Imaging red fluorescence in model eye
    f = figure('position', [541         152        1386         719]);
        plot([561 561], [0 1], '-.', 'LineWidth', 2.5); hold on
        
        plot(Flurophore.RhodamineRed.WAVELENGTH, Flurophore.RhodamineRed.NORM_EXCITATION, 'LineWidth', 1.5); hold on
        plot(Flurophore.RhodamineRed.WAVELENGTH, Flurophore.RhodamineRed.NORM_EMISSION, 'LineWidth', 1.5); hold on
        plot(Flurophore.Alexa568.WAVELENGTH, Flurophore.Alexa568.NORM_EXCITATION, 'LineWidth', 1.5); hold on
        plot(Flurophore.Alexa568.WAVELENGTH, Flurophore.Alexa568.NORM_EMISSION, 'LineWidth', 1.5); hold on
        
        plot(Filter.FF01_630_92.WAVELENGTH, Filter.FF01_630_92.TRANSMISSION, '--'); hold on
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Alexa 568 vs Texas Red')
        xlim([350 700])
        ylim([0 1.01])

        l = legend('561 laser', 'RhodamineRed - excitation', 'RhodamineRed - emission', 'Alexa568 - excitation', 'Alexa568 - emission', 'FF01-630/92');
               
        set(l, 'Fontsize', 12, 'Box', 'off')   


%% -- plot --
% GCaMP marmoset experiment
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % LED stimuli
        area(LED.M420L3.WAVELENGTH, LED.M420L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#a566ff'), 'EdgeColor', 'none'); hold on
        area(LED.M590L3.WAVELENGTH, LED.M590L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#f9d59c'), 'EdgeColor', 'none'); hold on
        
        % photoreceptors
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.marmoset_scone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.marmoset_543,   'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.marmoset_556,   'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.marmoset_563,   'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.marmoset_rod,   'LineWidth', 1.5, 'Color', 'k'); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.marmoset_melanopsin, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '--'); hold on
        
        % calcium indicator
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on

        % barrier filter
        plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION, 'color', colorHex2RGB('#00aba9'), 'LineStyle', ':', 'LineWidth', 1.5); hold on
        
        % imaging laser
        plot(System.L488_system.WAVELENGTH, System.L488_system.NORM_INTENSITY, 'LineWidth', 2)
        
        % channelrhodopsins
        plot(Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
%         plot(Biology.chronos.WAVELENGTH, Biology.chronos.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Marmoset photoreceptor sensitivity')
        xlim([350 700])
        ylim([0 1.01])

        l = legend('LED - 420L3', 'LED - 590L3', 'marmoset s-cone', 'marmoset cone 543', 'marmoset cone 556', 'marmoset cone 563', 'marmoset rod', 'marmoset melanopsin',...
            'GCaMP5g - excitation', 'GCaMP5g - emission', '520/35 barrier filter', '488 nm imaging laser', 'Chrimson');
               
        set(l, 'Fontsize', 12, 'Box', 'off')     
        
        save_name = 'GCaMP marmoset vision restoration experiment.pdf';
        save_dir  = fullfile('.', 'Figures');
        save_path = fullfile(save_dir, save_name);
        exportpdf(f, save_path)        
        
%% -- plot --
% macaque GCaMP experiment
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % LED stimuli
        area(LED.M420L3.WAVELENGTH, LED.M420L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#a566ff'), 'EdgeColor', 'none'); hold on
        area(LED.M590L3.WAVELENGTH, LED.M590L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#f9d59c'), 'EdgeColor', 'none'); hold on
        
        % photoreceptors
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_scone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_mcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_lcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_rod, 'LineWidth', 1.5, 'Color', 'k'); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_melanopsin, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '--'); hold on
        
        % calcium indicator
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on

        % barrier filter
        plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION, 'color', colorHex2RGB('#00aba9'), 'LineStyle', ':', 'LineWidth', 1.5); hold on
        
        % imaging laser
        plot(System.L488_system.WAVELENGTH, System.L488_system.NORM_INTENSITY, 'LineWidth', 2)
        
        % channelrhodopsins
%         plot(Biology.ReaChR_max.WAVELENGTH, Biology.ReaChR_max.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
%         plot(Biology.CatCh.WAVELENGTH, Biology.CatCh.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
        plot(Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
%         plot(Biology.C1V1.WAVELENGTH, Biology.C1V1.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
%         plot(Biology.chronos.WAVELENGTH, Biology.chronos.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Macaque photoreceptor sensitivity')
        xlim([350 700])
        ylim([0 1.01])

        l = legend('LED - 420L3', 'LED - 590L3', 'macaque s-cone', 'macaque m-cone', 'macaque l-cone', 'macaque rod', 'macaque melanopsin',...
            'GCaMP5g - excitation', 'GCaMP5g - emission', '520/35 barrier filter', '488 nm imaging laser', 'Chrimson');
               
        set(l, 'Fontsize', 12, 'Box', 'off')     
        
        save_name = 'GCaMP macaque vision restoration experiment.pdf';
        save_dir  = fullfile('.', 'Figures');
        save_path = fullfile(save_dir, save_name);
        exportpdf(f, save_path)

%% -- plot --
% macaque 1P JRGECO1a experiment
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % stimuli
        area(LED.M455L3.WAVELENGTH, LED.M455L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#728ef4'), 'EdgeColor', 'none'); hold on
%         area(LED.M470L3.WAVELENGTH, LED.M470L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#728ef4'), 'EdgeColor', 'none'); hold on
        area(LED.M490L3.WAVELENGTH, LED.M490L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#81def2'), 'EdgeColor', 'none'); hold on
        
        % photoreceptors
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_scone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_mcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_lcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_rod, 'LineWidth', 1.5, 'Color', 'k'); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_melanopsin, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '--'); hold on
        
        % R-GECO1
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on

        % barrier filter
        plot(Filter.FF01_630_92.WAVELENGTH, Filter.FF01_630_92.TRANSMISSION), hold on
         
        % imaging laser (1070)
        plot(System.L561_system.WAVELENGTH, System.L561_system.NORM_INTENSITY, 'LineWidth', 2)
        
        % channelrhodopsins
        plot(Biology.CatCh.WAVELENGTH, Biology.CatCh.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
        plot(Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
        plot(Biology.chronos.WAVELENGTH, Biology.chronos.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on       

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Macaque jRGECO1a vision restoration experiment')
        xlim([350 700])
        ylim([0 1.01])

        l = legend('LED - 455L3', 'LED - 490L3', 'macaque s-cone', 'macaque m-cone', 'macaque l-cone', 'macaque rod', 'macaque melanopsin',...
            'R-GECO1 - excitation', 'R-GECO1 - emission', 'FF01-630/92 - emission filter', '561 nm imaging laser',...
            'CatCh', 'ChrimsonR', 'chronos');
               
        set(l, 'Fontsize', 12, 'Box', 'off')   
        
    save_name = 'JRGECO1a macaque vision restoration experiment.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)        
        
%% -- plot --
% macaque 2P JRGECO1a experiment
    f = figure('position', [541         152        1386         719]);
%         % stimuli
%         area(LED.M420L3.WAVELENGTH, LED.M420L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#a566ff'), 'EdgeColor', 'none'); hold on
        
        % photoreceptors
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_scone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_mcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_lcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_rod, 'LineWidth', 1.5, 'Color', 'k'); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_melanopsin, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '--'); hold on
        
        % R-GECO1
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EXCITATION_2P, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on

        % barrier filter
         
        % imaging laser (1070)
        plot([1070 1070], [0 1], 'LineWidth', 2)
        
%         % debug
%         plot(Flurophore.RGECO1.WAVELENGTH .* 2, Flurophore.RGECO1.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on        

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Macaque photoreceptor sensitivity')
        xlim([350 1100])
        ylim([0 1.01])

        l = legend('macaque s-cone', 'macaque m-cone', 'macaque l-cone', 'macaque rod', 'macaque melanopsin',...
            'R-GECO1 - excitation', 'R-GECO1 - 2P excitation', 'R-GECO1 - emission', '1070 nm imaging laser');
               
        set(l, 'Fontsize', 12, 'Box', 'off')   
        
    save_name = 'JRGECO1a 2P macaque experiment.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)
        

%% -- plot --
% LED power measurements with filters - 455
    f = figure('position', [541         152        1386         719]);
    
        % 455 - datasheet
        plot(LED.M455L3.WAVELENGTH, LED.M455L3.NORM_INTENSITY, 'LineStyle', ':', 'LineWidth', 1.5); hold on
        
        % 455 LED with 475/42 filter - measured
        plot(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        
        % 455 LED with 475/42 + 455/10 filter - measured
        plot(System.M455L3_475_42_455_10.WAVELENGTH, System.M455L3_475_42_455_10.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        
        % compute percentage of power measurement spectra from stimulus
        power_spectra_pct = trapz(System.M455L3_475_42_455_10.WAVELENGTH, System.M455L3_475_42_455_10.NORM_INTENSITY) / trapz(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY) * 100;
        text(0.01, 0.95, sprintf('Spectra area ratio %% = %0.2f', power_spectra_pct), 'units', 'normalized', 'fontsize', 12)
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Power Measurement - 455 LED')
        xlim([350 700])
        ylim([0 1.01])

        l = legend('LED 455 - datasheet', 'LED 455 - 475/42 - measured', 'LED 455 - 475/42 + 455/10 - measured');
               
        set(l, 'Fontsize', 12, 'Box', 'off') 
        
% LED power measurements with filters - 565
    f = figure('position', [541         152        1386         719]);
    
        plot(LED.M565L3.WAVELENGTH, LED.M565L3.NORM_INTENSITY, 'LineStyle', ':', 'LineWidth', 1.5); hold on
        
        % 565 LED with 571/72 filter
        plot(System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        
        % 565 LED with 571/72 + 568/10 filter
        plot(System.M565L3_571_72_568_10.WAVELENGTH, System.M565L3_571_72_568_10.NORM_INTENSITY, 'LineWidth', 1.5); hold on       
        
        % compute percentage of power measurement spectra from stimulus
        power_spectra_pct = trapz(System.M565L3_571_72_568_10.WAVELENGTH, System.M565L3_571_72_568_10.NORM_INTENSITY) / trapz(System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY) * 100;
        text(0.01, 0.95, sprintf('Spectra area ratio %% = %0.2f', power_spectra_pct), 'units', 'normalized', 'fontsize', 12)
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Power Measurement - 565 LED')
        xlim([350 700])
        ylim([0 1.01])

        l = legend('LED 565 - datasheet', 'LED 565 - 571/72 - measured', 'LED 565 - 571/72 + 568/10 - measured');
               
        set(l, 'Fontsize', 12, 'Box', 'off') 
                
%% -- plot --
% Measured system light sources
    f = figure('position', [541         152        1386         719]);
        % 365 LED no filter
        plot(System.M365L2_system.WAVELENGTH, System.M365L2_system.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        
        % 455 LED with 475/42 filter
        plot(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        
        % 565 LED with 571/72 filter
        plot(System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        
        % Toptica laser sources
        plot(System.L488_system.WAVELENGTH, System.L488_system.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        plot(System.L515_system.WAVELENGTH, System.L515_system.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        plot(System.L561_system.WAVELENGTH, System.L561_system.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        plot(System.L640_system.WAVELENGTH, System.L640_system.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        
        % SLD High Power - IR imaging laser
        plot(System.SLD_HP_system.WAVELENGTH, System.SLD_HP_system.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        
        % Wavefront diode - IR beacon
        plot(System.Wavefront_system.WAVELENGTH, System.Wavefront_system.NORM_INTENSITY, 'LineWidth', 1.5); hold on
        

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Mouse 1P Light Sources')
        xlim([350 900])
        ylim([0 1.01])

        l = legend('LED 365 - system', 'LED 455 - 475/42 - system', 'LED 565 - 571/72 - system', 'Toptica - 488', 'Toptica - 515', 'Toptica - 561', 'Toptica - 640', 'SLD HP - system', 'Wavefront - system');
               
        set(l, 'Fontsize', 12, 'Box', 'off') 
        
%% -- plot --
% GCaMP vs YFP - filter selection, signal isolation
f = figure('position', [541         152        1386         719]);
subplot(2,1,1)
    % excitation
    plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'LineWidth', 1.5); hold on
    plot(Flurophore.YFP.WAVELENGTH, Flurophore.YFP.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '--'); hold on

    % 488
    plot(System.L488_system.WAVELENGTH, System.L488_system.NORM_INTENSITY, 'LineWidth', 2)

    % 515
    plot(System.L515_system.WAVELENGTH, System.L515_system.NORM_INTENSITY, 'LineWidth', 2)

    % GCaMP filter
    plot(Filter.FF01_504_12.WAVELENGTH, Filter.FF01_504_12.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on
    plot(Filter.FF01_500_10.WAVELENGTH, Filter.FF01_500_10.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on

    box off
    set(gca, 'tickdir', 'out', 'FontSize', 12)
    xlabel('wavelength (nm)')
    ylabel('normalised sensitivity/intensity')
    title('excitation')
    xlim([450 600])
    ylim([0 1.01])

    l = legend('GCaMP5g excitation', 'YFP excitation', '488', '515', '504/12 filter', '500/10 filter');

    set(l, 'Fontsize', 12, 'Box', 'off')


subplot(2,1,2)
    % emission
    plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, 'LineWidth', 1.5); hold on
    plot(Flurophore.YFP.WAVELENGTH, Flurophore.YFP.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '--'); hold on

    % current green channel filter
    plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on

    % GCaMP filter
    plot(Filter.FF01_500_10.WAVELENGTH, Filter.FF01_500_10.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on
    plot(Filter.FF01_504_12.WAVELENGTH, Filter.FF01_504_12.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on

    plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on

    % YFP filter
    plot(Filter.FF01_550_49.WAVELENGTH, Filter.FF01_550_49.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on

        % calculate area under cuves for filters and emision spectra
        Flurophore.GCaMP5g.NORM_EMISSION(isnan(Flurophore.GCaMP5g.NORM_EMISSION)) = 0;
        area_GCaMP5g = trapz(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION); % a.u.
        area_YFP = trapz(Flurophore.YFP.WAVELENGTH, Flurophore.YFP.NORM_EMISSION); % a.u.

        GCaMP5g_500 = crossProduct(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, Filter.FF01_500_10.WAVELENGTH, Filter.FF01_500_10.TRANSMISSION);
        GCaMP5g_504 = crossProduct(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, Filter.FF01_504_12.WAVELENGTH, Filter.FF01_504_12.TRANSMISSION);
        GCaMP5g_550 = crossProduct(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, Filter.FF01_550_49.WAVELENGTH, Filter.FF01_550_49.TRANSMISSION);
        GCaMP5g_520 = crossProduct(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION);

        YFP_500 = crossProduct(Flurophore.YFP.WAVELENGTH, Flurophore.YFP.NORM_EMISSION, Filter.FF01_500_10.WAVELENGTH, Filter.FF01_500_10.TRANSMISSION);
        YFP_504 = crossProduct(Flurophore.YFP.WAVELENGTH, Flurophore.YFP.NORM_EMISSION, Filter.FF01_504_12.WAVELENGTH, Filter.FF01_504_12.TRANSMISSION);
        YFP_550 = crossProduct(Flurophore.YFP.WAVELENGTH, Flurophore.YFP.NORM_EMISSION, Filter.FF01_550_49.WAVELENGTH, Filter.FF01_550_49.TRANSMISSION);

        ratio_GCaMP5g_500_GCaMP5g = GCaMP5g_500.INTEGRAL/area_GCaMP5g;
        ratio_GCaMP5g_504_GCaMP5g = GCaMP5g_504.INTEGRAL/area_GCaMP5g;
        ratio_GCaMP5g_550_GCaMP5g = GCaMP5g_550.INTEGRAL/area_GCaMP5g;
        ratio_GCaMP5g_520_GCaMP5g = GCaMP5g_520.INTEGRAL/area_GCaMP5g;                

        ratio_YFP_500_YFP = YFP_500.INTEGRAL/area_YFP;
        ratio_YFP_504_YFP = YFP_504.INTEGRAL/area_YFP;
        ratio_YFP_550_YFP = YFP_550.INTEGRAL/area_YFP;

    text(0.01, 0.95, sprintf('500/10: GCaMP5g = %0.2f %%; YFP = %0.2f %%', ratio_GCaMP5g_500_GCaMP5g * 100, ratio_YFP_500_YFP * 100), 'units', 'normalized')
    text(0.01, 0.90, sprintf('504/12: GCaMP5g = %0.2f %%; YFP = %0.2f %%', ratio_GCaMP5g_504_GCaMP5g * 100, ratio_YFP_504_YFP * 100), 'units', 'normalized')
    text(0.01, 0.85, sprintf('550/49: GCaMP5g = %0.2f %%; YFP = %0.2f %%', ratio_GCaMP5g_550_GCaMP5g * 100, ratio_YFP_550_YFP * 100), 'units', 'normalized')
    text(0.01, 0.80, sprintf('520/35: GCaMP5g = %0.2f %%', ratio_GCaMP5g_520_GCaMP5g * 100), 'units', 'normalized')

    box off
    set(gca, 'tickdir', 'out', 'FontSize', 12)
    xlabel('wavelength (nm)')
    ylabel('normalised sensitivity/intensity')
    title('emission')
    xlim([450 600])
    ylim([0 1.01])

    l = legend('GCaMP5g emission', 'YFP emission', '520/35 filter', '500/10 filter', '504/12 filter', '520/70 filter', '550/49 filter');

    set(l, 'Fontsize', 12, 'Box', 'off')        

%% -- plot --
% GCaMP - filter selection
f = figure('position', [541         152        1386         719]);
    % excitation
    plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '--'); hold on

    % 488
    plot(System.L488_system.WAVELENGTH, System.L488_system.NORM_INTENSITY, 'LineWidth', 2)

    % emission
    plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, 'LineWidth', 1.5); hold on

    % GCaMP filter
%     plot(Filter.FF01_500_10.WAVELENGTH, Filter.FF01_500_10.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on
%     plot(Filter.FF01_504_12.WAVELENGTH, Filter.FF01_504_12.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on
    plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on
    plot(Filter.FF01_525_45.WAVELENGTH, Filter.FF01_525_45.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on
    plot(Filter.FF01_530_55.WAVELENGTH, Filter.FF01_530_55.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on
    plot(Filter.ET535_70.WAVELENGTH, Filter.ET535_70.TRANSMISSION, 'LineWidth', 1.5, 'LineStyle', ':'); hold on


        % calculate area under cuves for filters and emision spectra
        Flurophore.GCaMP5g.NORM_EMISSION(isnan(Flurophore.GCaMP5g.NORM_EMISSION)) = 0;
        area_GCaMP5g = trapz(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION); % a.u.

        GCaMP5g_500 = crossProduct(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, Filter.FF01_500_10.WAVELENGTH, Filter.FF01_500_10.TRANSMISSION);
        GCaMP5g_504 = crossProduct(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, Filter.FF01_504_12.WAVELENGTH, Filter.FF01_504_12.TRANSMISSION);
        GCaMP5g_520 = crossProduct(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION);
        GCaMP5g_525 = crossProduct(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, Filter.FF01_525_45.WAVELENGTH, Filter.FF01_525_45.TRANSMISSION);
        GCaMP5g_530 = crossProduct(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, Filter.FF01_530_55.WAVELENGTH, Filter.FF01_530_55.TRANSMISSION);
        GCaMP5g_535 = crossProduct(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, Filter.ET535_70.WAVELENGTH, Filter.ET535_70.TRANSMISSION);

      
        ratio_GCaMP5g_500_GCaMP5g = GCaMP5g_500.INTEGRAL/area_GCaMP5g;
        ratio_GCaMP5g_504_GCaMP5g = GCaMP5g_504.INTEGRAL/area_GCaMP5g;
        ratio_GCaMP5g_520_GCaMP5g = GCaMP5g_520.INTEGRAL/area_GCaMP5g;
        ratio_GCaMP5g_525_GCaMP5g = GCaMP5g_525.INTEGRAL/area_GCaMP5g; 
        ratio_GCaMP5g_530_GCaMP5g = GCaMP5g_530.INTEGRAL/area_GCaMP5g; 
        ratio_GCaMP5g_535_GCaMP5g = GCaMP5g_535.INTEGRAL/area_GCaMP5g; 

    
    text(0.01, 0.95, sprintf('500/10: GCaMP5g = %0.2f %%; vs 520/35 = %0.2f %%', ratio_GCaMP5g_500_GCaMP5g * 100, ratio_GCaMP5g_500_GCaMP5g / ratio_GCaMP5g_520_GCaMP5g * 100), 'units', 'normalized')
    text(0.01, 0.90, sprintf('504/12: GCaMP5g = %0.2f %%; vs 520/35 = %0.2f %%', ratio_GCaMP5g_504_GCaMP5g * 100, ratio_GCaMP5g_504_GCaMP5g / ratio_GCaMP5g_520_GCaMP5g * 100), 'units', 'normalized')
    text(0.01, 0.85, sprintf('520/35: GCaMP5g = %0.2f %%; vs 520/35 = %0.2f %%', ratio_GCaMP5g_520_GCaMP5g * 100, ratio_GCaMP5g_520_GCaMP5g / ratio_GCaMP5g_520_GCaMP5g * 100), 'units', 'normalized', 'color', 'r')
    text(0.01, 0.80, sprintf('525/45: GCaMP5g = %0.2f %%; vs 520/35 = %0.2f %%', ratio_GCaMP5g_525_GCaMP5g * 100, ratio_GCaMP5g_525_GCaMP5g / ratio_GCaMP5g_520_GCaMP5g * 100), 'units', 'normalized')
    text(0.01, 0.75, sprintf('530/55: GCaMP5g = %0.2f %%; vs 520/35 = %0.2f %%', ratio_GCaMP5g_530_GCaMP5g * 100, ratio_GCaMP5g_530_GCaMP5g / ratio_GCaMP5g_520_GCaMP5g * 100), 'units', 'normalized')
    text(0.01, 0.70, sprintf('535/70: GCaMP5g = %0.2f %%; vs 520/35 = %0.2f %%', ratio_GCaMP5g_535_GCaMP5g * 100, ratio_GCaMP5g_535_GCaMP5g / ratio_GCaMP5g_520_GCaMP5g * 100), 'units', 'normalized')
    

    box off
    set(gca, 'tickdir', 'out', 'FontSize', 12)
    xlabel('wavelength (nm)')
    ylabel('normalised sensitivity/intensity')
    title('emission')
    xlim([450 600])
    ylim([0 1.01])

    l = legend('GCaMP5g excitation', '488', 'GCaMP5g emission', '520/35 filter', '525/45 filter', '530/55 filter', '535/70 filter');

    set(l, 'Fontsize', 12, 'Box', 'off')        


%% -- calculate actual 365 power --
    % 365 power was measured at 400, the limit of the power meter
    % convert to 365 nm power measurements
    
    density_365 = power2density(PowerCurves.Power_400_20160103.POWER, 400);
    power_365   = density2power(density_365, 365);

    
%% -- Plot -- 
% micron imaging spectra
% evaluation of more optimal filters
% micron filters:
% excitation: 
%   - ch2 (green) - 498 SP (FF01_498_SP)
%       - ch2 (green) - 469/35 (FF01_469_35) % old
%   - ch3 (red)   - 525/40 (FF01_525_45)
% barrier:
%   - ch2 (green) - BLP01_488R.txt % AO 520/35
%   - ch3 (red)   - 650/150? (FF01-650/150-25)

f = figure('position', [917          20        1557         933], 'color', [1 1 1]);
    a = subplot(2,1,1); % CH2
        % Excitation filters
        plot(Filter.FF01_498_SP.WAVELENGTH, Filter.FF01_498_SP.TRANSMISSION, '--'); hold on % CONFIRMED
%         plot(Filter.FF01_469_35.WAVELENGTH, Filter.FF01_469_35.TRANSMISSION, '--'); hold on
%         plot(Filter.FF02_475_50.WAVELENGTH, Filter.FF02_475_50.TRANSMISSION, '--'); hold on
        
        % Barrier filters
        plot(Filter.BLP01_488R.WAVELENGTH, Filter.BLP01_488R.TRANSMISSION, '--'); hold on % CONFIRMED
        plot(System.micron_CH2_emission.WAVELENGTH, System.micron_CH2_emission.NORM_INTENSITY, '--'); hold on 
%         plot(Filter.FF03_525_50.WAVELENGTH, Filter.FF03_525_50.TRANSMISSION, '--'); hold on
%         plot(Filter.FF01_525_45.WAVELENGTH, Filter.FF01_525_45.TRANSMISSION, '--'); hold on
%         plot(Filter.FF02_525_40.WAVELENGTH, Filter.FF02_525_40.TRANSMISSION, '--'); hold on
        plot(Filter.ET525_50.WAVELENGTH, Filter.ET525_50.TRANSMISSION, '--'); hold on
%         plot(Filter.FF01_535_50.WAVELENGTH, Filter.FF01_535_50.TRANSMISSION, '--'); hold on
%         plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION, 'color', colorHex2RGB('#36ff00'), 'LineStyle', '-.'), hold on
        
        % Flurophores
%         plot(Flurophore.GFP.WAVELENGTH, Flurophore.GFP.NORM_EXCITATION, '-'); hold on
%         plot(Flurophore.GFP.WAVELENGTH, Flurophore.GFP.NORM_EMISSION, '-'); hold on
        plot(Flurophore.YFP.WAVELENGTH, Flurophore.YFP.NORM_EXCITATION, 'y-'); hold on
        plot(Flurophore.YFP.WAVELENGTH, Flurophore.YFP.NORM_EMISSION, 'y-.'); hold on
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EXCITATION, 'g-'); hold on
        plot(Flurophore.GCaMP5g.WAVELENGTH, Flurophore.GCaMP5g.NORM_EMISSION, 'g-.');  hold on
        plot(Flurophore.peredox.WAVELENGTH, Flurophore.peredox.NORM_EXCITATION, 'c-'); hold on
        plot(Flurophore.peredox.WAVELENGTH, Flurophore.peredox.NORM_EMISSION, 'c-.'); hold on
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('Micron imaging')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
%         '525/50 proposed CH2 filter', '525/45 proposed CH2 filter', '525/40 proposed CH2 filter'
        
        l = legend('FF01-498 SP CH2 Excitation', 'BLP01-488R CH2 Barrier', 'measured CH2 Barrier', 'ET525/50 proposed CH2 filter', ...
... %             'GFP - excitation', 'GFP - emission',...
            'YFP - excitation', 'YFP - emission',...
            'GCaMP5g - excitation', 'GCaMP5g - emission',... 
            'peredox - excitation', 'peredox - emission',... 
            'location', 'northeast');
        set(l, 'Fontsize', 12, 'Box', 'off') 
        set(a, 'color', [.8 .8 .8])
        
    a = subplot(2,1,2); % CH3
        % Excitation filters
        plot(Filter.FF01_525_45.WAVELENGTH, Filter.FF01_525_45.TRANSMISSION, '--'); hold on % CONFIRMED
        plot(Filter.FF01_538_40.WAVELENGTH, Filter.FF01_538_40.TRANSMISSION, '--'); hold on % proposed
%         plot(Filter.ET535_50.WAVELENGTH, Filter.ET535_50.TRANSMISSION, '--'); hold on % proposed
        
        
        % Barrier filters
        plot(System.micron_CH3_emission.WAVELENGTH, System.micron_CH3_emission.NORM_INTENSITY, '--'); hold on % CONFIRMED
%         plot(Filter.FF01_650_150.WAVELENGTH, Filter.FF01_650_150.TRANSMISSION, '--'); hold on DOES NOT MATCH THE INSTALLED FILTER
        
        % Flurophores
        plot(Flurophore.tdTomato.WAVELENGTH, Flurophore.tdTomato.NORM_EXCITATION, 'y-'); hold on
        plot(Flurophore.tdTomato.WAVELENGTH, Flurophore.tdTomato.NORM_EMISSION, 'y-.'); hold on
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EXCITATION, 'm-'); hold on
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EMISSION, 'm-.'); hold on
        plot(Flurophore.DsRed.WAVELENGTH, Flurophore.DsRed.NORM_EXCITATION, 'r-'); hold on
        plot(Flurophore.DsRed.WAVELENGTH, Flurophore.DsRed.NORM_EMISSION, 'r-.'); hold on        

        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('Micron imaging')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('525/45 CH3 Excitation', '538/40 proposed CH3 Excitation', 'measured CH3 Barrier',...
            'tdTomato - excitation', 'tdTomato - emission',...
            'RGECO1 - excitation', 'RGECO1 - emission',...
            'DsRed - excitation', 'DsRed - emission',...
            'location', 'northeast');
        set(l, 'Fontsize', 12, 'Box', 'off') 
        set(a, 'color', [.8 .8 .8])
    
        
    save_name = 'Micron imaging.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)    

%% -- Compute light power of micron --
    % compute FOV 
    FOV_mm  = 2.4893;        % diameter in mm
    FOV_um  = FOV_mm * 1000; % diameter in µm
    FOV_cm  = FOV_mm / 10;   % diameter in cm
    FOV_deg = FOV_um / 33;   % diameter in degrees
    
    FOV_area_mm = pi * (FOV_mm/2)^2; % area mm^2
    FOV_area_um = pi * (FOV_um/2)^2; % area µm^2
    FOV_area_cm = pi * (FOV_cm/2)^2; % area cm^2
    
    % light powers
    power_ch1 = [0.15 21]; % mW @ 625 nm
    power_ch2 = [0.092 7.5]; % mW @ 450 nm; filter: 498 SP (FF01_498_SP)
%     power_ch2 = [0.031 2.9]; % mW @ 625 nm; filter: 469/35 (FF01_469_35) % old
    power_ch3 = [0.031 3.9]; % mW @ 625 nm; filter: 525/40 (FF01_525_45)
    
    % power density (mW / unit area)
%     power_ch1 ./ FOV_area_mm
%     power_ch1 ./ FOV_area_um
    power_ch1 ./ FOV_area_cm
    
%     power_ch2 ./ FOV_area_mm
%     power_ch2 ./ FOV_area_um
    power_ch2 ./ FOV_area_cm    
    
%     power_ch3 ./ FOV_area_mm
%     power_ch3 ./ FOV_area_um
    power_ch3 ./ FOV_area_cm    
    
%% -- Plot -- 
% micron imaging of peredox FLIM sensor
f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
    % CH2
        % Excitation filters
        plot(Filter.FF01_498_SP.WAVELENGTH, Filter.FF01_498_SP.TRANSMISSION, '--'); hold on 

        % Barrier filters
        plot(Filter.BLP01_488R.WAVELENGTH, Filter.BLP01_488R.TRANSMISSION, '--'); hold on
        
        % Flurophores
%         plot(Flurophore.GFP.WAVELENGTH, Flurophore.GFP.NORM_EXCITATION, '-'); hold on
%         plot(Flurophore.GFP.WAVELENGTH, Flurophore.GFP.NORM_EMISSION, '-.'); hold on
        plot(Flurophore.peredox.WAVELENGTH, Flurophore.peredox.NORM_EXCITATION, '-'); hold on
        plot(Flurophore.peredox.WAVELENGTH, Flurophore.peredox.NORM_EMISSION, '-.'); hold on
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('Micron imaging')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        l = legend('FF01-498\_SP', 'BLP01-488R CH2 Barrier',...
...%             'GFP - excitation', 'GFP - emission',...
            'peredox - excitation', 'peredox - emission',...   
             'location', 'northeast');
        set(l, 'Fontsize', 12, 'Box', 'off') 
          
    save_name = 'Micron imaging of peredox FLIM sensor.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)   

%% -- Plot --
% evaluate GCaMP imaging filters with 565 stimulus led
    f = figure('position', [541         152        1386         719]);
        
        % 565 LED with 571/72 filter
        area(System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#f9d59c'), 'EdgeColor', 'none'); hold on
    
        % 565 led
%         plot(System.M565L3_system.WAVELENGTH, System.M565L3_system.NORM_INTENSITY); hold on
        
        % filters
        plot(Filter.FF01_571_72.WAVELENGTH, Filter.FF01_571_72.TRANSMISSION); hold on
        plot(Filter.FF01_580_60.WAVELENGTH, Filter.FF01_580_60.TRANSMISSION); hold on
%         plot(Filter.FF01_582_75.WAVELENGTH, Filter.FF01_582_75.TRANSMISSION); hold on
        
        % GCaMP emission filters
        plot(Filter.FF01_520_35.WAVELENGTH, Filter.FF01_520_35.TRANSMISSION), hold on
        
        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350 700])
        ylim([0 1.01])
        title('evaluation of filters for 565 led')
        
        l = legend('LED 565 - 571/72 - system', '575/72', '580/60', 'GCaMP barrier filter (520/35)');
        set(l, 'Fontsize', 12, 'Box', 'off')
 
%% -- plot --
% macaque 1P JRGECO1a experiment JEM no vision restoration
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % stimuli
        area(LED.M455L3.WAVELENGTH, LED.M455L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#728ef4'), 'EdgeColor', 'none'); hold on
%         area(LED.M470L3.WAVELENGTH, LED.M470L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#728ef4'), 'EdgeColor', 'none'); hold on
        area(LED.M490L3.WAVELENGTH, LED.M490L3.NORM_INTENSITY, 'FaceColor', colorHex2RGB('#81def2'), 'EdgeColor', 'none'); hold on
        
        % photoreceptors
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_scone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_mcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_lcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_rod, 'LineWidth', 1.5, 'Color', 'k'); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.macaque_melanopsin, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '--'); hold on
        
        % R-GECO1
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EXCITATION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
        plot(Flurophore.RGECO1.WAVELENGTH, Flurophore.RGECO1.NORM_EMISSION, 'LineWidth', 1.5, 'LineStyle', '-.'); hold on

        % barrier filter
        plot(Filter.FF01_630_92.WAVELENGTH, Filter.FF01_630_92.TRANSMISSION), hold on
         
        % imaging laser (1070)
        plot(System.L561_system.WAVELENGTH, System.L561_system.NORM_INTENSITY, 'LineWidth', 2)
        
        % channelrhodopsins
        plot(Biology.CatCh.WAVELENGTH, Biology.CatCh.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
        plot(Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on
        plot(Biology.chronos.WAVELENGTH, Biology.chronos.NORM_RESPONSE, 'LineStyle', '--', 'LineWidth', 1.5); hold on       

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        title('Macaque jRGECO1a vision restoration experiment')
        xlim([350 700])
        ylim([0 1.01])

        l = legend('LED - 455L3', 'LED - 490L3', 'macaque s-cone', 'macaque m-cone', 'macaque l-cone', 'macaque rod', 'macaque melanopsin',...
            'R-GECO1 - excitation', 'R-GECO1 - emission', 'FF01-630/92 - emission filter', '561 nm imaging laser',...
            'CatCh', 'ChrimsonR', 'chronos');
               
        set(l, 'Fontsize', 12, 'Box', 'off')   
        
    save_name = 'JRGECO1a macaque vision restoration experiment.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)
    %% -- Plot --       
% VEP brain recording LED stimuli
    f = figure('position', [541         152        1386         719], 'color', [1 1 1]);
        % LEDs
%         plot(System.M365L2_system.WAVELENGTH, System.M365L2_system.NORM_INTENSITY, 'Color', colorHex2RGB('#a200ff'), 'LineStyle', '-'); hold on
%         plot(System.M455L3_475_42.WAVELENGTH, System.M455L3_475_42.NORM_INTENSITY, 'Color', colorHex2RGB('#728ef4'), 'LineStyle', '-'); hold on
%         plot(System.M565L3_571_72.WAVELENGTH, System.M565L3_571_72.NORM_INTENSITY, 'Color', colorHex2RGB('#f9d59c'), 'LineStyle', '-'); hold on
%         plot(System.M595L3_BLP_594R.WAVELENGTH, System.M595L3_BLP_594R.NORM_INTENSITY, 'Color', colorHex2RGB('#ff9999'), 'LineStyle', '-'); hold on
        
        plot(LED.M365L2.WAVELENGTH, LED.M365L2.NORM_INTENSITY, '--'); hold on
        plot(LED.M455L3.WAVELENGTH, LED.M455L3.NORM_INTENSITY, '--'); hold on
        plot(LED.M565L3.WAVELENGTH, LED.M565L3.NORM_INTENSITY, '--'); hold on
        plot(LED.MCWHL5.WAVELENGTH, LED.MCWHL5.NORM_INTENSITY, '--'); hold on
        plot(LED.MWWHL4.WAVELENGTH, LED.MWWHL4.NORM_INTENSITY, '--'); hold on
        
        % photoreceptor sensitivites
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_mcone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_scone, 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_rods , 'LineWidth', 1.5); hold on
        plot(Photoreceptor.WAVELENGTH, Photoreceptor.mouse_melanopsin, 'LineWidth', 1.5); hold on
        
        % channelrhodopsins
        plot(Biology.chrimson.WAVELENGTH, Biology.chrimson.NORM_RESPONSE, '-'); hold on
        

        box off
        set(gca, 'tickdir', 'out', 'FontSize', 12)
        title('VEP mouse experiment')
        xlabel('wavelength (nm)')
        ylabel('normalised sensitivity/intensity')
        xlim([350  700])
        ylim([0   1.01])
        
        % '365 - measured', '455 - 475/42 - measured', '565 - 571/72 - measured'
        l = legend('365 led', '455 led', '565 led', 'MCWHL5 cool white led', 'MWWHL4 warm white led',... 
             'mcone', 'scone', 'rod', 'melanopsin', 'chrimson');
        
        set(l, 'Fontsize', 12, 'Box', 'off') 
        
    save_name = 'mouse VEP experiment.pdf';
    save_dir  = fullfile('.', 'Figures');
    save_path = fullfile(save_dir, save_name);
    exportpdf(f, save_path)

        
        
        
        
end

%% supporting functions
function output_struct = loadTablesFromFolder(table_dir)
    % load spectra .txt files from folder as datasets and return in structure

    fileobj = dir(table_dir);
    filelist = {fileobj(:).name}';
    filelist_txt_ind = cellfun(@(x) ~isempty(strfind(x, '.txt')), filelist, 'UniformOutput', true);
    filelist = filelist(filelist_txt_ind); % filter out only txt files

    for qq = 1:length(filelist) % for each file, load
        temp.field_name = filelist{qq}(1:end-4);
        output_struct.(temp.field_name) = import_dataset(fullfile(table_dir, filelist{qq}));
    end
end

function data_norm = normaliseCurves(data)  %#ok<*DEFNU>
    data_norm = data - min(data);
    data_norm = data_norm ./ max(data_norm);
end

function rgb = colorHex2RGB(hex_color)    
    % converts hex color codes to ratio RGB values [0 1]
    % disigned for specifying plot colours in hex format
    % Input hex_color as string eg. '#00ffab', '00ffab'
    
    % remove hash if present
    hash_ind = strfind(hex_color, '#'); 
    if ~isempty(hash_ind)
        hex_color(hash_ind) = [];
    end
    
    hex_color = upper(hex_color); % convert to upper case
    
    % red 
    red_hex = hex_color(1 : 2);
    red_uint8 = hex2dec(red_hex);
    r = red_uint8/255;
    
    % green
    green_hex = hex_color(3 : 4);
    green_uint8 = hex2dec(green_hex);
    g = green_uint8/255;
    
    % blue
    blue_hex = hex_color(5 : 6);
    blue_unit8 = hex2dec(blue_hex);
    b = blue_unit8/255;
    
    %
    rgb = [r g b];
end

function output_struct = crossProduct(spectra_a_wavelength, spectra_a_spectra, spectra_b_wavelength, spectra_b_spectra)
    % Computes cross product of 2 curves.
    % Uses linear interpolation to make the wavelength values of b quivalent to a 
    %   before calculating the cross product.
    spectra_b_spectra_interp = interp1(spectra_b_wavelength, spectra_b_spectra, spectra_a_wavelength);
    
    spectra_product = spectra_a_spectra .* spectra_b_spectra_interp;
    spectra_product(isnan(spectra_product)) = 0; % convet nans to 0
    spectra_product_intergral = trapz(spectra_a_wavelength, spectra_product); % integral under curve
    
    output_struct.WAVELENGTH = spectra_a_wavelength;
    output_struct.SPECTRA    = spectra_product;
    output_struct.INTEGRAL   = spectra_product_intergral;
end

function n_photons = power2density(power, lambda)
    % calculating number of photons (photon density) from light power
    % inputs:
    %   power  (uW)
    %   lambda (nm)
    % E = (n*h*c)/ lamda

    % assuming t = 1 sec, thefore watts is equivelent to joules
    power  = power  .* 10^-6; % convert to W
    lambda = lambda * 10^-9; % convert to m

    % constants
    h = 6.626 * 10 ^-34; % planks constant
    c = 2.998 * 10^8; % speed of light

    n_photons = (power .* lambda) ./ (h * c);
end

function power = density2power(density, lambda)
    % calculate light power from light density (number of photos)
    % inputs:
    %   density (n photons)
    %   lambda  (nm)
    % output:
    %   power (uW)
    % E = (n*h*c)/ lamda    
    
    lambda = lambda * 10^-9; % convert to m
    
    % constants
    h = 6.626 * 10 ^-34; % planks constant
    c = 2.998 * 10^8; % speed of light
    
    power = (density .* h .* c) ./ lambda;
    power = power .* 10^6;
end

function exportpdf(handle, save_name)
    % my own save figure as pdf function that exports the figure as you see it
    % corrects for stupid colour changing and page dimension changing
    
    set(handle, 'PaperPositionMode', 'auto')
    paperPos = get(handle, 'PaperPosition');
    set(handle, 'PaperSize', paperPos(3:4))
    set(handle, 'inverthardcopy', 'off')      % prevent chaging of background colours
    saveas(handle, save_name, 'pdf')
end



