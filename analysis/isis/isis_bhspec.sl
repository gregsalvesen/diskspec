% NAME:
% -----
% isis_bhspec.sl
%
% PURPOSE:
% --------
% Fit the GRO J1655-40 spectrum with the model "TBabs*bhspec".
% Only two free parameters: "log(L/Ledd)" and "fcol"
%
% OUTPUTS: <-- To the directory ../results
% --------
% bhspec.par - Parameter file for the model fit "TBabs*bhspec"
% bhspec.fit - Bestfit parameters and uncertainties for model "TBabs*bhspec"

%====================================================================================================
% DEFINE SOME VARIABLES

% Output filenames
variable dout = "/Users/salvesen/research/fcolabs/results/";
variable fpars_bhspec = dout + "bhspec.par";
variable fbfit_bhspec = dout + "bhspec.fit";

% Data filenames
variable ddata = "/Users/salvesen/research/fcolabs/data/";
variable fdata = ddata + "source_annulus_grp10.pha";

% Energy range for fitting
variable Emin = 0.5;   % [keV]
variable Emax = 10.0;  % [keV]

% Galactic value for nH (Dickey & Lockman 1990)
variable nHGal = 0.74;  % [10^22 atoms cm^-2]

% GROJ1655-40 bhspec parameters
% Need to calculate norm, logmass, inc <-- cos(inc)
variable norm     = 1.0;
variable logmass  = 1.0;
variable loglumin = -1.0;  % Unknown (free)
variable inc      = 0.5;
variable spin     = 0.92;

%====================================================================================================
% LOAD IN THE DATASET

% Workaround for the Swift RMF file not conforming to the OGIP standards.
Rmf_OGIP_Compliance=0;

% Load the data...and the associated BACKFILE, RESPFILE, ANCRFILE
variable data_id = load_data(fdata);

% Group the data (20 counts/bin equates to a S/N of 4.47 --> 4.47 = 20/sqrt(20))
group(data_id; min_sn=4.47, bounds=Emin, unit="keV");

% Only notice energies in the range 0.5-10 keV
notice_values(data_id, Emin, Emax; unit="keV");

%====================================================================================================
% (1) FIT MODEL "TBabs*bhspec"

% Add bhspec as an additive table model
%add_atable_model ("bhspec_spin2_0.1.fits", "bhspec");
add_atable_model ("bhspec2.fits", "bhspec");

% Define the fit function as TBabs*bhspec
fit_fun("TBabs(1)*bhspec(1)");

% Set the parameters: value, freeze, min, max (Note: 0=thaw, 1=freeze)
variable thw = 0;  % Thaw
variable frz = 1;  % Freeze
set_par("TBabs(1).nH",        nHGal,    frz);
set_par("bhspec(1).norm",     norm,     frz);
set_par("bhspec(1).logmass",  logmass,  frz);
set_par("bhspec(1).loglumin", loglumin, thw);
set_par("bhspec(1).inc",      inc,      frz);
set_par("bhspec(1).spin",     spin,     frz);

% Fit (Note: There is no free norm parameter to perform a renorm on)
() = fit_counts;

% Collect parameter uncertainties (level=0 for 68% confidence)
variable loglumin_lo;  variable loglumin_hi;
(loglumin_lo, loglumin_hi) = conf_loop("bhspec(1).loglumin", 0);

% Collect the best fit parameters
variable loglumin = get_par("bhspec(1).loglumin");

% Collect the best-fit chi-square and degrees of freedom
variable stats  = eval_stat_counts();
variable chisqr = stats.statistic;
variable dof    = stats.num_bins - stats.num_variable_params;

% Save the TBabs*bhspec fit parameter file
save_par(fpars_bhspec);

% Save the TBabs*bhspec best fit parameters and their 1-sigma uncertainties
fp = fopen(fbfit_bhspec, "w");
() = fprintf(fp, "Best-fit parameters and 1-sigma range for model: TBabs*bhspec\n");
() = fprintf(fp, "loglumin (loglumin_lower, loglumin_upper) chisqr     dof\n");
() = fprintf(fp, "%8.6f (%8.6f, %8.6f) %8.6f (%8.6f, %8.6f) %12.6f %3.0f\n",
             loglumin, loglumin_lo[0], loglumin_hi[0], chisqr, dof);
() = fclose(fp);

% Exit ISIS
quit;
