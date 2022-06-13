% NAME:
% -----
% isis_diskbb_powerlaw.sl
%
% PURPOSE:
% --------
% Fit the GRO J1655-40 spectrum with the model "TBabs*(diskbb+powerlaw)" w/ nH as a fixed parameter.
%
% OUTPUTS: <-- To the directory /Users/salvesen/research/fcolabs/results/
% --------
% isis_diskbb_powerlaw.par - Parameter file for the model fit
% isis_diskbb_powerlaw.fit - Bestfit parameters and uncertainties
%
% NOTES:
% ------
% Calculate the unabsorbed disk flux --> http://space.mit.edu/CXC/isis/archive/2011/0777.html

%====================================================================================================
% DEFINE SOME VARIABLES

% Output filenames
variable dout  = "/Users/salvesen/research/fcolabs/results/";
variable fpars = dout + "isis_diskbb_powerlaw.par";
variable fbfit = dout + "isis_diskbb_powerlaw.fit";

% Data filenames
variable ddata = "/Users/salvesen/research/fcolabs/data/";
variable fdata = ddata + "source_annulus_grp10.pha";

% Energy range for fitting
variable Emin = 0.5;   % [keV]
variable Emax = 10.0;  % [keV]

% Galactic value for nH (Dickey & Lockman 1990)
variable nHGal = 0.74;  % [10^22 atoms cm^-2]

% Powerlaw index does not converge, so we must (somewhat arbitrarily) set it ourselves
variable Gamma = -2.3;

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
% FIT MODEL "TBabs*(diskbb+powerlaw)" W/ nH *FIXED* AND OUTPUT RESULTS

% Define the fit function as TBabs*(diskbb+powerlaw)
fit_fun("TBabs(1)*(diskbb(1)+powerlaw(1))");

% Freeze nH at its Galactic value
variable frz = 1;  % Freeze
set_par("TBabs(1).nH", nHGal, frz);

% powerlaw: Set the powerlaw index
set_par("powerlaw(1).PhoIndex", Gamma, frz, -3, 9);

% Renorm and fit
() = renorm_counts;
() = fit_counts;

% Collect parameter uncertainties (level=0 for 68% confidence)
variable Kbb_lo; variable Kbb_hi;
variable Tin_lo; variable Tin_hi;
variable Kpl_lo; variable Kpl_hi;
(Kbb_lo, Kbb_hi) = conf_loop("diskbb(1).norm", 0);
(Tin_lo, Tin_hi) = conf_loop("diskbb(1).Tin", 0);
(Kpl_lo, Kpl_hi) = conf_loop("powerlaw(1).norm", 0);

% Collect the best fit parameters
variable Kbb = get_par("diskbb(1).norm");
variable Tin = get_par("diskbb(1).Tin");
variable Kpl = get_par("powerlaw(1).norm");

% Collect the best-fit chi-square and degrees of freedom
variable stats  = eval_stat_counts();
variable chisqr = stats.statistic;
variable dof    = stats.num_bins - stats.num_variable_params;

% Save the TBabs*(diskbb+powerlaw) fit parameter file
save_par(fpars);

% Save the TBabs*(diskbb+powerlaw) best fit parameters and their 1-sigma uncertainties
fp = fopen(fbfit, "w");
() = fprintf(fp, "Best-fit parameters for model: TBabs*(diskbb+powerlaw)\n");
() = fprintf(fp, "Kbb (Kbb_lower, Kbb_upper) Tin (Tin_lower, Tin_upper) Kpl (Kpl_lower, Kpl_upper) chisqr dof\n");
() = fprintf(fp, "%8.6f (%8.6f, %8.6f) %8.6f (%8.6f, %8.6f) %8.6f (%8.6f, %8.6f) %12.6f %3.0f\n",
                  Kbb, Kbb_lo[0], Kbb_hi[0], Tin, Tin_lo[0], Tin_hi[0],
                  Kpl, Kpl_lo[0], Kpl_hi[0],
                  chisqr, dof);
() = fclose(fp);

% Exit ISIS
quit;
