% NAME:
% -----
% isis_diskbb_compTT.sl
%
% PURPOSE:
% --------
% Fit the GRO J1655-40 spectrum with the model "TBabs*(diskbb+compTT)" w/ nH as a fixed parameter.
%
% OUTPUTS: <-- To the directory /Users/salvesen/research/fcolabs/results/
% --------
% isis_diskbb_compTT.par - Parameter file for the model fit
% isis_diskbb_compTT.fit - Bestfit parameters and uncertainties
%
% NOTES:
% ------
% Calculate the unabsorbed disk flux --> http://space.mit.edu/CXC/isis/archive/2011/0777.html

%====================================================================================================
% DEFINE SOME VARIABLES

% Output filenames
variable dout  = "/Users/salvesen/research/fcolabs/results/";
variable fpars = dout + "isis_diskbb_compTT.par";
variable fbfit = dout + "isis_diskbb_compTT.fit";

% Data filenames
variable ddata = "/Users/salvesen/research/fcolabs/data/";
variable fdata = ddata + "source_annulus_grp10.pha";

% Energy range for fitting
variable Emin = 0.5;   % [keV]
variable Emax = 10.0;  % [keV]

% Galactic value for nH (Dickey & Lockman 1990)
variable nHGal = 0.74;  % [10^22 atoms cm^-2]

% Electron temperature does not converge, so we must (somewhat arbitrarily) set it ourselves
variable kT = 20;

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
% FIT MODEL "TBabs*(diskbb+compTT)" W/ nH *FIXED* AND OUTPUT RESULTS

% Define the fit function as TBabs*(diskbb+compTT)
fit_fun("TBabs(1)*(diskbb(1)+compTT(1))");

% Freeze nH at its Galactic value
variable frz = 1;  % Freeze
set_par("TBabs(1).nH", nHGal, frz);

% compTT: Set the electron temperature
% compTT: Set the redshift to zero
set_par("compTT(1).kT", kT, frz);
set_par("compTT(1).Redshift", 0, frz);

% compTT: tie comptTT.T0 to diskbb.Tin (need to adjust min/max to avoid hard limit error)
variable thw = 0;  % Thaw
set_par("diskbb(1).Tin", 1.0, thw, 0.001, 100);
set_par("compTT(1).T0", 1.0, thw, 0.001, 100);
tie("diskbb(1).Tin", "compTT(1).T0");

% Renorm and fit
() = renorm_counts;
() = fit_counts;

% Collect parameter uncertainties (level=0 for 68% confidence)
variable Kbb_lo;  variable Kbb_hi;
variable Tin_lo;  variable Tin_hi;
variable KTT_lo;  variable KTT_hi;
variable taup_lo; variable taup_hi;
(Kbb_lo, Kbb_hi)   = conf_loop("diskbb(1).norm", 0);
(Tin_lo, Tin_hi)   = conf_loop("diskbb(1).Tin", 0);
(KTT_lo, KTT_hi)   = conf_loop("compTT(1).norm", 0);
(taup_lo, taup_hi) = conf_loop("compTT(1).taup", 0);

% Collect the best fit parameters
variable Kbb  = get_par("diskbb(1).norm");
variable Tin  = get_par("diskbb(1).Tin");
variable KTT  = get_par("compTT(1).norm");
variable taup = get_par("compTT(1).taup");

% Collect the best-fit chi-square and degrees of freedom
variable stats  = eval_stat_counts();
variable chisqr = stats.statistic;
variable dof    = stats.num_bins - stats.num_variable_params;

% Save the TBabs*(diskbb+compTT) fit parameter file
save_par(fpars);

% Save the TBabs*(diskbb+compTT) best fit parameters and their 1-sigma uncertainties
fp = fopen(fbfit, "w");
() = fprintf(fp, "Best-fit parameters for model: TBabs*(diskbb+compTT)\n");
() = fprintf(fp, "Kbb (Kbb_lower, Kbb_upper) Tin (Tin_lower, Tin_upper) KTT (KTT_lower, KTT_upper) taup (taup_lower, taup_upper) chisqr dof\n");
() = fprintf(fp, "%8.6f (%8.6f, %8.6f) %8.6f (%8.6f, %8.6f) %8.6f (%8.6f, %8.6f) %8.6f (%8.6f, %8.6f) %12.6f %3.0f\n",
                  Kbb, Kbb_lo[0], Kbb_hi[0], Tin, Tin_lo[0], Tin_hi[0],
                  KTT, KTT_lo[0], KTT_hi[0], taup, taup_lo[0], taup_hi[0],
                  chisqr, dof);
() = fclose(fp);

% Exit ISIS
quit;
