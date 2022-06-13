% THIS STEP OF ADDING POWERLAW, SIMPL, COMPTT TO KERRBB STILL NEEDS TO BE FINISHED

% NAME:
% -----
% isis_kerrbb_simpl.sl
%
% PURPOSE:
% --------
% Fit the GRO J1655-40 spectrum with the model "TBabs*simpl*(kerrbb)" w/ nH frozen.
% Only two free parameters: "Mdd" and "hd"
%
% OUTPUTS: <-- To the directory /Users/salvesen/research/fcolabs/results/
% --------
% isis_kerrbb.par - Parameter file for the model fit
% isis_kerrbb.fit - Bestfit parameters and uncertainties

%====================================================================================================
% DEFINE SOME VARIABLES

% Output filenames
variable dout  = "/Users/salvesen/research/fcolabs/results/";
variable fpars = dout + "isis_kerrbb_simpl.par";
variable fbfit = dout + "isis_kerrbb_simpl.fit";

% Data filenames
variable ddata = "/Users/salvesen/research/fcolabs/data/";
variable fdata = ddata + "source_annulus_grp10.pha";

% Energy range for fitting
variable Emin = 0.5;   % [keV]
variable Emax = 10.0;  % [keV]

% Galactic value for nH (Dickey & Lockman 1990)
variable nHGal = 0.74;  % [10^22 atoms cm^-2]

% GROJ1655-40 kerrbb parameters
variable norm  = 1.0;
variable eta   = 0.0;
variable a     = 0.92;
variable i     = 68.65;
variable Mbh   = 5.4;
variable Mdd   = 1.0;  % Unknown (free)
variable Dbh   = 3.2;
variable hd    = 1.7;  % Unknown (free)
variable rflag = 1;
variable lflag = 1;

% Powerlaw index does not converge, so we must (somewhat arbitrarily) set it ourselves
variable Gamma = 2.3;

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
% FIT MODEL "TBabs*simpl*(kerrbb)"

% Define the fit function as TBabs*simple*(kerrbb)
fit_fun("TBabs(1)*simpl(1, kerrbb(1))");

% Set the parameters: value, freeze, min, max (Note: 0=thaw, 1=freeze)
variable thw = 0;  % Thaw
variable frz = 1;  % Freeze
set_par("TBabs(1).nH",     nHGal, frz);
set_par("kerrbb(1).norm",  norm,  frz);
set_par("kerrbb(1).eta",   eta,   frz);
set_par("kerrbb(1).a",     a,     frz);
set_par("kerrbb(1).i",     i,     frz);
set_par("kerrbb(1).Mbh",   Mbh,   frz);
set_par("kerrbb(1).Mdd",   Mdd,   thw);
set_par("kerrbb(1).Dbh",   Dbh,   frz);
set_par("kerrbb(1).hd",    hd,    thw);
set_par("kerrbb(1).rflag", rflag, frz);
set_par("kerrbb(1).lflag", lflag, frz);

% simpl: Set the powerlaw index
% simpl: Allow for both up/down scattering
set_par("simpl(1).Gamma", Gamma, frz);
set_par("simpl(1).UpScOnly", 0, frz);

% Fit (Note: There is no free norm parameter to perform a renorm on)
() = fit_counts;

list_par;

% Collect parameter uncertainties (level=0 for 68% confidence)
variable Mdd_lo;   variable Mdd_hi;
variable hd_lo;    variable hd_hi;
variable fScat_lo; variable fScat_hi;
(Mdd_lo, Mdd_hi)     = conf_loop("kerrbb(1).Mdd", 0);
(hd_lo,  hd_hi)      = conf_loop("kerrbb(1).hd", 0);
(fScat_lo, fScat_hi) = conf_loop("simpl(1).FracSctr", 0);

% Collect the best fit parameters
variable Mdd   = get_par("kerrbb(1).Mdd");
variable hd    = get_par("kerrbb(1).hd");
variable fScat = get_par("simpl(1).FracSctr");

% Collect the best-fit chi-square and degrees of freedom
variable stats  = eval_stat_counts();
variable chisqr = stats.statistic;
variable dof    = stats.num_bins - stats.num_variable_params;

% Save the TBabs*simpl*(kerrbb) fit parameter file
save_par(fpars);

% Save the TBabs*simpl*(kerrbb) best fit parameters and their 1-sigma uncertainties
fp = fopen(fbfit, "w");
() = fprintf(fp, "Best-fit parameters and 1-sigma range for model: TBabs*simpl*(kerrbb)\n");
() = fprintf(fp, "Mdd (Mdd_lower, Mdd_upper) hd (hd_lower, hd_upper) fScat (fScat_lower, fScat_upper) chisqr dof\n");
() = fprintf(fp, "%8.6f (%8.6f, %8.6f) %8.6f (%8.6f, %8.6f) %8.6f (%8.6f, %8.6f) %12.6f %3.0f\n",
             Mdd, Mdd_lo[0], Mdd_hi[0], hd, hd_lo[0], hd_hi[0],
             fScat, fScat_lo[0], fScat_hi[0],
             chisqr, dof);
() = fclose(fp);

% Exit ISIS
quit;
