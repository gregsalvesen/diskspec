% NAME:
% -----
% isis_mcmc_kerrbb.sl
%
% PURPOSE:
% --------
% Use emcee to fit the GRO J1655-40 spectrum with the model "TBabs*kerrbb" w/ nH frozen.
%
% OUTPUTS: <-- To the directory /Users/salvesen/research/fcolabs/results/
% --------
% isis_mcmc_kerrbb.par        - parameter file used to initialize the walkers
% isis_chain_burn_kerrbb.fits - burn-in chain file
% isis_chain_mcmc_kerrbb.fits - emcee chain file

%====================================================================================================
% DEFINE SOME VARIABLES

% Output files
variable dout  = "/Users/salvesen/research/fcolabs/results/";
variable finit = dout + "isis_mcmc_kerrbb.par";         % Parameter file used to initialize the walkers
variable fburn = dout + "isis_chain_burn_kerrbb.fits";  % Burn-in chain file
variable fmcmc = dout + "isis_chain_mcmc_kerrbb.fits";  % emcee chain file

% Input parameter file from running the script: kerrbb.sl
variable fpars = dout + "isis_kerrbb.par";

% MCMC inputs
variable nwalk = 10;    % Number of walkers per free parameter
variable nburn = 500;   % Number of iterations, or steps, for each walker for the burn-in
variable nmcmc = 1000;  % Number of iterations, or steps, for each walker for the production run

% Data filenames
variable ddata = "/Users/salvesen/research/fcolabs/data/";
variable fdata = ddata + "source_annulus_grp10.pha";

% Energy range for fitting
variable Emin = 0.5;   % [keV]
variable Emax = 10.0;  % [keV]

% GRO J1655-40 known parametes and their 1-sigma uncertainties for MCMC priors
variable mean_a=0.92;                             % Adopt a uniform prior
variable mean_i=68.65;  variable sigma_i=0.91;    % Adopt a Gaussian prior
variable mean_Mbh=5.40; variable sigma_Mbh=0.18;  % Adopt a Gaussian prior
variable mean_Dbh=3.20; variable sigma_Dbh=0.2;   % Adopt a Gaussian prior

% GRO J1655-40 parameter ranges
variable a_min=0.9;   variable a_max=0.9999;
variable i_min=0.0;   variable i_max=85.0;
variable Mbh_min=0.0; variable Mbh_max=50.0;
variable Mdd_min=0.0; variable Mdd_max=10.0;
variable Dbh_min=0.0; variable Dbh_max=50.0;
variable hd_min=1.0;  variable hd_max=10.0;

%====================================================================================================
% LOAD IN THE DATASET

% Load ISISscripts
variable isisscripts="/Users/salvesen/soft/isis-scripts/share/isisscripts";
require(isisscripts);

% Workaround for the Swift RMF file not conforming to the OGIP standards.
Rmf_OGIP_Compliance=0;

% Load the data...and the associated BACKFILE, RESPFILE, ANCRFILE
variable data_id = load_data(fdata);

% Group the data (20 counts/bin equates to a S/N of 4.47 --> 4.47 = 20/sqrt(20))
group(data_id; min_sn=4.47, bounds=Emin, unit="keV");

% Only notice energies in the range 0.5-10 keV
notice_values(data_id, Emin, Emax; unit="keV");

%====================================================================================================
% emcee FIT TO THE MODEL "TBabs*kerrbb"

% Define the fit function as TBabs*kerrbb
fit_fun("TBabs(1)*kerrbb(1)");

% Load the best fit parameter file returned by the script: kerrbb.sl
load_par(fpars);

% Thaw parameters and restrict their range
variable thw = 0;  % Thaw
variable frz = 1;  % Freeze
set_par("kerrbb(1).a",   , thw, a_min,   a_max);
set_par("kerrbb(1).i",   , thw, i_min,   i_max);
set_par("kerrbb(1).Mbh", , thw, Mbh_min, Mbh_max);
set_par("kerrbb(1).Mdd", , thw, Mdd_min, Mdd_max);
set_par("kerrbb(1).Dbh", , thw, Dbh_min, Dbh_max);
set_par("kerrbb(1).hd",  , thw, hd_min,  hd_max);

% Save the parameter file used to initialize the walkers
save_par(finit);

% Define the priors
% - Uniform: a, Mdd, hd (default)
% - Gaussian: i, Mbh, Dbh (need to be specified)
variable parlist   = ["kerrbb(1).i", "kerrbb(1).Mbh", "kerrbb(1).Dbh"];
variable prior_i   = struct{mean=mean_i, sigma=sigma_i};
variable prior_Mbh = struct{mean=mean_Mbh, sigma=sigma_Mbh};
variable prior_Dbh = struct{mean=mean_Dbh, sigma=sigma_Dbh};
variable priorlist = [prior_i, prior_Mbh, prior_Dbh];
variable priors    = struct{parlist=parlist, priorlist=priorlist};

% Burn-in <-- ??? Specify init_chain=fwalk0 ??? NEED TO FIGURE OUT HOW TO DO THIS !!!
emcee(nwalk, nburn; pfile=finit, scale=0.0001, gaussian, priors=priors, output=fburn, sim_count=50, serial);

% Production run <-- Need to reinitialize walkers according to some scheme (Jo'j)
%emcee(nwalk, nmcmc; scale=0.01, gaussian, priors=priors, output=fmcmc, sim_count=50, serial);

% Exit ISIS
quit;
