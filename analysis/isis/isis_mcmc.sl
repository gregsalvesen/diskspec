% Run this script from the command line as follows:
% >> isis script.sl

% Data filenames
variable dir   = "/Users/salvesen/research/fcolabs/data/";
variable fdata = dir + "source_annulus_grp10.pha";
variable fback = dir + "back.pha";
variable frmf  = dir + "swxwt0to2s0_20010101v012.rmf";
variable farf  = dir + "xrt_annulus.arf";

% Energy range for fitting
variable Emin  = 0.5;   % [keV]
variable Emax  = 10.0;  % [keV]

% GRO J1655-40 parameters
variable nH=1.0;  variable nH_min=0.0;  variable nH_max=100000.0;
variable eta=0.0; variable eta_min=0.0; variable eta_max=1.0;
variable a=0.92;  variable a_min=-1.0;  variable a_max=0.9999;
variable i=68.65; variable i_min=0.0;   variable i_max=85.0;
variable Mbh=5.4; variable Mbh_min=0.0; variable Mbh_max=100.0;
variable Mdd=1.0; variable Mdd_min=0.0; variable Mdd_max=1000.0;
variable Dbh=3.2; variable Dbh_min=0.0; variable Dbh_max=100.0;
variable hd=1.7;  variable hd_min=1.0;  variable hd_max=10.0;
variable rflag=1;
variable lflag=1;
variable norm=1.0;

% Output parameter file
variable fpars = dir + "parameter_file.par";

% GRO J1655-40 known parametes and their 1-sigma uncertainties
variable mean_a=0.92;   variable sigma_a=0.02;
variable mean_i=68.65;  variable sigma_i=1.50;
variable mean_Mbh=5.40; variable sigma_Mbh=0.3;
variable mean_Dbh=3.20; variable sigma_Dbh=0.2;

% MCMC inputs
variable nw    = 100;   % Number of walkers per free parameter
variable nsim  = 1000;  % Number of iterations, or steps, for each walker
variable fmcmc = dir + "emcee-chain.fits";

% Load ISISscripts
variable isisscripts="/Users/salvesen/soft/isis-scripts/share/isisscripts";
require(isisscripts);

% The RMF file does not adhere closely to the OGIP standard format.
% The workaround is to set the global variable Rmf_OGIP_Compliance=0.
% This reduces the required level of standards compliance for all input RMF files.
Rmf_OGIP_Compliance=0;

% Load the data...and the associated BACKFILE, RESPFILE, ANCRFILE
variable data_id = load_data(fdata);

% We could load the background, RMF, and ARF explicitly, then assign them to the data set
% Note: The "() = " syntax allows us to ignore the return value
%define_back(data_id, fback);
%variable rmf_id = load_rmf(frmf);
%variable arf_id = load_arf(farf);
%assign_rmf(data_id, rmf_id);
%assign_arf(data_id, arf_id);

% Group the data (20 counts/bin equates to a S/N of 4.47 --> 4.47 = 20/sqrt(20))
% Careful: ISIS ignores any default binning specified in the FITS header unless Isis_Use_PHA_Grouping=1.
%Isis_Use_PHA_Grouping=1;
group(data_id; min_sn=4.47, bounds=Emin, unit="keV");

% Only notice energies in the range 0.5-10 keV
notice_values(data_id, Emin, Emax; unit="keV");

% Define the fit function as TBabs*kerrbb
fit_fun("TBabs(1)*kerrbb(1)");

% Set the parameters: value, freeze, min, max (Note: 0=thaw, 1=freeze)
variable thw=0;  % Thaw
variable frz=1;  % Freeze
set_par("TBabs(1).nH",   nH,  thw,   nH_min,  nH_max);
set_par("kerrbb(1).eta", eta, frz, eta_min, eta_max);
set_par("kerrbb(1).a",   a,   frz, a_min,   a_max);
set_par("kerrbb(1).i",   i,   frz, i_min,   i_max);
set_par("kerrbb(1).Mbh", Mbh, frz, Mbh_min, Mbh_max);
set_par("kerrbb(1).Mdd", Mdd, thw,   Mdd_min, Mdd_max);
set_par("kerrbb(1).Dbh", Dbh, frz, Dbh_min, Dbh_max);
set_par("kerrbb(1).hd",  hd,  thw,   hd_min,  hd_max);
set_par("kerrbb(1).rflag", rflag, frz);
set_par("kerrbb(1).lflag", lflag, frz);
set_par("kerrbb(1).norm",  norm,  frz);

% Renorm
%() = renorm_counts;

% Fit
() = fit_counts;

% Collect the best fit parameters
%nH  = get_par("TBabs(1).nH");
%Mdd = get_par("kerrbb(1).Mdd");
%hd  = get_par("kerrbb(1).hd");

% Collect parameter uncertainties (level=0 for 68% confidence)
%variable nH_lo;  variable nH_hi;
%variable Mdd_lo; variable Mdd_hi;
%variable hd_lo;  variable hd_hi;
%(nH_lo,  nH_hi)  = vconf("TBabs(1).nH", 0);
%(Mdd_lo, Mdd_hi) = vconf("kerrbb(1).Mdd", 0);
%(hd_lo,  hd_hi)  = vconf("kerrbb(1).hd", 0);

% Thaw the static system parameters and the normalization
thaw("kerrbb(1).a");
thaw("kerrbb(1).i");
thaw("kerrbb(1).Mbh");
thaw("kerrbb(1).Dbh");
thaw("kerrbb(1).norm");

% Save the fit parameters
save_par(fpars);

% Load the fit parameters
%load_par(fpars);

%====================================================================================================

% Fit using emcee (start with uniform initial walkers at the moment)
variable parlist   = ["kerrbb(1).a", "kerrbb(1).i", "kerrbb(1).Mbh", "kerrbb(1).Dbh"];
variable prior_a   = struct{mean=mean_a, sigma=sigma_a};
variable prior_i   = struct{mean=mean_i, sigma=sigma_i};
variable prior_Mbh = struct{mean=mean_Mbh, sigma=sigma_Mbh};
variable prior_Dbh = struct{mean=mean_Dbh, sigma=sigma_Dbh};
variable priorlist = [prior_a, prior_i, prior_Mbh, prior_Dbh];
variable priors    = struct{parlist=parlist, priorlist=priorlist};
emcee(nw, nsim; scale=0.1, gaussian, priors=priors, output=fmcmc, serial);

% Use the minimum chi^2 best fit parameters and uncertinties to initialize the MCMC walkers
% Use the current parameter values (not the output fpars file, which has some things frozen still)
% Specify init_chain=fwalk0
% Do the MCMC analysis with UNIFORM priors and with GAUSSIAN priors. Compare!

% Exit ISIS
quit;
