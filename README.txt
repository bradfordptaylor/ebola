README:

This document explains how to recreate the figures in BP Taylor et al. (2015) from the provided MATLAB code. Short descriptions of all dependent functions are provided. All model simulations are stochastic utilizing the gillespie algorithm. Sampling of the dynamics are done in regular intervals that start at the beginning of dynamics (e.g., weekly case counts for country data). In this readme, function dependencies are listed only once for each occurrence of the relevant data. Function dependencies are presented in a nested fashion to display subdependencies.

Figure 1:
Main file: 
plot_fig1.m -- automatically saves the 3 subfigures composing figure 1 in subdirectory 'figures'.

Dependent data: 
fig1_data.mat and allcasecnts.mat in subdirectory 'data' as created by process_fig1.m

Dependent Files: 
process_fig1.m -- runs the 10^4 simulations of the SEIRD model with a trigger population as explained in the text. outputs: fig1_data.mat. Also runs 5 trajectories of the model to recreate the left panel. outputs: allcasecnts.mat

	SEIRD_gill_fig1.m -- outputs the model dynamics for 42 days after a trigger population is hit

		getrxns.m -- outputs vector of changes to population structure given a list of reactions defined by the model.

		gilll_onestep.m -- runs one step of the gillespie algorithm. outputs the updated time after the reaction and the updated population.

			SEIRD.m -- defines the rate of each reaction in the SEIRD model given the current population.


	SEIRD_gill_fig1_traj.m -- simulates the model and outputs the full dynamics of up to 42 days after a defined trigger population is reached.
	
fitdata_pois.m -- outputs the best fit and standard error of the best fit exponential to data assuming poisson distributed errors.

diststats95.m -- outputs the median and upper and lower CI for an input distribution (vector)


Figure 2 --
Main File A: 
plot_simci_cum_trigger.m -- outputs the figure of the \tau_c CI for the simulated SEIRD dynamics with mock data with \tau = 20 days as explained in the text. This uses the cumulative case count.

Dependent data: 
ci_sim_cum_trig50.mat -- obtained from process_trigger.m

Dependent files:
process_trigger.m -- loads prior data of dynamics up to 42 days after a trigger population of 50 was reached. Curates and saves this data as 42 days of dynamics after varying trigger populations set at [10,20,30,40,50]. This is done to compare methods given different trigger populations. The 50 cumulative infected is used for Figure 2.

	Dependent data:
	simdata_final_L#_r#.mat--dynamics where parameters are specified by the first #  corresponding to a list of different characteristic growth rates of the underlying model. The r # refer to statistically independent ensembles of runs of the data. data obtained from cluster_simdyn_final.m

	Dependent files:
	cluster_simdyn_final.m -- given lambda and run number saves data of an ensemble of 100 runs. This is done to allow parallelization on a cluster.

		SEIRD_gill_fulldyn -- outputs the stochastic SEIRD dynamics sampled at daily intervals up until 42 days after a given trigger population set to 50.

			getrxns.m -- outputs vector of changes to population structure given a list of reactions defined by the model.

			gilll_onestep.m -- runs one step of the gillespie algorithm. outputs the updated time after the reaction and the updated population.

				SEIRD.m -- defines the rate of each reaction in the SEIRD model given the current population.

Main File B: plot_simciR0_cum_trigger.m -- outputs the figure of the R0 CI for the simulated SEIRD dynamics with mock data with \tau = 20 days as explained in the text. This uses the cumulative case count. It is completely dependent on the previous figure and the corresponding data and dependencies.

Figure 3 -- 
Main File A:
plot_country_R0direct.m -- automatically generates figures for the left hand column of Figure 3 by specifying the country. Input: 'guinea' = Guinea, 'sleone' = Sierra Leone 'liberia' = Liberia

Dependent data:
ci_country_R0direct.mat -- tau_c confidence intervals for different specified R0 as obtained by process_country_R0direct.m

Dependent files:
fitcasedata_pois_data.m -- fits exponential curves assuming poisson distributed errors to the cumulative case data from the specified country to output the measured inverse characteristic time.

	fitdata_pois_data.m -- fits exponential curves assuming poisson distributed errors to a given dynamics.

	Dependent data:
	casedata_gountries_gls.csv -- curated case data for all 3 countries obtained from the WHO

Dependent files (cont):
process_country_R0direct.m -- combines data of dynamics to obtain tau_c CI for various specified R0 for each country. The R0 considered are hardcoded as a vector and is analyzed depending on if the corresponding data exists in the 'data' subdirectory

	fitdata_pois.m -- outputs the best fit and standard error of the best fit exponential to data assuming poisson distributed errors.

	diststats95.m -- outputs the median and upper and lower CI for an input distribution (vector)

	Dependent Data: stochsim_country_R0#_r1.mat where # is the index for the R0 vector specifying the underlying R0 of the deterministic models. obtained from clusterfun_september.m

		Dependent Files:
		clusterfun_september.m -- Input: number specifying country, guinea=1,liberia=2,sleone=3; R0seq, index specifying R0 according to hardcoded vector; runseq = number to organize statistically independent ensembles of runs; numruns, number of runs in each saved ensemble. This file saves dynamics for stochastic simulations of SEIRD model with parameters specific to each country and the input underlying R0. A trigger population of 50 is hard-coded in for each country and the case data specifies what the trigger population is for the simulations.

			SEIRD_gill_datacum -- outputs dynamics for 42 days after the trigger population is specified.

				getrxns.m -- outputs vector of changes to population structure given a list of reactions defined by the model.

				gilll_onestep.m -- runs one step of the gillespie algorithm. outputs the updated time after the reaction and the updated population.

					SEIRD.m -- defines the rate of each reaction in the SEIRD model given the current population.

Main File B:
plot_country_projection_ccc.m -- automatically plots the case data for 42 days after the trigger population is reached and then projects forward based on the CI estimates. This function depends on data obtained from the previous part of the figure.

Dependent files:
R02lambda.m: given R0, rhoD, and other parameters not changed between sims this function gives the underlying theoretical characteristic time \tau_c.

Main File C:
plot_country_projection_icc.m -- Same as above but uses incident case data in the plot.


Figure S1:
Main File:
plot_trigger_compare.m -- automatically generates a plot comparing the R0 CI using our mock data but considering different trigger populations on the same datasets. All dependencies are contained in those of Fig 2

Figure S2:
Mainfile:
plot_sim_distcomp.m -- automatically generate plot comparing \tau_c distributions between fits to cumulative data and incident data. The data is the trigger 50 simulated data with the underlying characteristic time of 20 days. All dependencies are contained in those of Fig 2.

Figure S3:
Mainfile:
plot_errorhist.m -- automatically generate plot comparing the distributions of standard errors of \tau_c fits to cumulative data and incident data. The data is the trigger 50 simulated data with the underlying characteristic time of 20 days. All dependencies are contained in those of Fig 2.

Figure S4:
Mainfile:
figibm_lambda_bias_large -- automatically generates figure comparing the R0 identifiability problem as discussed in Weitz & Dushoff (2014). This function is self-contained given prior data figibm_lambda_bias_large.mat.
