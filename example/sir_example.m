%% SIR Model

% Add location to where epimodel is
%   or move <epimodel.m> file into the work folder
addpath('../src/');

%% SIR model
% Building the model
sir.infectious.recovered = 100;
sir2.susceptible.infectious.infectious = 105;
sir_model = epimodel(sir, sir2);

% Set Initial distribution
% If you do not fill in the values, it defaults to 0
init_dist.infectious = 1e-6;
init_dist.susceptible = 1 - init_dist.infectious;
sir_model.set_initial_dist(init_dist);

% Simulate. This is (end_time)
% Note: There is a time "discretization" parameter that one can set.
sir_model.simulate(10);

% The simulation results are in 'results'
figure(1);
subplot(1,3,1);
plot(sir_model.time_knots, sir_model.results.susceptible);
title('susceptible');
subplot(1,3,2);
plot(sir_model.time_knots, sir_model.results.infectious);
title('infectious');
subplot(1,3,3);
plot(sir_model.time_knots, sir_model.results.recovered);
title('recovered');

