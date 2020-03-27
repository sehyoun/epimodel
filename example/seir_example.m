%% SEIR model

% Add location to where epimodel is
%   or move <epimodel.m> file into the work folder
addpath('../src/')

%% SEIR model
% Building the model
seir.exposed.infectious = 35;
seir.exposed.recovered = 15;
seir.infectious.recovered = 15;
seir2.susceptible.infectious.exposed = 25;


seir_model = epimodel(seir, seir2);

% Set initial distribution
init_dist.infectious = 1e-6;
init_dist.susceptible = 1-init_dist.infectious;
seir_model.set_initial_dist(init_dist);

% Simulate
seir_model.simulate(15);

% The simulation results are in 'results'
figure(2);
subplot(1,4,1);
plot(seir_model.time_knots, seir_model.results.susceptible);
title('susceptible');
subplot(1,4,2);
plot(seir_model.time_knots, seir_model.results.exposed);
title('exposed');
subplot(1,4,3);
plot(seir_model.time_knots, seir_model.results.infectious);
title('infectious');
subplot(1,4,4);
plot(seir_model.time_knots, seir_model.results.recovered);
title('recovered');
