# Compartment Model Helper

Due to the COVID-19, many people would write different iteration of the compartment model. This is a MATLAB class written to reduce the repetitive coding. For example, the SEIR model can be written as

```octave
linear.exposed.infectious = 1;
linear.infectious.recovered = 1;
interaction.susceptible.infectious.exposed = 1;
seir_model = epimodel(linear, interaction);

init_dist.infectious = 1e-4;
init_dist.susceptible = 1 - 1e-4;
seir_model.set_initial_dist(init_dist);

seir_model.simulate(10);

plot(seir_model.time_knots, seir_model.results.susceptible);
title('susceptible');
```

Check the [documentation](https://sehyoun.com/epimodel) for more details on syntax.

