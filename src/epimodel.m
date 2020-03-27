classdef epimodel < handle
  % by SeHyoun Ahn, March 2020
  % BSD 2-Clause
  %

  properties
    simulated
    time_knots

    iii
    jjj
    kkk
    vvv

    name2loc
    loc2name
    n_states

    linear_dynamics
    init_dist

    results
  end

  methods
    function obj = epimodel(struct_linear, struct_interaction)
      % Parses struct into matrices
      %
      % Args:
      %   struct_linear (struct): struct corresponding to linear flows. The syntax is ::
      %
      %    struct.origin.target = flowrate
      %
      %   struct_interaction (struct): struct corresponding to interacted flows. The syntax is ::
      %
      %    struct.origin.interaction.target = flowrate
      %
      % Example:
      %   ::
      %
      %     % SEIR model
      %     struct_linear.exposed.infectious = 1;
      %     struct_linear.infectious.recovered = 1;
      %     struct_interaction.susceptible.infectious.exposed = 1;
      %     seir_model = epimodel(struct_linear, struct_interaction);
      %

      ii = [];
      jj = [];
      vv = [];
      fields_left = fieldnames(struct_linear);
      loc_curr = 1;
      obj.name2loc = struct();
      obj.loc2name = cell(50,1);
      for iter_left = 1:length(fields_left)
        field_left = fields_left{iter_left};
        if ~isfield(obj.name2loc, field_left)
          obj.name2loc.(field_left) = loc_curr;
          obj.loc2name{loc_curr} = field_left;
          loc_curr = loc_curr + 1;
        end
        inner_struct = struct_linear.(field_left);
        fields_right = fieldnames(inner_struct);
        for iter_right = 1:length(fields_right)
          field_right = fields_right{iter_right};
          if ~isfield(obj.name2loc, field_right)
            obj.name2loc.(field_right) = loc_curr;
            obj.loc2name{loc_curr} = field_right;
            loc_curr = loc_curr + 1;
          end
          ii = [ii; obj.name2loc.(field_left)];
          jj = [jj; obj.name2loc.(field_right)];
          vv = [vv; inner_struct.(field_right)];
        end
      end

      iii = [];
      jjj = [];
      kkk = [];
      vvv = [];
      fields_left = fieldnames(struct_interaction);
      for iter_left = 1:length(fields_left)
        field_left = fields_left{iter_left};
        if ~isfield(obj.name2loc, field_left)
          obj.name2loc.(field_left) = loc_curr;
          obj.loc2name{loc_curr} = field_left;
          loc_curr = loc_curr + 1;
        end
        inner_struct = struct_interaction.(field_left);
        fields_right = fieldnames(inner_struct);
        for iter_right = 1:length(fields_right)
          field_right = fields_right{iter_right};
          if ~isfield(obj.name2loc, field_right)
            obj.name2loc.(field_right) = loc_curr;
            obj.loc2name{loc_curr} = field_right;
            loc_curr = loc_curr + 1;
          end
          val_right = inner_struct.(field_right);
          fields_inner = fieldnames(val_right);
          for iter_right = 1:length(fields_inner)
            field_inner = fields_inner{iter_right};
            if ~isfield(obj.name2loc, field_inner)
              obj.name2loc.(field_inner) = loc_curr;
              obj.loc2name{loc_curr} = field_inner;
              loc_curr = loc_curr + 1;
            end
            iii = [iii; obj.name2loc.(field_left)];
            jjj = [jjj; obj.name2loc.(field_inner)];
            kkk = [kkk obj.name2loc.(field_right)];
            vvv = [vvv; val_right.(field_inner)];
          end
        end
      end

      loc_curr = loc_curr - 1;
      obj.n_states = loc_curr;
      obj.loc2name = obj.loc2name(1:loc_curr);
      obj.linear_dynamics = sparse(jj, ii, vv, loc_curr, loc_curr);

      obj.iii = iii;
      obj.jjj = jjj;
      obj.kkk = kkk;
      obj.vvv = vvv;
    end

    function set_initial_dist(obj, struct_dist)
      % Sets initial distribution
      %
      % Args:
      %   struct_dist (struct): struct of initial distribution. Unfilled values default to zero
      %
      % Example:
      %   ::
      %
      %     struct_dist.infectious = 1e-4;
      %     struct_dist.susceptible = 1 - 1e-4;
      %     seir_model.set_initial_dist(struct_dist);
      %


      obj.init_dist = zeros(obj.n_states, 1);
      fields = fieldnames(struct_dist);
      for iter = 1:length(fields)
        field = fields{iter};
        obj.init_dist(obj.name2loc.(field)) = struct_dist.(field);
      end
    end

    function simulate(obj, end_time, time_step)
      % Runs simulation/ODE forward
      %
      % Args:
      %   end_time(double): end time of simulation
      %   time_step(double): [default 1e-3] time step of discretization
      %
      % Example:
      %   ::
      %
      %     seir_model.simulate(10);
      %
      % See Also:
      %   * :attr:`results <epimodel.results>`
      %   * :attr:`set_initial_dist <epimodel.epimodel.set_initial_dist>`
      %
      % Note:
      %   Note that The values are computed using an "explicit" udpate, so `time_step` needs to small enough to satisfy the CFL condition.
      %


      if nargin < 3
        time_step = 1e-3;
      end
      n_simul = ceil(end_time/time_step);
      obj.time_knots = (0:n_simul).*time_step;
      A = obj.linear_dynamics.*time_step;

      x_curr = obj.init_dist;
      obj.simulated = zeros(obj.n_states, n_simul+1);
      obj.simulated(:, 1) = x_curr;
      for iter = 2:n_simul+1
        A_curr = A + time_step.*sparse(obj.jjj, obj.iii, obj.vvv.*x_curr(obj.kkk), obj.n_states, obj.n_states);

        A_diag = sum(A_curr);
        A_curr = A_curr - spdiags(A_diag(:), 0, obj.n_states, obj.n_states);

        x_curr = x_curr + (A_curr)*x_curr;
        obj.simulated(:, iter) = x_curr;
      end
      fields = fieldnames(obj.name2loc);
      for iter_field = 1:length(fields)
        obj.results.(fields{iter_field}) = obj.simulated(obj.name2loc.(fields{iter_field}), :)';
      end
      % obj.results.time_knots = obj.time_knots;
    end
  end
end
