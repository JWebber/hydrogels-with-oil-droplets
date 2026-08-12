function particle_paths = agent_based_diffusion(particle_starts, timesteps, K, unit_size, drop_radius, n_units)
    % ---
    % SET UP THE GRID PATTERN
    % ---
    if(unit_size/2 == round(unit_size/2)) % If even grid size, add 1 to make it odd (so there is a centre)
        unit_size = unit_size + 1;
    end
    if (drop_radius*2 >= unit_size) % Check on size of droplet
        error('The droplet must be smaller than the cell size');
    end

    % Decide which points are in the droplet and which are in the water
    dropPoints = get_droplet_points(unit_size, drop_radius);
    gridPoints = zeros(unit_size^2, 2);

    for i = 1 : unit_size
        gridPoints((i-1)*unit_size+1:i*unit_size, :) = [repmat(i, [unit_size, 1]), (1:unit_size)'];
    end
    
    waterPoints = setdiff(gridPoints, dropPoints, 'rows');
    
    % Now create the full grid to hold all points
    grid = zeros(unit_size*n_units, unit_size*n_units);

    % ---
    % SET UP THE TIMESTEPPER
    % ---
    n_particles = size(particle_starts, 1);
    particle_pos = particle_starts;
    particle_paths = ones(timesteps, n_particles, 2);

    rnd = rand(timesteps, n_particles, 2);

    % ---
    % RUN THE LOOP
    % ---
    for ts = 1 : timesteps

        % Initialise grid
        grid = zeros(unit_size*n_units, unit_size*n_units);
        for p = 1 : n_particles
            if(~isnan(particle_pos(p, 1)) && ~isnan(particle_pos(p, 2)))
                grid(particle_pos(p, 1), particle_pos(p, 2)) = grid(particle_pos(p, 1), particle_pos(p, 2)) + 1;
            end
        end 

        % Loop over particles
        for p = 1 : n_particles
            % Has the particle left the grid?
            if(isnan(particle_pos(p, 1)) || isnan(particle_pos(p, 2)))
                % Do nothing
            else
                % Get the particle's position in its grid
                grid_pos = [mod(particle_pos(p, 1), unit_size), mod(particle_pos(p, 2), unit_size)];
                X = floor(particle_pos(p, 1)/unit_size);
                Y = floor(particle_pos(p, 2)/unit_size);
    
                % Where is water and where is droplet in the subunit?
                subunitWaterPoints = [X, Y]*unit_size + waterPoints;
                subunitDropPoints = [X, Y]*unit_size + dropPoints;
    
                % Where are the other solute particles here?
                grWPs = grid(sub2ind([unit_size*n_units, unit_size*n_units], subunitWaterPoints(:, 1), subunitWaterPoints(:, 2)));
                grDPs = grid(sub2ind([unit_size*n_units, unit_size*n_units], subunitDropPoints(:, 1), subunitDropPoints(:, 2)));
    
                % Get the local concentration in water and oil
                conc_water = sum(grWPs(:))/size(waterPoints, 1);
                conc_oil = sum(grDPs(:))/size(dropPoints, 1);

                adjacencies = [0, 1; 0, -1; 1, 0; -1, 0]; %; 1, 1; 1, -1; -1, -1; -1, 1];
                adjPoints = grid_pos + adjacencies;


                pWO = 0.5*(1+ tanh(1-(2*conc_oil + 1/size(dropPoints, 1))/(K*(2*conc_water - 1/size(waterPoints, 1)))));
                pOW = 1-0.5*(1+ tanh(1-(2*conc_oil - 1/size(dropPoints, 1))/(K*(2*conc_water + 1/size(waterPoints, 1)))));

                if(~ismember(grid_pos, dropPoints, 'rows')) % Starts in water
                    % Are we adjacent to oil?
                    nearby_oil = ismember(adjPoints, dropPoints, 'rows');

                    if(any(nearby_oil) && rnd(ts, p, 1) <= pWO)
                        potential_locations = adjacencies(nearby_oil, :);
                        particle_pos(p, :) = particle_pos(p, :) + potential_locations(ceil(rnd(ts, p, 2)*sum(nearby_oil)), :);
                    else
                        potential_locations = adjacencies(~nearby_oil, :);
                        particle_pos(p, :) = particle_pos(p, :) + potential_locations(ceil(rnd(ts, p, 2)*sum(~nearby_oil)), :);
                    end
                else % Starts in oil
                    % Are we adjacent to water?
                    nearby_water = ~ismember(adjPoints, dropPoints, 'rows');

                    if(any(nearby_water) && rnd(ts, p, 1) <= pOW)
                        potential_locations = adjacencies(nearby_water, :);
                        particle_pos(p, :) = particle_pos(p, :) + potential_locations(ceil(rnd(ts, p, 2)*sum(nearby_water)), :);
                    else
                        potential_locations = adjacencies(~nearby_water, :);
                        particle_pos(p, :) = particle_pos(p, :) + potential_locations(ceil(rnd(ts, p, 2)*sum(~nearby_water)), :);
                    end
                end
    
                % Leaving the grid? Just set position to NaN then.
                if(particle_pos(p, 1) == 0 || particle_pos(p, 1) == unit_size*n_units || particle_pos(p, 2) == 0 || particle_pos(p, 2) == unit_size*n_units)
                    particle_pos(p, 1) = NaN;
                    particle_pos(p, 2) = NaN;
                end
            end
        end

        particle_paths(ts, :, :) = particle_pos;
    end
end

