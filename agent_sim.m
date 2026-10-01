function [Deff, allpos, means] = agent_sim(K, N, R, cell_size, Ncells, steps, window)
    % ---
    % AGENT SIMULATION
    % K = partition coefficient
    % N = number of particles
    % R = radius of oil droplets
    % cell_size = size of each repeated cell with a droplet at its centre
    % Ncells = number of cells in each direction (i.e. Ncells x Ncells
    %   cells)
    % steps = number of timesteps
    % window = window over which to average the calculation of Deff
    % ---

    % If number of cells is even, throw an error; we need a well-defined midpoint
    if(mod(Ncells, 2) == 0 )
        error('The grid size needs to be odd.');
    end

    mp = (cell_size/2)*[1 1]; % The midpoint of each cell

    % We assume that the drops don't get too close to the edge of the cell
    if (R >= (cell_size/2) - 1 || sqrt(2*pi)*R >= cell_size)
        error('The droplets come too close to the cell boundary.');
    end

    rng('shuffle');

    % INITIALISATION
    % Start with solute dispersed randomly through a box enclosing the oil
    % droplet with equal areas water and oil.
    box_size = sqrt(2*pi)*R;
    positions = -(box_size/2) + box_size*rand(N, 2); % Initialise start positions
    gridpositions = ((Ncells+1)/2)*ones(N, 2); % Place all particles in the middle cell

    % Return variables
    allpos = NaN*zeros(steps, N, 2); % Hold all positions
    means = NaN*zeros(steps, 1); % Hold all mean radii

    % Oil and water concentrations
    cw = zeros(Ncells, Ncells);
    co = zeros(Ncells, Ncells);
    vo = pi*R^2;
    vw = cell_size^2 - vo;

    % Initialise the oil and water concentrations
    for i = 1 : N
        p = positions(i, :);
        within_cell_pos = [mod(p(1)+(cell_size/2), cell_size), mod(p(2)+(cell_size/2), cell_size)];
        r = norm(within_cell_pos-mp);

        if(r <= R)
            co((Ncells+1)/2, (Ncells+1)/2) = co((Ncells+1)/2, (Ncells+1)/2) + 1/vo;
        else
            cw((Ncells+1)/2, (Ncells+1)/2) = cw((Ncells+1)/2, (Ncells+1)/2) + 1/vw;
        end
    end

    % Random seeds for the simulation
    jumpRandoms = rand(steps, N);
    stepRandoms = rand(steps, N);
    jumpPosRandoms = rand(steps, N, 2);

    % Main loop
    for t = 1 : steps
        allpos(t, :, :) = positions;
        means(t) = mean(sqrt(positions(:, 1).^2+positions(:, 2).^2), 'omitnan');

        for i = 1 : N
            % Where is the current particle?
            p = positions(i, :);
            grid_point = gridpositions(i, :);

            % Only timestep if we've not left the grid
            if(~isnan(p))
                within_cell_pos = [mod(p(1)+(cell_size/2), cell_size), mod(p(2)+(cell_size/2), cell_size)];
                old_within = within_cell_pos;
                r = norm(within_cell_pos-mp);

                jumped = 0;
    
                % Start with redistribution
                if (r > R)
                    % Base jump probabilities off an average of initial and
                    % final states or else we'll run into problems when co
                    % or cw = 0.
                    mean_co = co(grid_point(1), grid_point(2)) + 0.5/vo;
                    mean_cw = max(0, cw(grid_point(1), grid_point(2)) - 0.5/vw);
                    jump = PWO(K, mean_co, mean_cw)*(r<(R+1)); % ONLY jump if adjacent to droplet

                    if (jumpRandoms(t, i) < jump)
                        within_cell_pos = mp + R*jumpPosRandoms(t, i, 1)*[cos(2*pi*jumpPosRandoms(t, i, 2)), sin(2*pi*jumpPosRandoms(t, i, 2))];
                        co(grid_point(1), grid_point(2)) = co(grid_point(1), grid_point(2)) + 1/vo;
                        cw(grid_point(1), grid_point(2)) = cw(grid_point(1), grid_point(2)) - 1/vw;
                        jumped = 1;
                    end
                else
                    % Base jump probabilities off an average of initial and
                    % final states or else we'll run into problems when co
                    % or cw = 0.
                    mean_co = max(0, co(grid_point(1), grid_point(2)) - 0.5/vo);
                    mean_cw = cw(grid_point(1), grid_point(2)) + 0.5/vw;
                    jump = (1-PWO(K, mean_co, mean_cw))*(r>(R-1)); % ONLY jump if adjacent to water

                    if (jumpRandoms(t, i) < jump)
                        within_cell_pos = mp + (R+jumpPosRandoms(t, i, 1))*[cos(2*pi*jumpPosRandoms(t, i, 2)), sin(2*pi*jumpPosRandoms(t, i, 2))];
                        co(grid_point(1), grid_point(2)) = co(grid_point(1), grid_point(2)) - 1/vo;
                        cw(grid_point(1), grid_point(2)) = cw(grid_point(1), grid_point(2)) + 1/vw;
                        jumped = 1;
                    end
                end

                new_grid = grid_point;
                p = p + (within_cell_pos-old_within);

                if (jumped < 1)
                    theta = 2*pi*stepRandoms(t, i);
                    movedcell = [floor((within_cell_pos(1) + cos(theta))/cell_size), floor((within_cell_pos(2) + sin(theta))/cell_size)];
                    
                    pN = p + [cos(theta), sin(theta)];

                    % Deal with cases where we might diffuse between
                    % phases. Stop these by reducing the jump distance
                    % (phase boundaries are hard).
                    if(r < R)
                        if(norm(within_cell_pos+[cos(theta), sin(theta)]-mp) > R)
                            pN = p + (R-r)*[cos(theta), sin(theta)]/2;
                        end
                    else
                        if(norm(within_cell_pos+[cos(theta), sin(theta)]-mp) < R)
                            pN = p + (r-R)*[cos(theta), sin(theta)]/2;
                        end
                    end

                    p = pN;
    
                    new_grid = grid_point + movedcell;
    
                    if(norm(movedcell) ~= 0)
                        cw(grid_point(1), grid_point(2)) = cw(grid_point(1), grid_point(2)) - 1/vw;
    
                        if (new_grid(1) < 1 || new_grid(1) > Ncells || new_grid(2) < 1 || new_grid(2) > Ncells)
                            p = NaN*[1 1];
                            fprintf('PARTICLE LEFT GRID\n');
                        else
                            cw(new_grid(1), new_grid(2)) = cw(new_grid(1), new_grid(2)) + 1/vw;
                        end
                    end
                end

                positions(i, :) = p;
                gridpositions(i, :) = new_grid;
            end
        end

        
        
    end
    
    % Compute the running effective diffusivity
    p = polyfit(log(steps-window:steps), log(means(steps-window:steps)), 1);
    pow = p(1);

    fitfun = @(D, t)(sqrt(4*D*t));
    ops = optimset('display', 'none');
    fitcurve = linspace(1/10, 10, window+1);
    Deff = lsqnonlin(@(x)(fitcurve.*(fitfun(x, steps-window:steps) - means(steps-window:steps))), 1, 0, [], [], [], [], [], [], ops);


    fprintf(strcat("Fit exponent is ", num2str(pow), ", should be 0.5. Effective diffusivity is ", num2str(Deff), '.\n'));

    % Put the 'confidence' in the output for Deff
    Deff = [Deff, max(0, 1-8*abs(pow-1/2))];
end

function prob = PWO(K, co, cw)
    prob = 1 - co/(K*cw + co);
end