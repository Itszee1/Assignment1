

Nm = 10;     % Number of Molecules
colours = ["b", "c", "g", "k", "m", "r", "y"];   %This are the different colors that will be represented in the plot

kb = 1.38064852; %This is the Boltzmann constant
T = 300;         %This is Temperature in Kelvin
m = 9.10938356E-31; %This is mass of an Electron


TimeCollision = 0.2;   % This is measured in Picoseconds
MFP = 0;          % This is the Mean Free Path 
timesteps = 1E-3;           
scattProb = 1 - exp(-timesteps/TimeCollision);  %This is the scattering Probability

% Position of the molecule
molecules(:, 2) = molecules(:, 2)*200;  % x-coordinates
molecules(:, 1) = molecules(:, 1)*100;  % y-coordinates
molecules = rand(Nm, 2);


% column 3 is x-velocity,while column 4 is x-velocity
molecules(:, 3:4) = randn(Nm, 2) *  sqrt(kb * T / m) / 1E15;  % scaled to femtoseconds
% subplot(2, 1, 1)  % subplots don't look good in publish

figure(1)
title(sprintf("Initial Avg Velocity: %s", mean(sqrt(molecules(:, 3).^2 + molecules(:, 4).^2))))
histogram(rssq(molecules(:, 3:4), 50));

position_of_last_scatter = molecules(:, 3:4);



% main time loop
for i = 0:50      
    previous_particles = molecules;
    
    % update positions
    molecules(:, 1) = molecules(:, 1) + molecules(:, 4);
    molecules(:, 2) = molecules(:, 2) + molecules(:, 3);
    
     % This is to make sure the particles do not leave the shape from the
    % right side
    x_boundary_changes_right = molecules(:, 2) > 200;
    if any(x_boundary_changes_right)
        molecules(:, 2) = molecules(:, 2) .* ~x_boundary_changes_right;
    end
 %This is to make sure the particles do not leave the shape from the
    % right side
    x_boundary_changes_left = molecules(:, 2) < 0;
    if any(x_boundary_changes_left)
        molecules(:, 2) = molecules(:, 2) + 200 * x_boundary_changes_left - abs(molecules(:, 2) .* x_boundary_changes_left);
    end
    %This is to make sure the particles do not leave the shape from the top 
     %of shape
    y_boundary_changes_upper = molecules(:, 1) > 100;
    if any(y_boundary_changes_upper)
        molecules(:, 4) = molecules(:, 4) - (2 * molecules(:, 4) .* y_boundary_changes_upper);
        overshoot = (molecules(:, 1) - 100) .* y_boundary_changes_upper;
        molecules(:, 1) = molecules(:, 1) - 2 * overshoot;
    end
     %This is to make sure the particles do not leave the bottom of the
    %shape
    y_boundary_changes_lower = molecules(:, 1) < 0;
    if any(y_boundary_changes_lower)
        molecules(:, 4) = molecules(:, 4) - (2 * molecules(:, 4) .* y_boundary_changes_lower);
        overshoot = abs(molecules(:, 1)) .* y_boundary_changes_lower;
        molecules(:, 1) = molecules(:, 1) + 2 * overshoot;
    end
    
    % Model of Scattering electrons
    Scattered = rand(Nm, 1) < scattProb;
    if any(Scattered)
        Scattered(:, 2) = Scattered;
        FP = abs(position_of_last_scatter - Scattered(:, 2) .* molecules(:, 3:4));%This is free path
        MFP = mean(FP, "all");
        angles = randn(Nm, 1) .* 2 * pi;
        rethermalized_velocities(:, 1:2) = randn(Nm, 2) * sqrt(kb * T / m) / 1E15;
        rethermalized_velocities = rethermalized_velocities .* Scattered;
        molecules(:, 3:4) = molecules(:, 3:4) .* ~Scattered + rethermalized_velocities;
        position_of_last_scatter = molecules(:, 3:4);
    end
    
    % Plot of Molecules 
    x_boundary_affected = x_boundary_changes_right | x_boundary_changes_left;
    temp_avg = mean(((sqrt(molecules(:, 3).^2 + molecules(:, 4).^2) .* 1E15).^2) .* m ./ kb);
    % subplot(2, 1, 2)
    hold on
    figure(2)
    title(sprintf("Avg Temperature: %s, Mean Free Path is: %s", temp_avg, MFP))
    for i = 1:length(molecules)
        if ~x_boundary_affected(i)
            plot([previous_particles(i, 2) molecules(i, 2)], [previous_particles(i, 1) molecules(i, 1)], colours(mod(i, length(colours)) + 1))
        end
    end
    axis([0 200 0 100])
    pause(0.1)
end
