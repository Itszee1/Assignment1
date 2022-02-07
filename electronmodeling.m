

Nm = 10;   %This is number of molecules
colors = ["b", "c", "g", "k", "m", "r", "y"];   %This are the different colors that will be represented in the plot

kb = 1.38064852; %This is the Boltzmann constant
T = 300;         %This is Temperature in Kelvin
m = 9.10938356E-31; %This is mass of an Electron
ThermVelo = sqrt(kb * T / m) / 1E15;  %This is the calculation for Thermal velocity

TimeCollision = 0.2;   % This is measured in Picoseconds
TimeSteps = 1E-3;            


% The positions of molecules
molecules = rand(Nm, 2);
molecules(:, 2) = molecules(:, 2)*200;  % x-coordinates
molecules(:, 1) = molecules(:, 1)*100;  % y-coordinates



% column 3 is x-velocity,while column 4 is -velocity
angles = randn(Nm, 1) .* 2 * pi;
molecules(:, 3) = randn(Nm, 1) + ThermVelo*cos(angles);
molecules(:, 4) = randn(Nm, 1) + ThermVelo*sin(angles);



% Loop for time 
for i = 0:50      
    OldMolecules = molecules;
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
    
    % plot position updates of particles 
    x_boundary_affected = x_boundary_changes_right | x_boundary_changes_left;
    AvgTemp = mean(((sqrt(molecules(:, 3).^2 + molecules(:, 4).^2) .* 1E15).^2) .* m ./ kb);
    title(sprintf("Average Temperature: %s", AvgTemp))
    
    for i = 1:length(molecules)
        if ~x_boundary_affected(i)
            plot([OldMolecules(i, 2) molecules(i, 2)], [OldMolecules(i, 1) molecules(i, 1)], colors(mod(i, length(colors)) + 1))
        end
    end
    axis([0 200 0 100])
    hold on
    pause(0.1)
end