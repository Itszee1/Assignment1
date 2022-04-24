

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
    XchangesRight = molecules(:, 2) > 200;
    if any(XchangesRight)
        molecules(:, 2) = molecules(:, 2) .* ~XchangesRight;
    end
    %This is to make sure the particles do not leave the shape from the
    % right side
    XchangesLeft = molecules(:, 2) < 0;
    if any(XchangesLeft)
        molecules(:, 2) = molecules(:, 2) + 200 * XchangesLeft - abs(molecules(:, 2) .* XchangesLeft);
    end
     %This is to make sure the particles do not leave the shape from the top 
     %of shape
    YchangesUp = molecules(:, 1) > 100;
    if any(YchangesUp)
        molecules(:, 4) = molecules(:, 4) - (2 * molecules(:, 4) .* YchangesUp);
        overshoot = (molecules(:, 1) - 100) .* YchangesUp;
        molecules(:, 1) = molecules(:, 1) - 2 * overshoot;
    end
    %This is to make sure the particles do not leave the bottom of the
    %shape
    YchangesDown = molecules(:, 1) < 0;
    if any(YchangesDown)
        molecules(:, 4) = molecules(:, 4) - (2 * molecules(:, 4) .* YchangesDown);
        overshoot = abs(molecules(:, 1)) .* YchangesDown;
        molecules(:, 1) = molecules(:, 1) + 2 * overshoot;
    end
    
    % plot position updates of particles 
    Xchanges = XchangesRight | XchangesLeft;
    AvgTemp = mean(((sqrt(molecules(:, 3).^2 + molecules(:, 4).^2) .* 1E15).^2) .* m ./ kb);
    title(sprintf("Average Temperature: %s", AvgTemp))
    
    for i = 1:length(molecules)
        if ~Xchanges(i)
            plot([OldMolecules(i, 2) molecules(i, 2)], [OldMolecules(i, 1) molecules(i, 1)], colors(mod(i, length(colors)) + 1))
        end
    end
    axis([0 200 0 100])
    hold on
    pause(0.1)
end