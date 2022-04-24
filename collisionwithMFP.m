

Nm = 20;     % Number of Molecules
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

oldMolecules = molecules(:, 3:4);



% main time loop
for i = 0:50      
    prMolecules = molecules;
    
    % update positions
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
    
    % Model of Scattering electrons
    Scattered = rand(Nm, 1) < scattProb;
    if any(Scattered)
        Scattered(:, 2) = Scattered;
        FP = abs(oldMolecules - Scattered(:, 2) .* molecules(:, 3:4));%This is free path
        MFP = mean(FP, "all");
        angles = randn(Nm, 1) .* 2 * pi;
        rethermVelo(:, 1:2) = randn(Nm, 2) * sqrt(kb * T / m) / 1E15; %This is rethermalized velocities
        rethermVelo = rethermVelo .* Scattered;
        molecules(:, 3:4) = molecules(:, 3:4) .* ~Scattered + rethermVelo;
        oldMolecules = molecules(:, 3:4);
    end
    
    % Plot of Molecules 
    Xchanges = XchangesRight | XchangesLeft;
    avgTemp = mean(((sqrt(molecules(:, 3).^2 + molecules(:, 4).^2) .* 1E15).^2) .* m ./ kb);
    % subplot(2, 1, 2)
    hold on
    figure(2)
    title(sprintf("Average Temperature: %s, The Mean Free Path is: %s", avgTemp, MFP))
    for i = 1:length(molecules)
        if ~Xchanges(i)
            plot([prMolecules(i, 2) molecules(i, 2)], [prMolecules(i, 1) molecules(i, 1)], colours(mod(i, length(colours)) + 1))
        end
    end
    axis([0 200 0 100])
    pause(0.1)
end
