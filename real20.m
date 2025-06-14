clear; clc;
%====================== INPUT DATA ==============================

r=0.30; % length of the crank PA in m
m=0.70; % length of the ink CD in m
p=0.55; % length of the distance OP in m
d=0.85; % length of the distance OB in m
s=0.65; % length of the distance CB in m
w=2; % angular velocity of the crank in rad/s
n=0.80; % length of the link EF, m
t=0.20; % length of the distance DE, m
L=2.0; % length of the distance L in m

%======================= COMPUTATION ==================================

% Define the range for x with steps of 0.001 radians
phi = (0:0.001:2*pi);

% define temporary variables and the function y
tmp1 = r * cos(phi);
tmp2 = r * sin(phi);
tmp3 = p - tmp1;
ang1 =  atan(tmp3./(d + tmp2));
ang2 = (pi/2) - ang1;
tmp4 = s * cos(ang2);
tmp5 = s * sin(ang2);
tmp6 = sqrt(m^2 - tmp5.^2);
tmp7 = L - t - tmp6 - tmp4 - d;

% Calculate the function y
yF = sqrt(n^2 - tmp7.^2);


%===================== Minimum and Maximum Values =====================

% Initialize variables to store maximum and minimum values
Maximum = yF(1); % Initialize variable
Minimum = yF(1);  % Initialize variable

% Iterate through yF values and find maximum and minimum values
for i = 1:length(phi)
    if yF(i) > Maximum
       Maximum = yF(i); % Update max value if current yF value is larger
       Maximum_x = phi(i); % Store corresponding x value
    end
    
   if yF(i) < Minimum
        Minimum = yF(i); % Update min value if current yF value is smaller
        Minimum_x = phi(i); % Store corresponding x value
   end
end

%===================== Velocity Calculation =====================

% Initialize array to store velocity values
vF = zeros(size(phi));

% Calculate velocity using central difference method
for i = 2:length(phi)-1
    vF(i) = (yF(i+1) - yF(i-1)) / (phi(i+1) - phi(i-1));
end

% Handle boundary points (endpoints)
vF(1) = (yF(2) - yF(1)) / (phi(2) - phi(1));
vF(end) = (yF(end) - yF(end-1)) / (phi(end) - phi(end-1));

% Multiply by angular velocity to get the actual velocity
vF = w * vF;

%===================== Acceleration Calculation =====================

% Initialize array to store acceleration values
aF = zeros(size(phi));
% Calculate acceleration using central difference method, excluding endpoints
for i = 3:length(phi)-2
    aF(i) = (vF(i+1) - vF(i-1)) / (phi(i+1) - phi(i-1));
end

%calculate the two endpoints for this function 
aF(2) = (vF(3) - vF(2)) / (phi(3) - phi(2));
aF(end-1) = (vF(end-1) - vF(end-2)) / (phi(end-1) - phi(end-2));

%exclude startpoint and last endpoint as there isn't enough information to
%calculate them
aF([1, end]) = NaN;

aF = aF * w;

%-------------------velocity 0 points-------------------
% Initialize array to store indices of zero points
velocity_0 = [];

% Detect sign changes in vF array
for i = 1:length(phi)-1
    if (vF(i) < 0 && vF(i+1) > 0) || (vF(i) > 0 && vF(i+1) < 0)
        velocity_0(end+1) = i; % Store index of zero point
    end
end

% Convert zeros to corresponding phi values
velocity_0 = phi(velocity_0);

%======================= OUTPUT ==================================

%plot the graph with all 3 functions
figure;
hold on
plot(phi, yF, 'black-', 'LineWidth', 0.8);
plot(phi, vF, 'b-', 'LineWidth', 0.8);
plot(phi, aF, 'r', 'LineWidth', 0.8);
plot(velocity_0, zeros(size(velocity_0)), 'ko', 'MarkerSize', 8);
hold off

axis ([0 6.5 -1 1]);
title ('Displacement, velocity and acceleration of slider F')
xlabel('Crank angle, rad');
ylabel('Displacement [m], velocity [m/s] and acceleration [m/s^2]');

grid on;

% Highlight points on the graph
legend('Displacement', 'Velocity', 'Acceleration', 'Location', 'SouthWest');

% Format x-axis ticks to show values in radians
xticks(0:pi/4:2*pi)
xticklabels({'0', '\pi/4','\pi/2', '3\pi/4', '\pi', '5\pi/4','3\pi/2', '7\pi/4','2\pi'})

% Output all calculated values

disp("The maximum height of Slider F (in m):");
disp(Maximum);
disp("The corresponding crank angle (in rad):");
disp(Maximum_x);
disp("The minimum height of slider F (in m):");
disp(Minimum);
disp("The corresponding crank angle (in rad):");
disp(Minimum_x);

disp("The points of zero velocity (in rad):")
disp(velocity_0)


