%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jackson L. Cole
% M01250797
% Middle Tennessee State University
% PHYS 3200 :: Spring 2017 :: Dr. Robertson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scientific Modeling Project
% Topic: Basic n-body simulation using the Euler-Cromer method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information about the program for the user
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Gravitational n-Body Simulation using the Euler-Cromer Method')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Author		: Jackson L. Cole')
disp('Affiliation	: Middle Tennessee State University')
disp('Semester		: Spring 2017')
disp('Course		: PHYS 3200: Scientific Modeling')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Details:')
disp('This program makes use of the Euler-Cromer method to simulate gravitational interactions among a specified number of bodies. In other words, this is just an n-Body simulation. The program does not seek to simulate any particular real system, but a starting point for the range of masses of the bodies included was the recent LIGO discovery of merging black holes. The range of possible masses of the bodies in the simulation is from less than 1 solar mass to upwards of 50 solar masses.')
disp('###')
disp('The user will be prompted to input a number of bodies for the simulation. The choice of the number of bodies is left up to the user, of course, but one may be advised to choose LOW NUMBERS...')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Asking for user input
numberofbodies = input('Please enter a number of bodies to be included in the simulation: ');
time = input('Please enter a duration of observation (Be advised that 50 [years] seems to be reasonable to see something neat without being too excessive...): ');
disp(['The simulation will simulate ' num2str(numberofbodies) ' bodies over a period of ' num2str(time) ' years.'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing the time and time step
lowertbound = 0;		% [years]
uppertbound = time;		% [years]
stepsize=0.0001;		% [years]

t=lowertbound:stepsize:uppertbound;	% Time values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants and basics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Units of AU, Solar Masses, and years seemed to play much more nicely with MATLAB
G = (6.67408e-11)*(1.989e30)*(3.154e+7)^2/(1.496e+11)^3; % [AU^3 msol^-1 yr^-2]
n = numberofbodies;				% Number of bodies in the simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocating arrays for position vectors and velocity vectors
x = zeros(n,numel(t));		% x-position
y = zeros(n,numel(t));		% y-position
z = zeros(n,numel(t));		% z-position
vx = zeros(n,numel(t));		% x velocity
vy = zeros(n,numel(t));		% y velocity
vz = zeros(n,numel(t));		% z velocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of n masses
masses = zeros(1,n);				% Preallocating an array for masses
lowmass = 1;	% [solar masses]	% Lowest possible mass

% The following loop generates random masses
for i=1:n
    masses(i)=(lowmass)*rand*50;    % Generates random masses in units of msol
end
masses(1) = lowmass * 35;			% Quite massive mass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of mass products (explanation below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following generates every possible combination of mass products
% for use in the equation for the force due to gravity
% (Fg = Gm(i)m(j)/r^3). The resulting array is multidimensional, and
% several further operations are performed prior to its use.

massproducts = [];			% Initializing the array without preallocating
							%%a maximum size

% The following loop iterates through each possible combination of
% masses to create an array where each possible combination of
% products is represented. Naturally, this results in a scenario
% where the upper diagonal and the lower diagonal contain essentially
% the same data. This is important for later.
for i=1:n
    for k=1:n
        massproducts(i,k)=masses(i) * masses(k);
    end
end

% At this point, it is wise to acknowledge that in the above loop,
% we end up with products that are essentially masses(i) * masses
% (i), which would correspond in the primary loop of the program
% to an introduction of a force between mass(i) and mass(i), which
% is a force due to a body applied on THAT body. This is of course
% problematic, so the below operation removes the innappropriate
% products from the massproducts array (which occur on the diagonal),
% and replaces all of these masses with 0.

massproducts(1:n+1:n*n) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determining initial conditions (which are largely randomized)

% The simulation is setup to take place in a space with a radius of
% roughly 30 AU (this has the most fun results for such large masses.

lim = 100;						% Limiting radius
vmax = 100;						% Maximum allowed velocity initial velocity
								%%in any direction
v = round(sqrt(((vmax)^2)/3));	% Finding a component of the velocity
								%%that would give us that velocity;
								%%%from Keplers third law

% The following loop sets up the initial conditions using random integers
for i=1:n
    x(i,1) = randi([-lim lim]);	% Random initial x position
	y(i,1) = randi([-lim lim]);	% Random initial y position
	z(i,1) = randi([-lim lim]);	% Random initial z position

	vx(i,1) = rand*randi([-v v]);	% Random initial x velocity
	vy(i,1) = rand*rand*randi([-v v]);	% Random intiial y velocity
	vz(i,1) = rand*rand*rand*randi([-v v]);	% Random initial z velocity

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% The folling conditional statements serve to reduce the initial velocity
	% to more of a glancing velocity with respect to the center of mass or
	% the system. This is accomplished by simply using separate intervals in
	% which only one of the velocities is allowed to remain as it was.
	temp = [x(i,1) y(i,1) z(i,1)];
	if max(temp) == temp(1)
		vx(i,1) == 0.001*vx(i,1);
	elseif max(temp) == temp(2)
		vy(i,1) == 0.001*vy(i,1);
	elseif max(temp) == temp(3)
		vz(i,1) == 0.001*vz(i,1);
	end
	rr=rand;						% Generates a random fraction
	if rr > 0.66 && rr < 0.99		% Reduces the velocity to the z component
		if rr > 0.66 && rr < 0.8249
			vx(i,1) = 0.001*vx(i,1);
		elseif rr > 0.8250 && rr < 0.99
			vy(i,1) = 0.001*vy(i,1);
		end
	elseif rr > 0.33 && rr < 0.65	% Reduces the velocity to the x component
		if rr > 0.33 && rr < 0.49
			vy(i,1) = 0.001*vy(i,1);
		elseif rr > 0.491 && rr < 0.65
			vz(i,1) = 0.001*vz(i,1);
		end
	elseif rr > 0.01 && rr < 0.32	% Reduces the velocity to the y component
		if rr > 0.01 && rr < 0.1649
			vx(i,1) = 0.001*vx(i,1);
		elseif rr > 0.1650 && rr < 0.32
			vz(i,1) = 0.001*vz(i,1);
		end
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% The following specific set of initial conditions ensures that the massive
% mass specifically specified above is located at the origin. This allows
% us to have somewhat of a "control" body to observe.
x(1,1) = 0;
y(1,1) = 0;
z(1,1) = 0;

vx(1,1) = rand;
vy(1,1) = rand;
vz(1,1) = rand;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Loop of Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop overview: This loop accomplishes several things which will be mentioned
% briefly here. First, the loop populates a multidimensional array of radii
% between bodies in the simulation. These radii are then used in the
% calculation of the force between bodies, which is used in the Euler-Cromer
% loop below. That description is somewhat vague, but the procedure I came up
% with will be explained in further detail below.
for k = 1:numel(t)
    clear r			% Making sure nothing is left from the last iteration

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% The following nested loop generates an array of radii between each body.
	% That is, there are n rows, where each row essentially represents a list
	% of radii between the nth body and each other body (notice that this
	% mirrors the mass product generation performed previously). Also, notice
	% that this generates a radius between the nth body and the nth body,
	% which returns an innappropriate value for that element. These values are
	% removed with an identical procedure to the one used to remove the
	% inappropriate mass products from the relevant diagonal.
    for i = 1:n
        for j = 1:n
            r(i,j) = sqrt( (x(j,k)-x(i,k)).^2 + (y(j,k)-y(i,k)).^2 + (z(j,k)-z(i,k)).^2 );
        end
    end
    r(1:n+1:n*n) = 0;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% The following nested loop is where the actual calculations happen, for
	% lack of a better introduction. Here, an element by element operation is
	% performed using the nth (or ith) row of the massproduct array and the nth
	% (or ith) row of the radii array to compute the mm/r^3 quanitity relevent
	% to each pair of masses in the simulation. Notice that at this point, the
	% innappropriate values removed above just return a 0 quantity.
	% Given that the distributive property of multiplication applies (because
	% it does), a summation is performed of each mm/r^3 "row" that is generated,
	% which is then multiplied as G*(sum(mm/r^3))*position*dt to determine the
	% new velocities.
    for i=1:n
        clear massproductoverradiuscubed	% Clearing this out
		clear summm							% Clearing this out
        massproductoverradiuscubed = (1/masses(i))*massproducts(i,:)./(r(i,:).^3);
        massproductoverradiuscubed(isnan(massproductoverradiuscubed)) = 0;
        summm=sum(massproductoverradiuscubed);
        vx(i,k+1) = vx(i,k) - ( ((G*x(i,k)))*stepsize)*summm;
        vy(i,k+1) = vy(i,k) - ( ((G*y(i,k)))*stepsize)*summm;
        vz(i,k+1) = vz(i,k) - ( ((G*z(i,k)))*stepsize)*summm;
    end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Standard Euler-Cromer procedure
    x(:,k+1) = x(:,k) + vx(:,k+1) * stepsize;
    y(:,k+1) = y(:,k) + vy(:,k+1) * stepsize;
    z(:,k+1) = z(:,k) + vz(:,k+1) * stepsize;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Operations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The individual rows (i.e. the position v. time data) for each body must
% be plotted using a for loop to work in all cases. I ran into several
% issues getting the loop to work with the multidimensional array, so I
% just converted each row into a cell array for plotting.
xcells = cell(1,numel(t));
ycells = cell(1,numel(t));
zcells = cell(1,numel(t));
for i=1:n
	xcells{i} = x(i,:);
	ycells{i} = y(i,:);
	zcells{i} = z(i,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determining the best limits on the 3d plot to maximize the data viewed
xabs = mean(abs(x(:)));				% Absolute values for entire x array
yabs = mean(abs(y(:)));				% Absolute values for entire y array
zabs = mean(abs(z(:)));				% Absolute values for entire z array

maxes = [0 0 0];			% Preallocating space for the maximum values
							%%in x,y,z
maxes(1) = mode(10*max(xabs));	% Returns a single value for the max(abs(x))
maxes(2) = mode(10*max(yabs));	% Returns a single value for the max(abs(y))
maxes(3) = mode(10*max(zabs));	% Returns a single value for the max(abs(z))
maxi = max(maxes);			% Finds the maximum absolute value in ALL
							%%position arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
figure
hold all			% Holds all plots to ensure everything is plotted
ax = gca;
ax.XGrid = 'on';	% XGrid
ax.YGrid = 'on';	% YGrid
ax.ZGrid = 'on';	% ZGrid
Az = 47;			% Viewing options: azimuth
El = 6;			% Viewing options: elevation
view(Az,El);		% Sets viewing options
scatter3(x(:,1),y(:,1),z(:,1))	% Scatter plot of starting points
								%%for the n-bodies
scatter3(x(:,n),y(:,n),z(:,n))	% Scatter plot of the ending points
								%%for the n-bodies

% The following loop plots the cartesian coordinate data for each body
for i=1:numel(xcells)
	plot3(xcells{i},ycells{i},zcells{i})
end
hold on
xlim([-maxi maxi])
ylim([-maxi maxi])
zlim([-maxi maxi])
title(['$n$-Body Simulation: $n=' num2str(numberofbodies) '$, $t=[0,' num2str(time) ']$ [years]'],'interpreter','latex')
xlabel('$x$ [AU]','interpreter','latex')
ylabel('$y$ [AU]','interpreter','latex')
zlabel('$z$ [AU]','interpreter','latex')
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notifying me when my program is done running.
load handel
sound(y,Fs)
