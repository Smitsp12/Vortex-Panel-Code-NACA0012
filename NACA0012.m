clc;
%VORTEX PANEL CODE NACA 0012

%Defining initial condition for the problem
Rho = 1.1225;           %freestream density
Vinf = 5.5;             %free stream velocity
alpha = 6*pi/180;       %ngle of attack in radians.

%creating fine grid at locations of high curvature (Leading and Trailing
%edge) and coarse grid for remaining geometry.

x1 = linspace(0,0.25,80);       %from leading edge to quater chord (c/4)
x2 = linspace(0.25,0.75,50);    %(c/4 to 3c/4)
x2(1) = [];                     %to avoid repeating last element of x1 which is same as first element of x2
x3 = linspace(0.75,1,80);       %(3c/4) to trailing edge
x3(1) = [];                     %to avoid repeating last element of x1 which is same as first element of x2
x = [x1';x2';x3'];              %conocatenating all x coordinates
x = flip(x);

%calculating y coordinates from equation for the lower surface of NACA 0012
y = -0.594689181*(0.298222773*sqrt(x) - 0.127125232*x - 0.357907906*x.^2 + 0.291984971*x.^3 - 0.105174606*x.^4); 

X = flip(x);                    %reversing  the order of elements to calculate y coordinates for upper surface.

X(1,:) = [];                    %removing the coordinate (0,0) as it is already considered in lower surface.

%calculating y coordinates from equation for the upper surface of NACA 0012
Y = 0.594689181*(0.298222773*sqrt(X) - 0.127125232*X - 0.357907906*X.^2 + 0.291984971*X.^3 - 0.105174606*X.^4);
y = [y; Y];                     %conocatenating x and y coordinates for upper and lower surfaces.
x = [x; X];

%panels are indexed starting from lower surface near trailing edge and
%moving in counter-clockwise direction
N = length(x) - 1;              %Defining number of panels 


s = zeros(N,1);                 %matrix [s] for length of panel
for i = 1:N                     
    s(i) =  sqrt((x(i+1) - x(i))^2 + (y(i+1) - y(i))^2);
end

phi = zeros(N,1);               %matrix [phi] for angle w.r.t +X-axis matrix
for i = 1:N/2                   
    phi(i) = pi + atan((y(i+1)-y(i))/(x(i+1)-x(i)));    %formula is modified to get angle value w.r.t +x-axis
end

for i = (N/2)+1:N
    phi(i) = atan((y(i+1)-y(i))/(x(i+1)-x(i)));
    if phi(i) < 0               %formula is modified to get angle value w.r.t +x-axis
        phi(i) = 2*pi + phi(i); 
    end
end

xm = zeros(N,1);                %matrix [xm, ym] for control point of each panel.
ym = zeros(N,1);
for i = 1:N
    xm(i) = (x(i+1)+x(i))/2;
    ym(i) = (y(i+1)+y(i))/2;
end

%evaluating the coefficients of the integral term in set of equations
A = zeros(N,N);
B = zeros(N,N);
C = zeros(N,N);
D = zeros(N,N);
E = zeros(N,N);
I = zeros(N,N);

for i = 1:N
    for j = 1:N
        A(i,j) = -(xm(i) - x(j))*cos(phi(j)) - (ym(i) - y(j))*sin(phi(j)); 
        B(i,j) = (xm(i) - x(j))^2 + (ym(i) - y(j))^2;
        C(i,j) = -cos(phi(i) - phi(j));
        D(i,j) = (xm(i) - x(j))*cos(phi(i)) + (ym(i) - y(j))*sin(phi(i));
        E(i,j) = sqrt(abs(B(i,j) - A(i,j)^2));
     end
end

for i = 1:N
    for j = 1:N
        I(i,j) = (C(i,j)/2)*log((s(j)^2 + 2*A(i,j)*s(j) + B(i,j))/B(i,j)) + ((D(i,j)-A(i,j)*C(i,j))/sqrt(B(i,j)-A(i,j)^2))*(atan((s(j) + A(i,j))/E(i,j)) - atan(A(i,j)/sqrt(B(i,j)-A(i,j)^2)));
        if i == j               %condition to ensure that diagonal element of I are always 0 and not NaN.
            I(i,j) = 0;
        end
    end
end

%Imposing Kutta condition by removing equation for the Nth panel and adding
%gamma(1) + gamma(N) = 0

I(N,:) = zeros(1,N);
I(N,1) = 1;
I(N,N) = 1;

v = 2*pi*Vinf*sin(phi-alpha);   %matrix [v] for RHS 2*pi*Vinf*sin(phi-alpha) terms. 
v(N,1) = 0;                     %Last element of [v] will be zero due to kutta condition as discussed above.

gamma = I\v;                    %solving gamma = inverse[I] * [v]
Gamma = s .* gamma;             %total circulation = summation of (gamma * panel length) for all panels.

cl = 2*sum(Gamma)/Vinf;         %sectional lift coefficient = 2*total circulation/Vinf. sum(Gamma) will add all elements of Gamma.
disp("sectional lift coefficient = "+cl)
disp("number of panel used = "+N)