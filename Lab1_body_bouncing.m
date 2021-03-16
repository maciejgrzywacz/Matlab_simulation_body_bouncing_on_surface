close all;
clc;

syms z(x,y);
syms x_t(t);
syms y_t(t);
syms z_t(t);    
syms t;

%inputs
p = [0.5, 0.5, 1];        % starting positions
v = [-3, 0, 0.5];         % starting velocities
f(x,y) = x.^2+y.^2-x.*y;  % plane function
k = 1;                    % energy losses
%------

%constants
iterations = 5;
g = 10;
a = [0,0,-g];
m = 1;

syms f_dx_substituted(t);
f_dx = diff(f,x);

syms f_dy_substituted(t);
f_dy = diff(f,y);

% program outputs
times = zeros(iterations, 1);
positions = zeros(iterations, 3);
E_kins = zeros(iterations, 1);
E_pots = zeros(iterations, 1);
%-------

% draw plane
[X,Y] = meshgrid(-5:.5:5);
mesh(X, Y, double(f(X,Y)))
hold on
view([3, 0, 3]);

% iterate incidents
for i = 1:iterations
  
  % find time to incident
  % based on equation 5 in instruction
  left=[-1/2.*g, v(3), p(3)];
  z_substituted = f((p(1)+v(1).*t), (p(2)+v(2).*t));
  
  full_equation_coeffs = left - coeffs(z_substituted, 'All');
  ti = roots(full_equation_coeffs);    % find equation roots
  ti = ti(imag(ti)==0);         % get only real roots
  ti = max(ti);                % get greater value
  times(i, 1) = ti;
  
  %draw trajectory
  plot_dt = ti / 20;
  x_t(t) = p(1) + v(1) * t;
  y_t(t) = p(2) + v(2) * t;
  z_t(t) = p(3) + v(3) * t - 1/2*g*t.^2;  
  for t_w = 0:plot_dt:ti
      plot3(double(x_t(t_w)), double(y_t(t_w)), double(z_t(t_w)), '-or');
      drawnow
  end
  
  % new starting positions
  p = p + v .* ti + 1/2*a*ti.^2;
  positions(i,:) = p;  
  
  % new starting velocities and energies
  E_pots(i,:) = m .* g .* p(3);
  
  v(3) = v(3)-g*ti;           % velocity at point of collision before incident
  N = [-double(f_dx(p(1),p(2))),-double(f_dy(p(1),p(2))), 1];
  n = N / norm(N);
  E_kins(i,:) = m .* norm(v)^2/2;    
  v = sqrt(k) * (v - 2*(dot(v,n)).*n);    % velocity after incident
end

results_table_header = ["dt", "x", "y", "z", "E kin", "E pot", "E total"]; 
results = cat(2,times,positions, E_kins, E_pots, E_kins+E_pots);
cat(1, results_table_header, results)

hold off;
