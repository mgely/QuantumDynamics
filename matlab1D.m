% The spatial grid

N = 1024; % number of interior grid points
L = 100; % system extends from x=0 to x=L
h = L / N; % grid spacing
tau = 0.1; % time step
simulation_time = 100; % simulation time
x = h:h:L; % coordinates of grid points






% The potential V(x)

V_0 = 0.5; % height of potential barrier
V_width = 5; % width of potential barrier
V_center = L/2; % center of potential barrier
half_width = abs(0.5 * V_width);
V = zeros(N,1);

%%%%%%%%%%%% CHOOSE BETWEEN TUNELING AND HARMONIC POTENTIAL %%%%%%%%%%%%%%%%%%
%Tunelling
% for j = 1:N
%     if (abs(j*h - V_center) <= half_width)
%         V(j)  = V_0;
%     end
% end

%Harmonic
%V = 0.1*(x-V_center).^2;






% Initial wave packet
x_0 = L/2-5; % location of center
E = 1; % average energy
sigma_0 = L / 100; % width of wave packet
psi_norm = 1 / sqrt(sigma_0 * sqrt(pi)); % norm of psi
k_0 = 0; % average wavenumber
psi_exp_factor = exp(-(x-x_0).^2./ (2 * sigma_0 * sigma_0));

psi  = zeros(N,1);
for j = 1:N
    psi(j) = complex(cos(k_0*x(j)),sin(k_0*x(j)))*psi_exp_factor(j); % complex wavefunction
end





% initialize the phase rotation factors
T_exp_factor = zeros(N,1);
V_exp_factor = zeros(N,1);
for j = 1:N
% kinetic factor exp[-iT/h_bar tau]
    if j < N / 2
        p = j;
    else
        p = j - N;
    end
    
    p = p * 2 * pi / L;
    theta = - p * p / 2 * tau;
    T_exp_factor(j) = complex(cos(theta), sin(theta));

% potential factor exp[-iV(x)/(2h_bar) tau]
    theta = - V(j) / 2 * tau;
    V_exp_factor(j) = complex(cos(theta), sin(theta));
end




% Simulation
t = 0; % time
j = 1;%frame index


while t<simulation_time
    % first half potential phase rotation
    for i = 1:N
            psi(i) = psi(i) * V_exp_factor(i);
    end
    
    
    % FFT to momentum space
    psi = fft(psi);
%     psi = fftshift(psi);
    
    % kinetic phase rotation
    for i = 1:N
            psi(i) = psi(i) * T_exp_factor(i);
    end
    
    % FFT back to position space
%     psi = ifftshift(psi);
    psi = ifft(psi);
    
    % second half potential phase rotation
    for i = 1:N
            psi(i) = psi(i) * V_exp_factor(i);
    end
    
    % Store the frame
    plot(x,abs(psi),x,real(V))
    axis([0 L -1.5 1.5])
    M(j)=getframe(gcf); % leaving gcf out crops the frame in the movie.
    
    j = j +1;
    t = t + tau;
end


