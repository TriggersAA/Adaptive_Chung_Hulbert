clc;
clear;
%%----------------------------ASSIGNMENT 1----------------------%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The assingment has been broken down into 6 functions. 
%-------------------------------------------
% 1) just advances chung - hulbert in time. 
% 2) Just does the ZX error.
% 3) Implements the constant time-step 
% 4) Further implements the adptive method 
% 5) Implements the adaptive CH method 
% 6) implements plot function
%----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PARAMETERS

%Section and material parameters
EA = 500; L = 1;
rho_A = 1; g = 9.8;

% define M, Mass Matrix and K, stiffness matrix

M =[rho_A*(L^3)/3 0; 0 rho_A*L/3];
K = [rho_A*g*(L^2)/2 0; 0 EA/L];

%Second Task
%%%%%%%%%%%%
%    TIME STEP FOR BOTH ADATIVE AND NON-ADApTIVE CASE. THIS FUNCTION JUST 
%    ADVANCES IN A TIME STEP USING THE GENERALIZED-ALPHA SCHEME USING CHUNG
%    AND HULBERT
%%%%%%%%%%%%

% A function to calculate response for next time step

% Inputs are;
% u -> displacement vector, velocity -> velocity vector acc -> acceleration,
% a_f, a_m,beta, gamma -> C-Hulbert constants
% non_adaptive_inv -> inverse of k in the non-adative case to speed
% computation, sol_type -> string to distinguish both cases.

function [u_next, velocity_next, acc_next] = GA_Step(u_n, velocity_n, acc_n, r_na_np1, ...
    C,M,K,dt, a_f, a_m, beta, gamma,non_adaptive_inv, sol_type)

    % If adaptive method, compute effective stiffness inverse
    if strcmp(sol_type, "adaptive")
        k_eff = M * ((1 - a_m) / (beta * dt^2)) + C * gamma * (1 - a_f) / (beta * dt) + K * (1 - a_f);

        
        
    end

    % residual computation step
    r_eff = r_na_np1 - K * a_f * u_n ...
        + M * ( ((1 - a_m)/(beta * dt^2)) * u_n + ((1 - a_m)/(beta * dt))  * velocity_n ...
        + ((1 - a_m - 2 * beta)/(2 * beta)) * acc_n ) + ...
        C * ( (gamma * (1 - a_f)/(beta * dt)) * u_n + ((gamma - gamma * a_f - beta)/beta) * velocity_n ...
        + ((gamma - 2 * beta) * (1 - a_f) * dt / (2 * beta)) * acc_n );
    %uses solver type to compute displacement
    if strcmp(sol_type, "adaptive")
        u_next = k_eff\ r_eff;
        
    else

        u_next = non_adaptive_inv * r_eff;
    end
   %computes nre velocity
    velocity_next = (gamma / (beta * dt)) * (u_next - u_n) - ((gamma - beta) / beta) * velocity_n ...
             - ((gamma - 2 * beta) / (2 * beta)) * dt * acc_n;

    %computes new acceleration 
    acc_next = (1 / (beta * dt^2)) * (u_next - u_n) - (1 / (beta * dt)) * velocity_n ...
              - ((1 - 2 * beta) / (2 * beta)) * acc_n;

    
    %u_next     = (1 - a_f) * u_next     + a_f * u_n;
    %uDot_next  = (1 - a_f) * uDot_next  + a_f * uDot_n;
    %uDDot_next = (1 - a_m) * uDDot_next + a_m * uDDot_n;
end

% Task 4
%%----------------------------ZIENKIEWICZ- XIE ERROR FUNCTION----------------------%%
function [ZX_error, error_norm, rel_error] = ZX_Error(beta,u_n, u_next,acc_n, acc_next,dt)


    err_factor = (6*beta -1)/6;    % Error factor
    ZX_error = err_factor *(acc_next -acc_n) *dt^2;   % computes error vector
    error_norm = norm(ZX_error);         %computes error norm
    rel_error = error_norm/(norm(u_next - u_n)); % computes elative error
    
   

end

% Chung-Hulbert solver with a constant time march

function [u, velocity, acc,timings,ZX_container,times] = Const_Chung_Hulbert_Solver( rho_Inf,alpha_1, alpha_2, initial_dt,final_dt, vel_vec0, disp_vec0,M,K)
    
    % Time definitions
    t_0 = 0;
    t_f = final_dt;
    
    dt = initial_dt;
    N = ((t_f - t_0)/dt);   % Step number, initial

    timings = [t_0:dt:t_f];
    times = zeros(1,N); % stores time at each stage
    % Parameter definitions

    C = alpha_1*M + alpha_2*K;   % damping matrix C calculation
    % Generalized-Alpha parameters for Chung_Hulbert

    a_f = rho_Inf/(rho_Inf+1);
    a_m = (2*rho_Inf -1)/(rho_Inf+1);
    beta =0.25*(1-a_m+a_f)^2;
    gamma = (0.5 - a_m + a_f);
    

    %Data containers -> displacement,vel, acc with initial space allocation
    % no change due to constant time

    u = zeros(2,N+1);
    velocity = zeros(2,N+1);
    acc = zeros(2,N+1);
    
    %inital conditions -> u_0, uDot_0, uDDot_0
    u_0 = disp_vec0;
    uDot_0 =vel_vec0;
    uDDot_0 = M\(-C*uDot_0 - K*u_0);

    % Inverse of effective stiffness matrix calculation to time gains

    k_eff = M*((1-a_m)/(beta*dt^2)) + C*gamma*(1-a_f)/(beta*dt) + K*(1-a_f);

    % Initial conditions allocation to first entry in response parameters
    u(:,1) = u_0;
    velocity(:,1) = uDot_0;

    acc(:,1) =uDDot_0;
    
    K_inv = inv(k_eff);
    r_n_a = zeros(2,N+1); % input data for force. can be initailized to a value for forced conditions

    % ZX container; First two rows store the error vector
    % Next row stores the error _norm
    %forth row stores the relative error norm
    % fifth row stores the cummulative error
    ZX_container = zeros(5,N);
    % Loop,
    
    for n=1:N   %march foward in time
        [u(:,n+1), velocity(:,n+1), acc(:,n+1)] = GA_Step(u(:,n), velocity(:,n),...
            acc(:,n), r_n_a(:,n+1), C,M,K, dt, a_f, a_m, beta, gamma, K_inv,"");
        times(n) = dt;
        %compute ZX_error and stores it
        [ZX_container(1:2,n),ZX_container(3,n),ZX_container(4,n)]= ZX_Error(beta,u(:,n),u(:,n+1),acc(:,n),acc(:,n+1),dt);
        ZX_container(5,n) = ZX_container(3,n) + ZX_container(5,max(1,n-1));
  
    end
end

%Task 5
%%----------------------------FUNCTION FOR ADAPTIVE TIME UPDATE----------------------%%
% dt_old -> old data for time, eta_acceptable -> value given for error
% calculation, v1 and v2 for error bounds
%returns new time
function dt_new = Adaptive_time_update(dt_old,eta_acceptable, rel_error_new,v1,v2 )

    if (v1*eta_acceptable <= rel_error_new) && (rel_error_new <= v2*eta_acceptable) % updates based on error bounds
        dt_new = dt_old;
    else
        dt_new = dt_old *(sqrt(eta_acceptable/rel_error_new));
    end

end

%%----------------------------ADAPTIVE CH SCHEME----------------------%%

function [u, velocity, acc,cummulative_time,ZX_container,times] = Adaptive_Chung_Hulbert_Solver(rho_Inf,alpha_1, alpha_2, initial_dt,final_t, vel_vec0, disp_vec0,M,K)
    % Time definitions
    t_0 = 0;
    t_f = final_t;

    % allocates initial conditions to variables
    dt = initial_dt;
    v1 = 1; v2 = 10;
    eta_accepted = 1e-3;
    
    
    dt_0 = dt; %for time run

    % intial pre-allocation of data containers
    % % Step number, initial
    N = floor((t_f - t_0)/dt_0); % ensures a whole number is gotten

    times = zeros(1,N+1); % stores time at each stage
    cummulative_time = zeros(1,N+1); %sums time up

    % Parameter definitions
    C = alpha_1*M + alpha_2*K;
    a_f = rho_Inf/(rho_Inf+1);
    a_m = (2*rho_Inf -1)/(rho_Inf+1);
    beta = 0.25*(1-a_m+a_f)^2;
    gamma = 0.5 - a_m + a_f;
    

    %Data containers -> displacement,vel, acc

    u = zeros(2,N+1);
    velocity = zeros(2,N+1);
    acc = zeros(2,N+1);
     
    % ZX container; First two rows store the error vector
    % Next row stores the error _norm
    %forth row stores the relative error norm
    % fifth row stores the cummulative error

    ZX_container = zeros(5,N);
    
    %inital conditions -> u_0, uDot_0, uDDot_0
    u_0 = disp_vec0;
    uDot_0 =vel_vec0;

    uDDot_0 = M\(-C*uDot_0 - K*u_0);
    
    %k_eff = M*((1-a_m)/(beta*dt^2)) + C*gamma*(1-a_f)/(beta*dt) + K*(1-a_f);
    %Assigns initial onditions to response parameters
    u(:,1) = u_0;
    velocity(:,1) = uDot_0;
    acc(:,1) =uDDot_0;
    
    %K_inv = inv(k_eff);
    r_n_a = zeros(2,N+1);


    n = 1; % tracks number of loops
    while (cummulative_time(n) <= t_f)  %stops by time final
        times(n) = dt;
        
        
        % uses the adaptive case of GA_step. It takes static K_inf as
        % zero as it is not needed and stores all needed response data when called. 
        
        [u(:,n+1), velocity(:,n+1), acc(:,n+1)] = GA_Step(u(:,n), velocity(:,n),...
            acc(:,n), r_n_a(:,n+1), C,M,K, dt, a_f, a_m, beta, gamma, 0,"adaptive");
        
        %Executes the error calculations
        [ZX_container(1:2,n), ZX_container(3,n), ZX_container(4,n)]= ZX_Error(beta,u(:,n),u(:,n+1),acc(:,n),acc(:,n+1),dt);
        ZX_container(5,n) = ZX_container(3,n) +ZX_container(5,max(1,n-1)); % Accumulates error
        
        cummulative_time(n+1) = cummulative_time(n) + dt;  % Accumulates time
        %checks error bounds
        if (ZX_container(4,n) <= v1*eta_accepted) || (ZX_container(4,n) >= v2*eta_accepted)
            
            dt = Adaptive_time_update(dt,eta_accepted,ZX_container(4,n),v1,v2 );
            
                          
        end
        
        % Below logic block just resizes the container if needed dynmically at
        % runtime by multiplying data storage by 2 and putting all data in
        % again.
        if n+1 >= N
            new_N_max = 2 * N;
            
            u(:, N+1:new_N_max+1)     = 0;
            velocity(:, N+1:new_N_max+1)  = 0;
            acc(:, N+1:new_N_max+1) = 0;
            r_n_a(:, N+1:new_N_max+1) = 0;
            cummulative_time(1, N+1:new_N_max+1) = 0;
            
            N = new_N_max;
        end
        n = n + 1;  
   
    end
    
    % Below logic block just trims the excessive zeros after done with with
    % the ounter loop
    u = u(:,1:n);
    velocity = velocity(:,1:n);
    acc = acc(:,1:n);
    cummulative_time = cummulative_time(:,1:n);
    ZX_container = ZX_container(:,1:n-1);

    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------
%           PLOTS
%----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fig1, fig2] = CH_Plots(u,uDot, uDDot,time_cummulative,time, ZX_results,solver_type,params)
    
    % plot data

    alpha1 =params(2); alpha2=params(3); rho_inf =params(1);
    %Determine plot case
    if strcmp(solver_type , "adaptive")
        adaptive_status = 'ON';
    else
        adaptive_status = 'OFF';
    end
    % print a title
    title_str = sprintf(['Dynamic response of the system after 5 seconds\n' ...
                         'α₁ = %.2f, α₂ = %.2f, ρ∞ = %.2f | Adaptive time step: %s'], ...
                         alpha1, alpha2, rho_inf, adaptive_status);

    %%% ----------- FIGURE 1: HISTORY PLOTS ------------------------%%
    % Plots displacement vel and acceleration in that order

    fig1 = figure('Color','w','Name','History Plots');
    sgtitle(title_str, 'FontSize', 9);      % Figure title
    % Plot 1: Rotation and Linear Displacement plot, both together
    subplot(2,3,1);
    
    plot(time_cummulative,u(1,:), "LineWidth",1.3);
    hold on 
    plot(time_cummulative,u(2,:),LineWidth=1);
    hold off
    ylabel('θ (rad) and u (m)');
    xlabel('Time (s)');
    legend('$\theta$','u','Interpreter', 'latex')
    title('Rotation and Linear displacement plot');
    grid on
    subplot(2,3,2);


    plot(time_cummulative,uDot(1,:), "LineWidth",1);
    ylabel('Angular velocity rad/s ');
    xlabel('Time (s)');
    legend('$\dot{\theta}$','Interpreter', 'latex')
    title('Angular velocity plot'); grid on

    subplot(2,3,3);

    plot(time_cummulative,uDot(2,:), "LineWidth",1);
    ylabel('Axial Velocity (m/s)');
    xlabel('Time (s)');
    legend('$\dot{u}$','Interpreter', 'latex')
    title('Linear velocity plot');grid on

    % Plot 2: Angular Velocity separately
    subplot(2,3,4);

    plot(time_cummulative,uDDot(1,:), "LineWidth",1);
    ylabel('Angular acceleration\ $\mathrm{(rad/s^2)}$', 'Interpreter', 'latex');
    xlabel('Time (s)');
    legend('$\ddot{\theta}$','Interpreter', 'latex')
    title('Angular acceleration plot');grid on
    % Plot 3: Axial Velocity separately
    subplot(2,3,5);
    % Plot 4: Angular Acceleration separately
    plot(time_cummulative,uDDot(2,:), "LineWidth",1);
    ylabel('Acceleration\ $\mathrm{(m/s^2)}$','Interpreter', 'latex');
    xlabel('Time (s)');
    legend('$\ddot{u}$', 'Interpreter', 'latex')
    title('Linear acceleration plot');grid on

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%% ----------- FIGURE 2: ERROR & TIME STEP EVOLUTION --------%%

    % Plots time evolution, rel error norm, cuumulative error norm 
    fig2 = figure('Color','w','Name','Evolution Plots');
    sgtitle(title_str, 'FontSize', 9);      

    subplot(2,2,1);
    % Plot 6: Time Step Evolution
    plot(time, "LineWidth",1);
    ylabel('$\Delta t$ (s)','Interpreter', 'latex');
    xlabel('Steps');
    legend('$\Delta t$','Interpreter', 'latex')
    title('Time Evolution');grid on

    subplot(2,2,2);
    % Plot 7: Relative Error Evolution with bounds for adaptive case and
    % none for linear case
    % Parameters
           
    eta = 0.001;        
    v1 = 0;         
    v2 = 10;  
    % Bounding lines
    lower_bound = v1 * eta;
    upper_bound = v2 * eta;
    
    if solver_type ~= "adaptive"
        % For fixed time step, just show relative error norm
        plot(time_cummulative(:,2:end),ZX_results(4,:), "LineWidth",1);
        ylabel('Relative Error Norm');
        xlabel('Time March (s)');
        legend('$\eta$','Interpreter', 'latex');grid on
    else
        
        % Draw bound lines and relative error for adaptive case
        yline(lower_bound, 'DisplayName', 'v_1 \cdot \eta');
        hold on;
        yline(upper_bound, 'r--', 'LineWidth', 1, 'DisplayName', 'v_2 \cdot \eta');

        % Plot relative error
        plot(time_cummulative(:,2:end),ZX_results(4,:), 'LineWidth', 1, 'Color', [1 0.6 0], 'DisplayName', '$\eta_{rel}$');
        hold off;
        legend('$\eta$','Interpreter', 'latex'); ylabel('Relative Error Norm');grid on
    end
    title('Relative Error Evolution');

    subplot(2,2,3);
    % Plot 8: Cumulative Error Norm plots
    plot(time_cummulative(:,2:end),ZX_results(5,:), "LineWidth",1);
    ylabel('Cummulative Error Norm');
    xlabel('Time (s)');
    legend('$\sum_{i} \| e \|$','Interpreter', 'latex')
    title('Cummulative Error Norm Evolution');grid on

end


%%% DATA LOADING
u_0 = [0;-L/5]; uDot_0 = [sqrt(g/(6*L));0];

%Const_Chung_Hulbert_Solver( rho_Inf,alpha_1, alpha_2, initial_dt,final_dt, vel_vec0, disp_vec0,M,K)
% Adaptive_Chung_Hulbert_Solver(rho_Inf,alpha_1, alpha_2, initial_dt,final_t, vel_vec0, disp_vec0,M,K)

param1 =[1,1,0]; param2=[0.1,0,0];  %rho infinity, alpha1, aplha2
[u1, uDot1, uDDot1,timings1,ZX_container1,time1] = Const_Chung_Hulbert_Solver(1,1,0,0.01,5,uDot_0,u_0,M,K);
[u2, uDot2, uDDot2,timings2,ZX_container2,time2] = Const_Chung_Hulbert_Solver(0.1,0,0,0.01,5,uDot_0,u_0,M,K);

fprintf('Value of theta and u after 5 seconds are %.5f and %.5f for 3(a)\n', u1(1,end), u2(2,end));
fprintf('Value of theta and u after 5 seconds are %.5f and %.5f for 3(b)\n', u2(1,end), u2(2,end));

[u3, uDot3, uDDot3,timings3,ZX_container3,time3] = Adaptive_Chung_Hulbert_Solver(1,1,0,0.01,5,uDot_0,u_0,M,K);
CH_Plots(u1,uDot1,uDDot1,timings1,time1,ZX_container1,"",param1);
CH_Plots(u2,uDot2,uDDot2,timings2,time2,ZX_container2,"",param2);
CH_Plots(u3,uDot3,uDDot3,timings3,time3,ZX_container3,"adaptive",param1);

%%% DATA display
%disp([u1(:,end), u2(:,end)])

fprintf('Constant-step case: cumulative error = %.5e, steps = %d\n', ZX_container1(5,end), size(ZX_container1,2));
fprintf('Adaptive-step case: cumulative error = %.5e, steps = %d\n', ZX_container3(5,end), size(ZX_container3,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%code for comparison redacted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

