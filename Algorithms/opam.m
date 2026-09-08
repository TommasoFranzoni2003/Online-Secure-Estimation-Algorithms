clc; clear all; close all;

%% SETTINGS

% loading of the dataset => A, C, x0, a, q, n, h 
load('dynamic_CPS_data.mat')

% defining the parameters linked to the system dynamically
n = size(A, 1); 
q = size(C, 1); 
h = 3;
k_limit = 500;

% initialization of arrays for the measures to be evaluated
state_est_error_x= zeros(k_limit,1); % k_limit = number of runs that we want to do
support_attack_error= zeros(k_limit,1);
false_positive= zeros(k_limit,1);
false_negative= zeros(k_limit,1);

% computation of the true-values for the support 
true_support= (a ~= 0);


%% O-PAM ALGORITHM

mu = 1;
lambda_pam = 0.5;

% initialization of x(0)
x_new= randn(n,1);
a_new= randn(q,1);
for k = 1:k_limit 
    
    % saving of the "previous" values to make the estimate
    x_old = x_new;
    a_old = a_new;
    
    % computing "true" y for the following iteration
    y_k= C*(A^k) * x0 + a;
    
    % step 1 -> fix a (old) and solve for x
    first_term_x_new= ((C*(A^k))'*(C*(A^k)) + mu*eye(n))^(-1);
    second_term_x_new= (mu*x_old + (C*(A^k))'*(y_k - a_old));
    
    x_new= first_term_x_new * second_term_x_new;
    
    % step 2 -> fix x (new) and solve for a
    S_app_val= y_k - (C*(A^k))*x_new + mu*a_old;
    a_new= (1/(1+mu)) * sign(S_app_val).*max(0, abs(S_app_val) - lambda_pam);
     
    % calculation of state estimation error for the k-th iteration 
    state_est_error_x(k) = ((norm((A^k)*x_new - (A^k)*x0)^2))/n;
    
    % calculation of support attack error for the k-th iteration 
    est_supp= (a_new~= 0);
    support_attack_error(k)= sum(abs(est_supp - true_support));
    
    % calculation of false positives for the k-th iteration 
    false_positive(k)= sum(max(0, (est_supp - true_support)));
    
    % calculation of false negatives for the k-th iteration 
    false_negative(k)= -sum(min(0, (est_supp - true_support)));
end


%% PLOTTING OF THE RESULTS

figure(1);

% state estimation error
subplot(2, 2, 1);
% use semilogy for log-scale
semilogy(1:k_limit, state_est_error_x, "g", 'LineWidth', 1.5);
title("State Estimation Error O-PAM");
xlabel("Time step (k)");
ylabel("Error (Log Scale)");
grid on;

% support attack error
subplot(2, 2, 2);
plot(1:k_limit, support_attack_error, "g", 'LineWidth', 1.5);
title("Support Attack Error O-PAM");
xlabel("Time step (k)");
grid on;

% false positives
subplot(2, 2, 3);
plot(1:k_limit, false_positive, "g", 'LineWidth', 1.5);
title("False Positive O-PAM");
xlabel("Time step (k)");
grid on;

% false negatives
subplot(2, 2, 4);
plot(1:k_limit, false_negative, "g", 'LineWidth', 1.5);
title("False Negative O-PAM");
xlabel("Time step (k)");
grid on;