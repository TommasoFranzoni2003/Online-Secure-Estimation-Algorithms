clear all; close all; clc;

%% SETTINGS

% loading of the dataset => A, C, x0, a, q, n, h 
load('dynamic_CPS_data.mat')

% defining the parameters linked to the system dynamically to prevent dimension mismatch
n = size(A, 1); 
q = size(C, 1); 
h = 3;
k_limit = 500;

% initialization of the arrays to collect the state estimation errors
state_est_error_x_opgd= zeros(k_limit,1); %k_limit = number of runs we want to do
state_est_error_x_opam= zeros(k_limit,1);
state_est_error_x_aopgd= zeros(k_limit,1);
state_est_error_x_aopam= zeros(k_limit,1);

% initialization of the arrays to collect the support attack errors
support_attack_error_opgd= zeros(k_limit,1);
support_attack_error_opam= zeros(k_limit,1);
support_attack_error_aopgd= zeros(k_limit,1);
support_attack_error_aopam= zeros(k_limit,1);

% initialization of the arrays to collect the false positives
false_positive_opgd= zeros(k_limit,1);
false_positive_opam= zeros(k_limit,1);
false_positive_aopgd= zeros(k_limit,1);
false_positive_aopam= zeros(k_limit,1);

% initialization of the arrays to collect the false negatives
false_negative_opam= zeros(k_limit,1);
false_negative_aopgd= zeros(k_limit,1);
false_negative_aopam= zeros(k_limit,1);

% computation of the "true" values for the support
true_support= (a ~= 0);


%% O-PGD ALGORITHM

nu = 0.01; % stepsize
lambda_pgd = 0.5;
% initialization of x(0)
x_new= randn(n,1);
a_new= randn(q,1);
for k = 1:k_limit 
    
    % computation of "true" y for the following iteration 
    y_k= C*(A^k) * x0 + a;
    
    % saving of the "previous" values to make the estimate
    x_old = x_new;
    a_old = a_new;
    
    % step 1
    % computation of the gradient of F respect to x
    gradient_x= (C*(A^k))'*((C*(A^k))*x_old + a_old - y_k);
    
    % calculation of the current estimate of x
    x_new = x_old - nu * gradient_x;
    
    % step 2
    % computation of the gradient of F respect to a
    gradient_a= (C*(A^k)) * x_old + a_old - y_k; 
    
    % calculation of S with nu and lamda
    S_app_val = a_old - nu*gradient_a;
    
    % computation of the current estimate of a
    a_new= sign(S_app_val).*max(0, abs(S_app_val) - nu*lambda_pgd);
    
    % computation of the state estimation error of the current iteration
    state_est_error_x_opgd(k) = ((norm((A^k)*x_new - (A^k)*x0)^2))/n;
    
    % computation of the support attack error of the current iteration
    est_supp_opgd= (a_new~= 0);
    support_attack_error_opgd(k)= sum(abs(est_supp_opgd - true_support));
    
    % computation of the false positives for the current iteration
    false_positive_opgd(k)= sum(max(0, (est_supp_opgd - true_support)));
    
    % computation of the false negatives for the current iteration
    false_negative_opgd(k)= -sum(min(0, (est_supp_opgd - true_support)));
end


%% O-PAM ALGORITHM

mu = 1;
lambda_pam = 0.5;
% initialization of x(0) and a = 0
x_new= randn(n,1);
a_new= randn(q,1);
for k = 1:k_limit 
   
    % saving of the previous values to make the estimate
    x_old = x_new;
    a_old = a_new;
    
    % computation of the "true" y for the following iteration
    y_k= C*(A^k) * x0 + a;
    
    % step 1 -> fix a
    first_term_x_new= ((C*(A^k))'*(C*(A^k)) + mu*eye(n))^(-1);
    second_term_x_new= (mu*x_old + (C*(A^k))'*(y_k - a_old));
    
    x_new= first_term_x_new * second_term_x_new;
    
    % step 2 -> fix x
    S_app_val= y_k - (C*(A^k))*x_new + mu*a_old;
    a_new= (1/(1+mu)) * sign(S_app_val).*max(0, abs(S_app_val) - lambda_pam);
     
    % computation of the state estimation error of the current iteration
    state_est_error_x_opam(k) = ((norm((A^k)*x_new - (A^k)*x0)^2))/n;
    
    % computation of the support attack error of the current iteration
    est_supp_opam= (a_new~= 0);
    support_attack_error_opam(k)= sum(abs(est_supp_opam - true_support));
    
    % computation of the false positives for the current iteration
    false_positive_opam(k)= sum(max(0, (est_supp_opam - true_support)));
    
    % computation of the false negatives for the current iteration
    false_negative_opam(k)= -sum(min(0, (est_supp_opam - true_support)));
end


%% AO-PGD ALGORITHM

nu = 0.01; % stepsize
lambda_aopgd = 0.5;
% initialization of x(0) and a = 0
x_new= randn(n,1);
a_new= randn(q,1);
for k = 1:k_limit 
    % saving of the previous values to make the estimate
    x_old = x_new;
    a_old = a_new;
    
    % step 1
    % computation of the gradient of F respect to x
    % in this case it is the sum of the gradients from j=1 to the current k
    
    gradient_x_sum= zeros(n,1);
    for j=1:k
        % computation of the "true" y for the following iteration
        y_j= C*(A^j) * x0 + a;
        
        % computation of the gradient
        gradient_x_sum= gradient_x_sum + (C*(A^j))'*((C*(A^j))*x_old + a_old - y_j);
    end
   
    % calculation of the current estimate of x
    x_new = x_old - nu * (1/(k+1)) * gradient_x_sum;
    
    % step 2
    % computation of the gradient of F respect to a
    % in this case it is the sum of the gradients from j=1 to the current k
    gradient_a_sum= zeros(q,1);
    
    for j=1:k
        % computation of the "true" y for the following iteration
        y_j= C*(A^j) * x0 + a;
        
        % computation of the gradient
        gradient_a_sum= gradient_a_sum + (C*(A^j)) * x_old + a_old - y_j; 
    end
    
    % calculation of S with nu and lamda
    S_app_val = a_old - nu*(1/(k+1))*gradient_a_sum;
    
    % computation of the current estimate of a
    a_new= sign(S_app_val).*max(0, abs(S_app_val) - nu*lambda_aopgd);
    
    % computation of the state estimation error of the current iteration
    state_est_error_x_aopgd(k) = ((norm((A^k)*x_new - (A^k)*x0)^2))/n;
    
    % computation of the support attack error of the current iteration
    est_supp_aopgd= (a_new~= 0);
    support_attack_error_aopgd(k)= sum(abs(est_supp_aopgd - true_support));
    
    % computation of the false positives for the current iteration
    false_positive_aopgd(k)= sum(max(0, (est_supp_aopgd - true_support)));
    
    % computation of the false negatives for the current iteration
    false_negative_aopgd(k)= -sum(min(0, (est_supp_aopgd - true_support)));
end


%% AO-PAM ALGORITHM

mu = 1;
lambda_aopam = 0.5;
% initialization of x(0) and a = 0
x_new= randn(n,1);
a_new= randn(q,1);
for k = 1:k_limit 
    % saving of the previous values to make the estimate
    x_old = x_new;
    a_old = a_new;
    
    % step 1 -> fix a
    % in this case I sum the terms from j=1 up to the current k
    sum_hessian_x= zeros(n,n);
    sum_grad_x= zeros(n,1);
    
    for j=1:k
        % computation of the "true" y for the following iteration
        y_j= C*(A^j) * x0 + a;
        
        % computation of the terms
        sum_hessian_x= sum_hessian_x + (C*(A^j))'*(C*(A^j));
        sum_grad_x= sum_grad_x + (C*(A^j))' * (y_j - a_old);
    end
    
    % computation of x_new by inverting the matrix at the end and using the mean (1/k)
    first_term_x_new= ((1/k)*sum_hessian_x + mu*eye(n))^(-1);
    second_term_x_new= (1/k)*sum_grad_x + mu*x_old;
    
    x_new= first_term_x_new*second_term_x_new;
    
    % step 2 -> fix x
    % in this case I sum the terms from j=1 up to the current k
    sum_err_a= zeros(q,1);
    
    for j=1:k
        % computation of the "true" y for the following iteration
        y_j= C*(A^j) * x0 + a;
        sum_err_a= sum_err_a + (y_j - C*(A^j)*x_new); 
    end
        
    % computation of the resulting value
    S_app_val= ((1/k)*sum_err_a + mu*a_old) / (1 + mu);
    a_new= sign(S_app_val) .* max(0, abs(S_app_val) - lambda_aopam/(1+mu));
     
    % computation of the state estimation error of the current iteration
    state_est_error_x_aopam(k) = ((norm((A^k)*x_new - (A^k)*x0)^2))/n;
    
    % computation of the support attack error of the current iteration
    est_supp_aopam= (a_new~= 0);
    support_attack_error_aopam(k)= sum(abs(est_supp_aopam - true_support));
    
    % computation of the false positives for the current iteration
    false_positive_aopam(k)= sum(max(0, (est_supp_aopam - true_support)));
    
    % computation of the false negatives for the current iteration
    false_negative_aopam(k)= -sum(min(0, (est_supp_aopam - true_support)));
end


%% PRINTING OF THE RESULTS OF THE ALGORITHMS COMPARED

figure(1);

% state estimation error
subplot(2, 2, 1);
semilogy(1:k_limit, state_est_error_x_opgd, "g", 'LineWidth', 1.5);
hold on;
semilogy(1:k_limit, state_est_error_x_opam, "r", 'LineWidth', 1.5);
semilogy(1:k_limit, state_est_error_x_aopgd, "b", 'LineWidth', 1.5);
semilogy(1:k_limit, state_est_error_x_aopam, "m", 'LineWidth', 1.5);
hold off;
title("State Estimation Error");
xlabel("Time step (k)");
ylabel("Error (Log Scale)");
grid on;
legend("O-PGD", "O-PAM", "AO-PGD", "AO-PAM");

% support attack error
subplot(2, 2, 2);
plot(1:k_limit, support_attack_error_opgd, "g", 'LineWidth', 1.5);
hold on;
plot(1:k_limit, support_attack_error_opam, "r", 'LineWidth', 1.5);
plot(1:k_limit, support_attack_error_aopgd, "b", 'LineWidth', 1.5);
plot(1:k_limit, support_attack_error_aopam, "m", 'LineWidth', 1.5);
hold off;
title("Support Attack Error");
xlabel("Time step (k)");
grid on;
legend("O-PGD", "O-PAM", "AO-PGD", "AO-PAM");

% false positive
subplot(2, 2, 3);
plot(1:k_limit, false_positive_opgd, "g", 'LineWidth', 1.5);
hold on;
plot(1:k_limit, false_positive_opam, "r", 'LineWidth', 1.5);
plot(1:k_limit, false_positive_aopgd, "b", 'LineWidth', 1.5);
plot(1:k_limit, false_positive_aopam, "m", 'LineWidth', 1.5);
hold off;
title("False Positive");
xlabel("Time step (k)");
grid on;
legend("O-PGD", "O-PAM", "AO-PGD", "AO-PAM");

% false negative
subplot(2, 2, 4);
plot(1:k_limit, false_negative_opgd, "g", 'LineWidth', 1.5);
hold on;
plot(1:k_limit, false_negative_opam, "r", 'LineWidth', 1.5);
plot(1:k_limit, false_negative_aopgd, "b", 'LineWidth', 1.5);
plot(1:k_limit, false_negative_aopam, "m", 'LineWidth', 1.5);
hold off;
title("False Negative");
xlabel("Time step (k)");
grid on;
legend("O-PGD", "O-PAM", "AO-PGD", "AO-PAM");