function [ f_comb ] = contact_forces( )
%This function minimizes the summation of all forces acting at points of
%contact on the hand. The cost function is subject to constraints 
%G_comb*f_comb = F_total where G_comb is = [G1, G2 ... G19]
%where each G is the grasp map of size 6x4. f_comb = [f1; f2;... f19] where 
%each f is the force vector = [fx; fy; fz; tau_z] 

close all;

G = orientations();

i = 1;

G_comb = zeros(6,76);
%Here we will manually concatenate all matrices along the 3rd dimension of
%G. So, G_comb is = [G1, G2 ... G19]
while(i <= 19)
    
    G_comb(:, 4*i-3:4*i) = G(:, :, i);
    
    i = i+1;
    
end

%This is the total torque that is imposed by the device. Note that only
%torque was acquired from simmechanics. Force was not measured under the
%impression that it would be zero
F_total = [0 5.2535 0 0.7 0 0.7]';

% %initial guess is given by (psuedo inverse of G_comb)*F_total
% %psuedo inverse of G_comb = (G_comb'*G_comb)^-1 * G_comb'
% G_inv = (G_comb' * G_comb) \ G_comb';
% f_comb0 = G_inv * F_total;

%initial guess is given by (psuedo inverse of G_comb)*F_total. The
%conventional approach for finding pseudo inverse gave negative values for
%normal force f_z and also G_inv was close to singularity. Hence the psuedo
%inverse is found using SVD. 
%G_comb = U*S*V', then pseudoInv(G_comb) = V*S_inv*U' where S_inv is found
%by taking the reciprocal of each non-zero element on the diagonal of S, 
%leaving the zeros in place, and then transposing the matrix
[U,S,V] = svd(G_comb);
S_inv = zeros(size(S));
for j = 1:6
    for k = 1:76
        if(S(j,k) ~= 0)
            S_inv(j,k) = 1/S(j,k);
        end
    end
end
S_inv = S_inv';

G_inv = V * S_inv * U';

f_comb0 = G_inv * F_total;

%f_z has to  be >=0. So we create an linear inequality constraint with
%A*f_comb = B such that A has ones in elements corresponding to f_z
%elements in f_comb
A = zeros(19,76);
for n = 1:19
    A(n, 4*n-1) = -1;
end
B = zeros(19,1);

options = optimoptions(@fmincon, 'MaxIterations', 1e6, 'MaxFunctionEvaluations', 1e6);

f_comb = fmincon(@myCost, f_comb0, A, B, G_comb, F_total, [], [], @confuneq, options);

con_force = zeros([4,19]);
for n = 1:19
    
    con_force(:,n) = f_comb(4*n-3 : 4*n);
    
end
%%
figure(1);
% subplot(3,4,1:2)
plot([46.22, 46.22+31.57, 46.22+31.57+21.67],...
    [con_force(3,1), con_force(3,2), con_force(3,3)],'k--o', 'LineWidth',1);
hold on;
plot([46.22, 46.22+31.57, 46.22+31.57+21.67],...
    [con_force(1,1), con_force(1,2), con_force(1,3)],'r--o', 'LineWidth',1);
hold on;
plot([46.22, 46.22+31.57, 46.22+31.57+21.67],...
    [con_force(2,1), con_force(2,2), con_force(2,3)],'b--o', 'LineWidth',1);
% ylim([-0.1,0.4]);
xlim([40,105]);
xlabel('Contact points along length of Digit I (mm)','FontSize',20, 'Color','k');
yyaxis left;
ylabel('Force (N)','FontSize',20, 'Color','k');
yyaxis right;
hold on;
plot([46.22, 46.22+31.57, 46.22+31.57+21.67],...
    [con_force(4,1), con_force(4,2), con_force(4,3)],'g--o', 'LineWidth',1);
ylabel('Torque (Nm)','FontSize',20, 'Color','k');
title('Digit I','FontSize',20, 'Color','k');
grid on;
leg = legend(' \it f_z',' \it f_x',' \it f_y',' \it \tau_z');
set(leg, 'FontSize', 20);
%%
figure(2);
% subplot(3,4,3:4)
plot([68.12, 68.12+39.78, 68.12+39.78+22.38, 68.12+39.78+22.38+15.82],...
    [con_force(3,4), con_force(3,5), con_force(3,6), con_force(3,7)],'k--o', 'LineWidth',1);
hold on;
plot([68.12, 68.12+39.78, 68.12+39.78+22.38, 68.12+39.78+22.38+15.82],...
    [con_force(1,4), con_force(1,5), con_force(1,6), con_force(1,7)],'r--o', 'LineWidth',1);
hold on;
plot([68.12, 68.12+39.78, 68.12+39.78+22.38, 68.12+39.78+22.38+15.82],...
    [con_force(2,4), con_force(2,5), con_force(2,6), con_force(2,7)],'b--o', 'LineWidth',1);
% ylim([-0.1,0.4]);
xlim([60,150]);
xlabel('Contact points along length of Digit II (mm)','FontSize',14, 'Color','k');
yyaxis left;
ylabel('Force (N)','FontSize',14, 'Color','k');
yyaxis right;
hold on;
plot([68.12, 68.12+39.78, 68.12+39.78+22.38, 68.12+39.78+22.38+15.82],...
    [con_force(4,4), con_force(4,5), con_force(4,6), con_force(4,7)],'g--o', 'LineWidth',1);
ylabel('Torque (Nm)','FontSize',14, 'Color','k');
title('Digit II','FontSize',14, 'Color','k');
grid on;
leg = legend(' \it f_z',' \it f_x',' \it f_y',' \it \tau_z');
set(leg, 'FontSize', 14);
%%
figure(3);
%subplot(3,4,5:6)
plot([64.60, 64.60+44.63, 64.60+44.63+26.33, 64.60+44.63+26.33+17.40],...
    [con_force(3,8), con_force(3,9), con_force(3,10), con_force(3,11)],'k--o', 'LineWidth',1);
hold on;
plot([64.60, 64.60+44.63, 64.60+44.63+26.33, 64.60+44.63+26.33+17.40],...
    [con_force(1,8), con_force(1,9), con_force(1,10), con_force(1,11)],'r--o', 'LineWidth',1);
hold on;
plot([64.60, 64.60+44.63, 64.60+44.63+26.33, 64.60+44.63+26.33+17.40],...
    [con_force(2,8), con_force(2,9), con_force(2,10), con_force(2,11)],'b--o', 'LineWidth',1);
% ylim([-0.1,0.4]);
xlim([60,160]);
xlabel('Contact points along length of Digit III (mm)','FontSize',14, 'Color','k');
yyaxis left;
ylabel('Force (N)','FontSize',14, 'Color','k');
yyaxis right;
hold on;
plot([64.60, 64.60+44.63, 64.60+44.63+26.33, 64.60+44.63+26.33+17.40],...
    [con_force(4,8), con_force(4,9), con_force(4,10), con_force(4,11)],'g--o', 'LineWidth',1);
ylabel('Torque (Nm)','FontSize',14, 'Color','k');
title('Digit III','FontSize',14, 'Color','k');
grid on;
leg = legend(' \it f_z',' \it f_x',' \it f_y',' \it \tau_z');
set(leg, 'FontSize', 14);
%%
figure(4)
%subplot(3,4,7:8)
plot([58.00, 58.00+41.37, 58.00+41.37+25.65, 58.00+41.37+25.65+17.30],...
    [con_force(3,12), con_force(3,13), con_force(3,14), con_force(3,15)],'k--o', 'LineWidth',1);
hold on;
plot([58.00, 58.00+41.37, 58.00+41.37+25.65, 58.00+41.37+25.65+17.30],...
    [con_force(1,12), con_force(1,13), con_force(1,14), con_force(1,15)],'r--o', 'LineWidth',1);
hold on;
plot([58.00, 58.00+41.37, 58.00+41.37+25.65, 58.00+41.37+25.65+17.30],...
    [con_force(2,12), con_force(2,13), con_force(2,14), con_force(2,15)],'b--o', 'LineWidth',1);
% ylim([-0.1,0.4]);
xlim([52,150]);
xlabel('Contact points along length of Digit IV (mm)','FontSize',14, 'Color','k');
yyaxis left;
ylabel('Force (N)','FontSize',14, 'Color','k');
yyaxis right;
hold on;
plot([58.00, 58.00+41.37, 58.00+41.37+25.65, 58.00+41.37+25.65+17.30],...
    [con_force(4,12), con_force(4,13), con_force(4,14), con_force(4,15)],'g--o', 'LineWidth',1);
ylabel('Torque (Nm)','FontSize',14, 'Color','k');
title('Digit IV','FontSize',14, 'Color','k');
grid on;
leg = legend(' \it f_z',' \it f_x',' \it f_y',' \it \tau_z');
set(leg, 'FontSize', 14);
%%
figure(5);
%subplot(3,4,10:11)
plot([53.69, 53.69+32.74, 53.69+32.74+18.11, 53.69+32.74+18.11+15.96],...
    [con_force(3,16), con_force(3,17), con_force(3,18), con_force(3,19)],'k--o', 'LineWidth',1);
hold on;
plot([53.69, 53.69+32.74, 53.69+32.74+18.11, 53.69+32.74+18.11+15.96],...
    [con_force(1,16), con_force(1,17), con_force(1,18), con_force(1,19)],'r--o', 'LineWidth',1);
hold on;
plot([53.69, 53.69+32.74, 53.69+32.74+18.11, 53.69+32.74+18.11+15.96],...
    [con_force(2,16), con_force(2,17), con_force(2,18), con_force(2,19)],'b--o', 'LineWidth',1);
% ylim([-0.1,0.4]);
xlim([50,125]);
xlabel('Contact points along length of Digit V (mm)','FontSize',14, 'Color','k');
yyaxis left;
ylabel('Force (N)','FontSize',14, 'Color','k');
yyaxis right;
hold on;
plot([53.69, 53.69+32.74, 53.69+32.74+18.11, 53.69+32.74+18.11+15.96],...
    [con_force(4,16), con_force(4,17), con_force(4,18), con_force(4,19)],'g--o', 'LineWidth',1);
ylabel('Torque (Nm)','FontSize',14, 'Color','k');
title('Digit V','FontSize',14, 'Color','k');
grid on;
leg = legend(' \it f_z',' \it f_x',' \it f_y',' \it \tau_z');
set(leg, 'FontSize', 14);

%%
%percentage contribution calculation toward total normal force
%calculate total normal force
normal_sum = sum(con_force(3,:));
conb = [sum(con_force(3,1:2)),mean(con_force(3,4:6)),mean(con_force(3,8:10)),...
    mean(con_force(3,12:14)),mean(con_force(3,16:18))]/normal_sum*100;
figure(6)
plot([1,2,3,4,5], conb,'k--','LineWidth',1);
plot([1,2,3,4,5], conb,'k--o');
ylabel('Percentage contribution to total force');
title('Percentage contribution of each finger to toal force');


function f_sum = myCost(f_comb)
f_sum = f_comb' * f_comb;



function [c, ceq] = confuneq(f_comb)
%coefficient of friction between skin and plastic is mu = 0.47
%torsional coefficient of friction between skin and plastic is assumed to 
%be same as mu
mu = 0.47;    
gamma = 0.47;   

c = zeros(38, 0);

%The nonlinear constraints are for each contact point f_x^2+f_y^2<=(mu*f_z)^2
%and abs(tau_z)<=gamma*f_z which are rewritten as f_x^2+f_y^2-(mu*f_z)^2<=0
%and tau_z^2-(gamma*f_z)^2<=0

n = 1;
while n <= 38

    c(n) = f_comb(2*n - 1)^2 + f_comb(2*n)^2 - (mu * f_comb(2*n + 1))^2;
    c(n+1) = f_comb(2*n + 2)^2 - (gamma * f_comb(2*n + 1))^2;
    
    n = n+2;
end

%There are no nonlinear equality constraint
ceq = [];