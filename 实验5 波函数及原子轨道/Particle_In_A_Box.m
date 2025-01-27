% MATLAB script to calculate the particle-in-a-box problem
%
%

figure

% Box length
L = 10;

% Interval
d = 0.1;

% X axis: 101 points from 0 to L
x = 0 : d : L;

% Two candidate wavefunctions
psi1 = sin(pi * x / L);
psi2 = sin(3 * pi * x / L);

plot(x, psi1, 'r', x, psi2, 'b');
title('Two candidate wavefunctions');
legend('\psi_1', '\psi_2','Location','SouthEast');

pause

% Probability = amplitude ^ 2
P1 = psi1 .^ 2;
P2 = psi2 .^ 2;

plot(x, psi1, 'r--', x, psi2, 'b--', x, P1, 'r', x, P2, 'b');

% Total probabilities
S1 = sum(P1) * d;  
S2 = sum(P2) * d; 

%title(sprintf('Probabilities S1 = %.2f, S2 = %.2f', S1, S2));
title('Probabilities S_1 = ?, S_2 = ?');
legend('\psi_1', '\psi_2','|\psi_1|^2', '|\psi_2|^2','Location','SouthEast');

pause

p1 = patch(x, P1, x*0);
set(p1,'FaceColor',[1 0.5 0.5],'EdgeColor',[1 0 0]);

pause

p2 = patch(x, P2, x*0);
set(p2,'FaceColor',[0.5 0.5 1],'EdgeColor',[0 0 1]);

pause

% Normalization
psi3 = psi1 / sqrt(S1);
psi4 = psi2 / sqrt(S2);

plot(x, psi1, 'r--', x, psi2, 'b--', x, psi3, 'r', x, psi4, 'b');
axis([0, L, -1, 1]);
title('Two normalized wavefunctions');
legend('\psi_1', '\psi_2','\psi_3', '\psi_4','Location','SouthEast');

pause

% Calculate probabilities again
P3 = psi3 .^ 2;
P4 = psi4 .^ 2;

plot(x, psi3, 'r--', x, psi4, 'b--', x, P3, 'r', x, P4, 'b');

% Total probabilities
S3 = sum(P3) * d;  
S4 = sum(P4) * d;

axis([0, L, -1, 1]);
%title(sprintf('Probabilities S3 = %.2f, S4 = %.2f', S3, S4));
title('Probabilities S_3 = ?, S_4 = ?');
legend('\psi_3', '\psi_4','|\psi_3|^2', '|\psi_4|^2','Location','SouthEast');

pause

p1 = patch(x, P3, x*0);
set(p1,'FaceColor',[1 0.5 0.5],'EdgeColor',[1 0 0]);

pause

p2 = patch(x, P4, x*0);
set(p2,'FaceColor',[0.5 0.5 1],'EdgeColor',[0 0 1]);

pause

% Calculate 1st and 2nd derivatives of psi3
psi3_1 = diff(psi3) / d;
psi3_2 = diff(psi3_1) / d;

h = subplot(3, 1, 1);
plot(h, x, psi3);
legend(h, '\psi_3', 'Location','NorthEast');
h = subplot(3, 1, 2);
plot(h, x(2:end), psi3_1);
legend(h, '\psi_3''', 'Location','NorthEast');
h = subplot(3, 1, 3);
plot(h, x(2:end-1), psi3_2);
xlabel(h, 'x');
legend(h, '\psi_3"', 'Location','SouthEast');

pause

% Energy
E = 0.00125;

h = subplot(1, 2, 1);
plot(h, x(2:end-1), -psi3_2 / 8 / pi^2);
xlabel(h, 'x');
legend(h, 'H\psi_3', 'Location','South');
h = subplot(1, 2, 2);
plot(h, x, psi3);
xlabel(h, 'x');
legend(h, '\psi_3', 'Location','South');

clearvars


