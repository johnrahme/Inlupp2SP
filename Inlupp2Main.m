T = 0.1;
q = 0.2;
sigma1 = 0.2;
sigma2 = 0.3;
C1 = [1 0];
C2 = [1 0];
Cm1 = [C1; C2];
Cm2 = C1+C2;

Q = [T^3/3 T^2/2; T^2/2 T]*q;
A = [1 T; 0 1];
B = [0 0]';
Dm1 = [0 0]';
Dm2 = 0;

Rm1 = [sigma1 0; 0 sigma2];
Rm2 = sigma1 + sigma2;
Nm1 = [0 0; 0 0];
Nm2 = 0;


system1 = ss(A, B, Cm1, Dm1, T);
system2 = ss(A, B, Cm2, Dm2, T);

k1 = kalman(system1, Q(1), Rm1, Nm1(1)); % feeeeeeeeeeeeeeeeeeel omfg we suck

