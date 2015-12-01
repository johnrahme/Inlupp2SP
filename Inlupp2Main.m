clear all

T = 0.1;
q = 0.2;
time = 10;
amount = time/T;
sigma1 = 0.2;
sigma2 = 0.3;
C1 = [1 0];
C2 = [1 0];
Cm1 = [C1; C2];
Cm2 = C1+C2;

Q = [T^3/3 T^2/2; T^2/2 T]*q;
A = [1 T; 0 1];
B = [0 0]';
R1 = Q;
R2_1 = [sigma1 0; 0 sigma2];
R2_2 = sigma1 + sigma2;
%--------------------------------------------------------------------------
%Q1

x=zeros([2 amount]);
x(:,1)=[0  1]';  % Initial state vector 

%simulation of system
for i = 1:amount
   x(:,i+1)=A*x(:,i)+sqrt(R1)*randn(2,1);
   y(:,i)= Cm1*x(:,i)+sqrt(R2_1)*randn(2,1);
  
end;


meanx = [0;0];
xk_k=zeros([2, amount]);
xk_k1=zeros([2 amount+1]);
xk_k1(:,1) = (mvnrnd(meanx, Q))'; 
Pkk1=zeros([2 2 amount]);
Pkk=zeros([2 2 amount]);
Pkk1(:,:,1)=eye(2);


% The Kalman Filter
for i = 1:amount
K=Pkk1(:,:,i)*Cm1'*inv(Cm1*Pkk1(:,:,i)*Cm1' + R2_1); %Kalman gain
Pkk(:,:,i)=Pkk1(:,:,i)-Pkk1(:,:,i)*Cm1'*inv(Cm1*Pkk1(:,:,i)*Cm1'+R2_1)*Cm1*Pkk1(:,:,i); %Covariance of aposteriori estimate
xk_k(:,i)=xk_k1(:,i)+K*(y(i)-Cm1*xk_k1(:,i)); %aposteriori estimate
xk_k1(:,i+1)=A*xk_k(:,i); %apriori estimate
Pkk1(:,:,i+1)=A*Pkk(:,:,i)*A'+R1; %covariance of apriori estimate
tracerk(i) = trace(Pkk(:,:,i));
tracerk_1(i) = trace(Pkk1(:,:,i));
end



t = T*[1:amount];
figure(1);
subplot(211), 
plot(t,x(1,1:amount),'b--',t,xk_k(1,1:amount),'g-',t,xk_k1(1,1:amount),'r-.')
grid
legend('True position','Aposteriori estimate','Apriori estimate')
ylabel('Position')
title('Position and velocity estimates')
subplot(212), 
plot(t,x(2,1:amount),'b--',t,xk_k(2,1:amount),'g-',t,xk_k1(2,1:amount),'r-.')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Apriori estimate')
xlabel('Time')

figure(2);
plot(t,tracerk,'b--',t,tracerk_1,'r-')
grid on
legend('Error covariance matrix for P_k|k-1','Error covariance matrix for P_k|k' )
disp(['Limiting value if the covariance matrix for P_k|k-1 :' , num2str(tracerk(amount))]);
disp(['Limiting value if the covariance matrix for P_k|k :' , num2str(tracerk_1(amount))]);
%--------------------------------------------------------------------------
%Q2

x_2=zeros([2 amount]);
x_2(:,1)=[0  1]';  % Initial state vector 
%simulation of system

for i = 1:amount
   x_2(:,i+1)=A*x_2(:,i)+sqrt(R1)*randn(2,1);
   y_2(i)= Cm2*x_2(:,i)+sqrt(R2_2)*randn;  
end;


xk_k_2=zeros([2, amount]);
xk_k1_2=zeros([2 amount+1]);
xk_k1_2(:,1) = (mvnrnd(meanx, Q))'; 
Pkk1_2=zeros([2 2 amount]);
Pkk_2=zeros([2 2 amount]);
Pkk1_2(:,:,1)=eye(2);

% The Kalman Filter
for i = 1:amount
K_2=Pkk1_2(:,:,i)*Cm2'*inv(Cm2*Pkk1_2(:,:,i)*Cm2' + R2_2); %Kalman gain
Pkk_2(:,:,i)=Pkk1_2(:,:,i)-Pkk1_2(:,:,i)*Cm2'*inv(Cm2*Pkk1_2(:,:,i)*Cm2'+R2_2)*Cm2*Pkk1_2(:,:,i); %Covariance of aposteriori estimate
xk_k_2(:,i)=xk_k1_2(:,i)+K_2*(y_2(i)-Cm2*xk_k1_2(:,i)); %aposteriori estimate
xk_k1_2(:,i+1)=A*xk_k_2(:,i); %apriori estimate
Pkk1_2(:,:,i+1)=A*Pkk_2(:,:,i)*A'+R1; %covariance of apriori estimate
tracerk_2(i) = trace(Pkk_2(:,:,i));
tracerk_1_2(i) = trace(Pkk1_2(:,:,i));
end

t = T*[1:amount];
figure(3);
subplot(211), 
plot(t,x_2(1,1:amount),'b--',t,xk_k_2(1,1:amount),'g-',t,xk_k1_2(1,1:amount),'r-.')
grid
legend('True position','Aposteriori estimate','Apriori estimate')
ylabel('Position')
title('Position and velocity estimates')
subplot(212), 
plot(t,x_2(2,1:amount),'b--',t,xk_k_2(2,1:amount),'g-',t,xk_k1_2(2,1:amount),'r-.')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Apriori estimate')
xlabel('Time')

figure(4);
plot(t,tracerk_2,'b--',t,tracerk_1_2,'r-')
grid on
legend('Error covariance matrix for P_k|k-1','Error covariance matrix for P_k|k' )
disp(['Limiting value if the covariance matrix for P_k|k-1 :' , num2str(tracerk_2(amount))]);
disp(['Limiting value if the covariance matrix for P_k|k :' , num2str(tracerk_1_2(amount))]);
%--------------------------------------------------------------------------
%Q3
sigmaN = 0.25;
figure(5)
plot(t,tracerk,'b--', t,tracerk_2,'r-' )
legend('Trace method 1' , 'Trace method 2')
xlabel('Time')

x_3=zeros([2 amount]);
x_3(:,1)=[0  1]';  % Initial state vector 

x_4=zeros([2 amount]);
x_4(:,1)=[0  1]';  % Initial state vector 
%simulation of system
% n = (mvnrnd([0,0], [0.25,0;0,0.25]))';
% n_2 = mvnrnd(0, 0.25);
n = eye(2)*sigmaN;
R2_1 = R2_1+n;
R2_2 = sigma1 + sigma2 + sigmaN;

for i = 1:amount
   x_3(:,i+1)=A*x_3(:,i)+sqrt(R1)*randn(2,1);
   x_4(:,i+1)=A*x_4(:,i)+sqrt(R1)*randn(2,1);
   y_3(:,i)= Cm1*x_3(:,i)+sqrt(R2_1)*randn(2,1);
   y_4(i)= Cm2*x_4(:,i)+sqrt(R2_2)*randn;  
end;


meanx = [0;0];
xk_k_3=zeros([2, amount]);
xk_k1_3=zeros([2 amount+1]);
xk_k1_3(:,1) = (mvnrnd(meanx, Q))'; 
Pkk1_3=zeros([2 2 amount]);
Pkk_3=zeros([2 2 amount]);
Pkk1_3(:,:,1)=eye(2);


% The Kalman Filter
for i = 1:amount
K_3=Pkk1_3(:,:,i)*Cm1'*inv(Cm1*Pkk1_3(:,:,i)*Cm1' + R2_1); %Kalman gain
Pkk_3(:,:,i)=Pkk1_3(:,:,i)-Pkk1_3(:,:,i)*Cm1'*inv(Cm1*Pkk1_3(:,:,i)*Cm1'+R2_1)*Cm1*Pkk1_3(:,:,i); %Covariance of aposteriori estimate
xk_k_3(:,i)=xk_k1_3(:,i)+K_3*(y_3(i)-Cm1*xk_k1_3(:,i)); %aposteriori estimate
xk_k1_3(:,i+1)=A*xk_k_3(:,i); %apriori estimate
Pkk1_3(:,:,i+1)=A*Pkk_3(:,:,i)*A'+R1; %covariance of apriori estimate
tracerk_3(i) = trace(Pkk_3(:,:,i));
tracerk_1_3(i) = trace(Pkk1_3(:,:,i));
end


figure(6);
subplot(211), 
plot(t,x_3(1,1:amount),'b--',t,xk_k_3(1,1:amount),'g-',t,xk_k1_3(1,1:amount),'r-.')
grid
legend('True position','Aposteriori estimate','Apriori estimate')
ylabel('Position')
title('Position and velocity estimates')
subplot(212), 
plot(t,x_3(2,1:amount),'b--',t,xk_k_3(2,1:amount),'g-',t,xk_k1_3(2,1:amount),'r-.')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Apriori estimate')
xlabel('Time')

figure(7);
plot(t,tracerk_3,'b--',t,tracerk_1_3,'r-')
grid on
legend('Error covariance matrix for P_k|k-1','Error covariance matrix for P_k|k' )
disp(['Limiting value if the covariance matrix for P_k|k-1 :' , num2str(tracerk_3(amount))]);
disp(['Limiting value if the covariance matrix for P_k|k :' , num2str(tracerk_1_3(amount))]);

%----------------------------------------------------

xk_k_4=zeros([2, amount]);
xk_k1_4=zeros([2 amount+1]);
xk_k1_4(:,1) = (mvnrnd(meanx, Q))'; 
Pkk1_4=zeros([2 2 amount]);
Pkk_4=zeros([2 2 amount]);
Pkk1_4(:,:,1)=eye(2);

% The Kalman Filter
for i = 1:amount
K_4=Pkk1_4(:,:,i)*Cm2'*inv(Cm2*Pkk1_4(:,:,i)*Cm2' + R2_2); %Kalman gain
Pkk_4(:,:,i)=Pkk1_4(:,:,i)-Pkk1_4(:,:,i)*Cm2'*inv(Cm2*Pkk1_4(:,:,i)*Cm2'+R2_2)*Cm2*Pkk1_4(:,:,i); %Covariance of aposteriori estimate
xk_k_4(:,i)=xk_k1_4(:,i)+K_4*(y_4(i)-Cm2*xk_k1_4(:,i)); %aposteriori estimate
xk_k1_4(:,i+1)=A*xk_k_4(:,i); %apriori estimate
Pkk1_4(:,:,i+1)=A*Pkk_4(:,:,i)*A'+R1; %covariance of apriori estimate
tracerk_4(i) = trace(Pkk_4(:,:,i));
tracerk_1_4(i) = trace(Pkk1_4(:,:,i));
end

figure(8);
subplot(211), 
plot(t,x_4(1,1:amount),'b--',t,xk_k_4(1,1:amount),'g-',t,xk_k1_4(1,1:amount),'r-.')
grid
legend('True position','Aposteriori estimate','Apriori estimate')
ylabel('Position')
title('Position and velocity estimates')
subplot(212), 
plot(t,x_4(2,1:amount),'b--',t,xk_k_4(2,1:amount),'g-',t,xk_k1_4(2,1:amount),'r-.')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Apriori estimate')
xlabel('Time')

figure(9);
plot(t,tracerk_4,'b--',t,tracerk_1_4,'r-')
grid on
legend('Error covariance matrix for P_k|k-1','Error covariance matrix for P_k|k' )
disp(['Limiting value if the covariance matrix for P_k|k-1 :' , num2str(tracerk_4(amount))]);
disp(['Limiting value if the covariance matrix for P_k|k :' , num2str(tracerk_1_4(amount))]);

figure(10)
plot(t,tracerk_3,'b--', t,tracerk_4,'r-' )
legend('Trace method 1' , 'Trace method 2')
xlabel('Time')
