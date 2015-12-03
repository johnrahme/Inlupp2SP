
clear all

T = 0.1;
q = 0.2;
runtime = 10;
amount = runtime/T;
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


mean_x = [0;0];
xk_km1=zeros([2, amount]);
xk_k=zeros([2 amount+1]);
xk_k(:,1) = (mvnrnd(mean_x, Q))'; 
Pkk=zeros([2 2 amount]);
Pkkm1=zeros([2 2 amount]);
Pkk(:,:,1)=eye(2);


% The Kalman Filter
for i = 1:amount
K=Pkk(:,:,i)*Cm1'*inv(Cm1*Pkk(:,:,i)*Cm1' + R2_1); %Kalman gain
Pkkm1(:,:,i)=Pkk(:,:,i)-Pkk(:,:,i)*Cm1'*inv(Cm1*Pkk(:,:,i)*Cm1'+R2_1)*Cm1*Pkk(:,:,i); %Covariance of aposteriori estimate
xk_km1(:,i)=xk_k(:,i)+K*(y(i)-Cm1*xk_k(:,i)); %aposteriori estimate
xk_k(:,i+1)=A*xk_km1(:,i); %apriori estimate
Pkk(:,:,i+1)=A*Pkkm1(:,:,i)*A'+R1; %covariance of apriori estimate
tracerk(i) = trace(Pkkm1(:,:,i));
tracerk_1(i) = trace(Pkk(:,:,i));
end



t = T*[1:amount];
figure(1);
subplot(211), 
plot(t,x(1,1:amount),'b--',t,xk_km1(1,1:amount),'k-',t,xk_k(1,1:amount),'r-.')
grid
legend('True position','Aposteriori estimate','Apriori estimate')
ylabel('Position')
title('Position and velocity estimates using method 1')
subplot(212), 
plot(t,x(2,1:amount),'b--',t,xk_km1(2,1:amount),'k-',t,xk_k(2,1:amount),'r-.')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Apriori estimate')
xlabel('Time')

figure(2);
plot(t,tracerk,'b--',t,tracerk_1,'r-')
grid on
title('Trace of error covariance matrix using method 1')
xlabel('Time')
legend('Error covariance matrix for P_k|k-1','Error covariance matrix for P_k|k' )
disp('The limiting values for method 1:');
disp(['Limiting value if the covariance matrix for P_k|k-1 :' , num2str(tracerk(amount))]);
disp(['Limiting value if the covariance matrix for P_k|k :' , num2str(tracerk_1(amount))]);
disp(' ');
%--------------------------------------------------------------------------
%Q2

x_2=zeros([2 amount]);
x_2(:,1)=[0  1]';  % Initial state vector 
%simulation of system

for i = 1:amount
   y_2(i)= Cm2*x(:,i)+sqrt(R2_2)*randn;  
end;


xk_km1_2=zeros([2, amount]);
xk_k_2=zeros([2 amount+1]);
xk_k_2(:,1) = (mvnrnd(mean_x, Q))'; 
Pkk_2=zeros([2 2 amount]);
Pkkm1_2=zeros([2 2 amount]);
Pkk_2(:,:,1)=eye(2);

% The Kalman Filter
for i = 1:amount
K_2=Pkk_2(:,:,i)*Cm2'*inv(Cm2*Pkk_2(:,:,i)*Cm2' + R2_2); %Kalman gain
Pkkm1_2(:,:,i)=Pkk_2(:,:,i)-Pkk_2(:,:,i)*Cm2'*inv(Cm2*Pkk_2(:,:,i)*Cm2'+R2_2)*Cm2*Pkk_2(:,:,i); %Covariance of aposteriori estimate
xk_km1_2(:,i)=xk_k_2(:,i)+K_2*(y_2(i)-Cm2*xk_k_2(:,i)); %aposteriori estimate
xk_k_2(:,i+1)=A*xk_km1_2(:,i); %apriori estimate
Pkk_2(:,:,i+1)=A*Pkkm1_2(:,:,i)*A'+R1; %covariance of apriori estimate
tracerk_2(i) = trace(Pkkm1_2(:,:,i));
tracerk_1_2(i) = trace(Pkk_2(:,:,i));
end

t = T*[1:amount];
figure(3);
subplot(211), 
plot(t,x(1,1:amount),'b--',t,xk_km1_2(1,1:amount),'k-',t,xk_k_2(1,1:amount),'r-.')
grid
legend('True position','Aposteriori estimate','Apriori estimate')
ylabel('Position')
title('Position and velocity estimates using method 2')
subplot(212), 
plot(t,x(2,1:amount),'b--',t,xk_km1_2(2,1:amount),'k-',t,xk_k_2(2,1:amount),'r-.')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Apriori estimate')
xlabel('Time')

figure(4);
plot(t,tracerk_2,'b--',t,tracerk_1_2,'r-')
grid on
title('Trace of error covariance matrix using method 1')
xlabel('Time')
legend('Error covariance matrix for P_k|k-1','Error covariance matrix for P_k|k' )
disp('The limiting values for method 2:');
disp(['Limiting value if the covariance matrix for P_k|k-1 :' , num2str(tracerk_2(amount))]);
disp(['Limiting value if the covariance matrix for P_k|k :' , num2str(tracerk_1_2(amount))]);
disp(' ');
%--------------------------------------------------------------------------
%Q3
sigmaN = 0.25;
figure(5)
plot(t,tracerk,'b--', t,tracerk_2,'r-' )
title('Comparison between error covariance matrices in method 1 and method 2')
legend('Trace method 1' , 'Trace method 2')
xlabel('Time')

x_3=zeros([2 amount]);
x_3(:,1)=[0  1]';  % Initial state vector 

x_4=zeros([2 amount]);
x_4(:,1)=[0  1]';  % Initial state vector 
%simulation of system

n = eye(2)*sigmaN;
R2_1 = R2_1+n;
R2_2 = sigma1 + sigma2 + sigmaN;

for i = 1:amount
   y_3(:,i)= Cm1*x(:,i)+sqrt(R2_1)*randn(2,1);
   y_4(i)= Cm2*x(:,i)+sqrt(R2_2)*randn;  
end;


xk_km1_3=zeros([2, amount]);
xk_k_3=zeros([2 amount+1]);
xk_k_3(:,1) = (mvnrnd(mean_x, Q))'; 
Pkk_3=zeros([2 2 amount]);
Pkkm1_3=zeros([2 2 amount]);
Pkk_3(:,:,1)=eye(2);


% The Kalman Filter
for i = 1:amount
K_3=Pkk_3(:,:,i)*Cm1'*inv(Cm1*Pkk_3(:,:,i)*Cm1' + R2_1); %Kalman gain
Pkkm1_3(:,:,i)=Pkk_3(:,:,i)-Pkk_3(:,:,i)*Cm1'*inv(Cm1*Pkk_3(:,:,i)*Cm1'+R2_1)*Cm1*Pkk_3(:,:,i); %Covariance of aposteriori estimate
xk_km1_3(:,i)=xk_k_3(:,i)+K_3*(y_3(i)-Cm1*xk_k_3(:,i)); %aposteriori estimate
xk_k_3(:,i+1)=A*xk_km1_3(:,i); %apriori estimate
Pkk_3(:,:,i+1)=A*Pkkm1_3(:,:,i)*A'+R1; %covariance of apriori estimate
tracerk_3(i) = trace(Pkkm1_3(:,:,i));
tracerk_1_3(i) = trace(Pkk_3(:,:,i));
end


figure(6);
subplot(211), 
plot(t,x(1,1:amount),'b--',t,xk_km1_3(1,1:amount),'k-',t,xk_k_3(1,1:amount),'r-.')
grid
legend('True position','Aposteriori estimate','Apriori estimate')
ylabel('Position')
title('Position and velocity estimates using method 1*')
subplot(212), 
plot(t,x(2,1:amount),'b--',t,xk_km1_3(2,1:amount),'k-',t,xk_k_3(2,1:amount),'r-.')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Apriori estimate')
xlabel('Time')

figure(7);
plot(t,tracerk_3,'b--',t,tracerk_1_3,'r-')
grid on
title('Trace of error covariance matrix using method 1*')
xlabel('Time')
legend('Error covariance matrix for P_k|k-1','Error covariance matrix for P_k|k' )
disp('The limiting values for method 1*:');
disp(['Limiting value if the covariance matrix for P_k|k-1 :' , num2str(tracerk_3(amount))]);
disp(['Limiting value if the covariance matrix for P_k|k :' , num2str(tracerk_1_3(amount))]);
disp(' ');

%----------------------------------------------------

xk_km1_4=zeros([2, amount]);
xk_k_4=zeros([2 amount+1]);
xk_k_4(:,1) = (mvnrnd(mean_x, Q))'; 
Pkk_4=zeros([2 2 amount]);
Pkkm1_4=zeros([2 2 amount]);
Pkk_4(:,:,1)=eye(2);

% The Kalman Filter
for i = 1:amount
K_4=Pkk_4(:,:,i)*Cm2'*inv(Cm2*Pkk_4(:,:,i)*Cm2' + R2_2); %Kalman gain
Pkkm1_4(:,:,i)=Pkk_4(:,:,i)-Pkk_4(:,:,i)*Cm2'*inv(Cm2*Pkk_4(:,:,i)*Cm2'+R2_2)*Cm2*Pkk_4(:,:,i); %Covariance of aposteriori estimate
xk_km1_4(:,i)=xk_k_4(:,i)+K_4*(y_4(i)-Cm2*xk_k_4(:,i)); %aposteriori estimate
xk_k_4(:,i+1)=A*xk_km1_4(:,i); %apriori estimate
Pkk_4(:,:,i+1)=A*Pkkm1_4(:,:,i)*A'+R1; %covariance of apriori estimate
tracerk_4(i) = trace(Pkkm1_4(:,:,i));
tracerk_1_4(i) = trace(Pkk_4(:,:,i));
end

figure(8);
subplot(211), 
plot(t,x(1,1:amount),'b--',t,xk_km1_4(1,1:amount),'k-',t,xk_k_4(1,1:amount),'r-.')
grid
legend('True position','Aposteriori estimate','Apriori estimate')
ylabel('Position')
title('Position and velocity estimates using method 2*')
subplot(212), 
plot(t,x(2,1:amount),'b--',t,xk_km1_4(2,1:amount),'k-',t,xk_k_4(2,1:amount),'r-.')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Apriori estimate')
xlabel('Time')

figure(9);
plot(t,tracerk_4,'b--',t,tracerk_1_4,'r-')
grid on
title('Trace of error covariance matrix using method 2*')
xlabel('Time')
legend('Error covariance matrix for P_k|k-1','Error covariance matrix for P_k|k' )
disp('The limiting values for method 2*:');
disp(['Limiting value if the covariance matrix for P_k|k-1 :' , num2str(tracerk_4(amount))]);
disp(['Limiting value if the covariance matrix for P_k|k :' , num2str(tracerk_1_4(amount))]);
disp(' ');

figure(10)
plot(t,tracerk_3,'b--', t,tracerk_4,'r-' )
title('Comparison between error covariance matrices in method 1* and method 2*')
legend('Trace method 1' , 'Trace method 2')
xlabel('Time')