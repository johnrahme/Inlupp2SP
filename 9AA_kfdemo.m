echo on
T=0.1; % Sampling interval
simtime=10; % Simulation time
samples=simtime/T;	%number of samples
x=zeros([2, samples+1]);
A=[1 T;0 1]; % System Matrix
D=[T^2/2 T]'; % System Noise Matrix
C=[1 0]; % Measurement Matrix
q=0.25;	% Process noise variance
R1=D*q*D';  % Noise covariance
r2=0.25; % Measurement noise 
x=zeros([2 samples]);
x(:,1)=[0  1]';  % Initial state vector 

%simulation of system
for i = 1:samples
   x(:,i+1)=A*x(:,i)+D*sqrt(q)*randn;
   y(i)= C*x(:,i)+sqrt(r2)*randn;
   echo off
end;

echo on

% Initialization of the Kalman Filter
pause 
xk_k=zeros([2, samples]);
xk_k1=zeros([2 samples+1]);
xk_k1(:,1) = [-2 0]'; % initial apriori state estimate 
Pkk1=zeros([2 2 samples]);
Pkk=zeros([2 2 samples]);
Pkk1(:,:,1)=eye(2);
% The Kalman Filter
pause
for i=1:samples,
   K(:,i)=Pkk1(:,:,i)*C'*inv(C*Pkk1(:,:,i)*C' + r2); %Kalman gain
   Pkk(:,:,i)=Pkk1(:,:,i)-Pkk1(:,:,i)*C'*inv(C*Pkk1(:,:,i)*C'+r2)*C*Pkk1(:,:,i); %Covariance of aposteriori estimate
   xk_k(:,i)=xk_k1(:,i)+K(:,i)*(y(i)-C*xk_k1(:,i)); %aposteriori estimate
   xk_k1(:,i+1)=A*xk_k(:,i); %apriori estimate
   Pkk1(:,:,i+1)=A*Pkk(:,:,i)*A'+R1; %covariance of apriori estimate
   echo off
end
echo on
%plotting of variables

pause
t = T*[1:samples];
figure(1);
subplot(211), 
plot(t,x(1,1:samples),'b--',t,y(1:samples),'m-',t,xk_k(1,:),'r-.')
grid
legend('True position','Measured position','Aposteriori estimate')
ylabel('Position')
title('Position and velocity estimates')
subplot(212), 
plot(t,x(2,1:samples),'b--',t,xk_k(2,:),'r-.')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate')
xlabel('Time')

%plotting of K(k)
%pause

figure(2)
plot(t,K(1,:),'-',t,K(2,:),'-.')
title('Kalman gain')
xlabel('Time')
legend('Position gain','Velocity gain')

%assumed variance of measurement noise smaller than it is

r2=0.0025;

pause

echo on
% Initialization of the Kalman Filter
%pause 
xk2_k=zeros([2, samples]);
xk2_k1=zeros([2 samples+1]);
xk2_k1(:,1) = [-2 0]'; % initial apriori state estimate 
Pkk21=zeros([2 2 samples]);
Pkk2=zeros([2 2 samples]);
Pkk21(:,:,1)=eye(2);
% The Kalman Filter
for i=1:samples,
   K2(:,i)=Pkk21(:,:,i)*C'*inv(C*Pkk21(:,:,i)*C' + r2); %Kalman gain
   Pkk2(:,:,i)=Pkk21(:,:,i)-Pkk21(:,:,i)*C'*inv(C*Pkk21(:,:,i)*C'+r2)*C*Pkk21(:,:,i); %Covariance of aposteriori estimate
   xk2_k(:,i)=xk2_k1(:,i)+K2(:,i)*(y(i)-C*xk2_k1(:,i)); %aposteriori estimate
   xk2_k1(:,i+1)=A*xk2_k(:,i); %apriori estimate
   Pkk21(:,:,i+1)=A*Pkk2(:,:,i)*A'+R1; %covariance of apriori estimate
   echo off
end

echo on

%plotting of estimates
pause
t = T*[1:samples];

figure(3);

subplot(211), 
plot(t,x(1,1:samples),'b--',t,y(1:samples),'m-',t,xk_k(1,:),'r-.',t,xk2_k(1,:),'g-')
grid
legend('True position','Measured position','Aposteriori estimate, correct variance','Aposteriori estimate, too small variance')
ylabel('Position')
title('Position and velocity estimates')
subplot(212), 
plot(t,x(2,1:samples),'b--',t,xk_k(2,:),'r-.',t,xk2_k(2,:),'g-')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate, correct variance','Aposteriori estimate, too small variance')
xlabel('Time')

%plotting of K(k)


figure(4)

subplot(111)
plot(t,K(1,:),'-',t,K(2,:),'-.',t,K2(1,:),'--',t,K2(2,:),'-')
title('Kalman gain')
xlabel('Time')
legend('Position gain','Velocity gain','Pos gain, too small variance','Vel gain, too small variance') 

%assumed variance of measurement noise larger than it is

r2=25;
pause

echo on
% Initialization of the Kalman Filter

xk2_k=zeros([2, samples]);
xk2_k1=zeros([2 samples+1]);
xk2_k1(:,1) = [-2 0]'; % initial apriori state estimate 
Pkk21=zeros([2 2 samples]);
Pkk2=zeros([2 2 samples]);
Pkk21(:,:,1)=eye(2);
% The Kalman Filter
for i=1:samples,
   K2(:,i)=Pkk21(:,:,i)*C'*inv(C*Pkk21(:,:,i)*C' + r2); %Kalman gain
   Pkk2(:,:,i)=Pkk21(:,:,i)-Pkk21(:,:,i)*C'*inv(C*Pkk21(:,:,i)*C'+r2)*C*Pkk21(:,:,i); %Covariance of aposteriori estimate
   xk2_k(:,i)=xk2_k1(:,i)+K2(:,i)*(y(i)-C*xk2_k1(:,i)); %aposteriori estimate
   xk2_k1(:,i+1)=A*xk2_k(:,i); %apriori estimate
   Pkk21(:,:,i+1)=A*Pkk2(:,:,i)*A'+R1; %covariance of apriori estimate
   echo off
end

echo on

%plotting of estimates
pause
t = T*[1:samples];

figure(5);

subplot(211), 
plot(t,x(1,1:samples),'b--',t,y(1:samples),'m-',t,xk_k(1,:),'r-.',t,xk2_k(1,:),'g-')
grid
legend('True position','Measured position','Aposteriori estimate','Too large variance')
ylabel('Position')
title('Position and velocity estimates. Assume variance of measurement noise too large.')
subplot(212), 
plot(t,x(2,1:samples),'b--',t,xk_k(2,:),'r-.',t,xk2_k(2,:),'g-')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Too large variance')
xlabel('Time')

%plotting of K(k)


figure(6)

subplot(111)
plot(t,K(1,:),'-',t,K(2,:),'-.',t,K2(1,:),'--',t,K2(2,:),'-')
title('Kalman gain. Assume variance of measurement noise too large')
xlabel('Time')
legend('Position gain','Velocity gain','Pos gain, too large variance','Vel gain, too large variance')



%assumed variance of initial state smaller than it is

r2=0.25;

echo on
% Initialization of the Kalman Filter
%pause 
xk2_k=zeros([2, samples]);
xk2_k1=zeros([2 samples+1]);
xk2_k1(:,1) = [-2 0]'; % initial apriori state estimate 
Pkk21=zeros([2 2 samples]);
Pkk2=zeros([2 2 samples]);
Pkk21(:,:,1)=0.01*eye(2);
pause
% The Kalman Filter
for i=1:samples,
   K2(:,i)=Pkk21(:,:,i)*C'*inv(C*Pkk21(:,:,i)*C' + r2); %Kalman gain
   Pkk2(:,:,i)=Pkk21(:,:,i)-Pkk21(:,:,i)*C'*inv(C*Pkk21(:,:,i)*C'+r2)*C*Pkk21(:,:,i); %Covariance of aposteriori estimate
   xk2_k(:,i)=xk2_k1(:,i)+K2(:,i)*(y(i)-C*xk2_k1(:,i)); %aposteriori estimate
   xk2_k1(:,i+1)=A*xk2_k(:,i); %apriori estimate
   Pkk21(:,:,i+1)=A*Pkk2(:,:,i)*A'+R1; %covariance of apriori estimate
   echo off
end

echo on

%plotting of estimates
pause

t = T*[1:samples];
figure(7);
subplot(211), 
plot(t,x(1,1:samples),'b--',t,y(1:samples),'m-',t,xk_k(1,:),'r-.',t,xk2_k(1,:),'g-')
grid
legend('True position','Measured position','Aposteriori estimate','Too small initial variance')
ylabel('Position')
title('Position and velocity estimates. Too small variance of initial state')
subplot(212), 
plot(t,x(2,1:samples),'b--',t,xk_k(2,:),'r-.',t,xk2_k(2,:),'g-')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Too small initial variance')
xlabel('Time')

%plotting of K(k)


figure(8)
subplot(111)
plot(t,K(1,:),'-',t,K(2,:),'-.',t,K2(1,:),'--',t,K2(2,:),'-')
title('Kalman gain, Too small initial variance')
xlabel('Time')
legend('Position gain','Velocity gain','Pos gain, Too small initial variance','Vel gain, Too small initial variance')



%Stationary Kalman filter
r2=0.25;

Kstat=K2(:,samples);
pause
echo on
% Initialization of the Kalman Filter
%pause 
xk2_k=zeros([2, samples]);
xk2_k1=zeros([2 samples+1]);
xk2_k1(:,1) = [-2 0]'; % initial apriori state estimate 
% The Kalman Filter
%pause
for i=1:samples,
   xk2_k(:,i)=xk2_k1(:,i)+Kstat*(y(i)-C*xk2_k1(:,i)); %aposteriori estimate
   xk2_k1(:,i+1)=A*xk2_k(:,i); %apriori estimate
   echo off
end

echo on

%plotting of estimates
pause

t = T*[1:samples];
figure(9);
subplot(211), 
plot(t,x(1,1:samples),'b--',t,y(1:samples),'m-',t,xk_k(1,:),'r-.',t,xk2_k(1,:),'g-')
grid
legend('True position','Measured position','Aposteriori estimate','Stationary KF')
ylabel('Position')
title('Position and velocity estimates. Stationary Kalman filter')
subplot(212), 
plot(t,x(2,1:samples),'b--',t,xk_k(2,:),'r-.',t,xk2_k(2,:),'g-')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate','Stationary KF')
xlabel('Time')


%plotting of K(k)


figure(10)
subplot(111)
plot(t,K(1,:),'-',t,K(2,:),'-.',t,Kstat(1)*ones([samples, 1]),'--',t,Kstat(2)*ones([samples, 1]),'-')
title('Kalman gain')
xlabel('Time')
legend('Position gain','Velocity gain','Pos gain, stationary KF','Vel gain, stationary KF')

echo off


