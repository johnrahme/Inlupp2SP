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


A = [1 T; 0 1];
B = [0 0]';
Dm1 = [0 0]';
Dm2 = 0;
R1 = [T^3/3 T^2/2; T^2/2 T]*q;
R2_1 = [sigma1 0; 0 sigma2];
R2_2 = sigma1 + sigma2;
Nm1 = [0 0; 0 0];
Nm2 = 0;

x=zeros([2 amount]);
x(:,1)=[0  1]';  % Initial state vector 

%simulation of system
for i = 1:amount
   x(:,i+1)=A*x(:,i)+Dm1*sqrt(q)*randn;
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
K=Pkk1(:,:,i)*Cm1'*inv(Cm1*Pkk1(:,:,i)*Cm1' + R1); %Kalman gain
Pkk(:,:,i)=Pkk1(:,:,i)-Pkk1(:,:,i)*Cm1'*inv(Cm1*Pkk1(:,:,i)*Cm1'+R1)*Cm1*Pkk1(:,:,i); %Covariance of aposteriori estimate
xk_k(:,i)=xk_k1(:,i)+K*(y(i)-Cm1*xk_k1(:,i)); %aposteriori estimate
xk_k1(:,i+1)=A*xk_k(:,i); %apriori estimate
Pkk1(:,:,i+1)=A*Pkk(:,:,i)*A'+R2_1; %covariance of apriori estimate
end



t = T*[1:amount];
figure(1);
subplot(211), 
plot(t,x(1,1:amount),'b--',t,y(1:amount),'m-',t,xk_k(1,:),'r-.')
grid
legend('True position','Measured position','Aposteriori estimate')
ylabel('Position')
title('Position and velocity estimates')
subplot(212), 
plot(t,x(2,1:amount),'b--',t,xk_k(2,:),'r-.')
grid
ylabel('Velocity')
legend('True velocity','Aposteriori estimate')
xlabel('Time')

