%% TODO Localization

clear;
clc;
load('Work_Space_Localization_Short_project.mat');

%% 1. Plot calculated Pose Estimation with Encoder data

for compress=1:1

global ini;
global auxCV;

pCalc = zeros(522,3);
V=zeros(2,2);
pkn=pk.signals.values(:,:,1);

ini = 1;
CoV = zeros(3,3,524);

%PoseEstimation calculus
for index=1:522    
    
    p=PoseEstimation(data_enc(index,6),data_enc(index,7),V,pkn);
    CoV(:,:,index) = auxCV;
    
    x=p(1);
    y=p(2);
    th=p(3);
    
    pCalc(index,1) = x;
    pCalc(index,2) = y;
    pCalc(index,3) = th;
end

%Display
for index=1:522 
    
    plot (pCalc(:,1), pCalc(:,2), 'r.','LineWidth',1.5) 
    
    % Plotting the covariance matrix
    hold on;
    plot_ellipse(CoV(1:2,1:2,index),[pCalc(:,1), pCalc(:,2)],'g'); 
    
    pause(0.1);
    %clf
end

end

%% 2 Ploting Land Marks in Robot Reference Frame 

for compress=1:1

t = 0: 2*pi/359 : 2*pi;
laserRx = zeros(522,360);
laserRy = zeros(522,360);


for i=1:522
    [X,Y] = pol2cart(t,lds_dis (i,2:361)); 
    
    laserRx(i,:)=X;
    laserRy(i,:)=Y;
    
    hold on
    sprintf("iter");
    scatter(X,Y);
end

end
    
%% 3 Plot the Land Marks in World Reference Frame

for compress=1:1
    
    TWR = zeros(4,4,522);
    laserW = zeros(522,2);
   
   for i=1:522
       
       x = pCalc(i,1);
       y = pCalc(i,2);
       
       th=pCalc(i,3);
       s =sin(th);
       c = cos(th);
       
       TWR(:,:,i) = ...
           [c -s 0 x;...
            s c  0 y;...
            0 0 1 0;...
            0 0 0 1];
        
        %Transformation
        % theta for all == 0
        
        Q = [laserRx(i) laserRy(i) 0 1]';
        R=TWR(:,:,i)*Q;
        laserW(i,:)=R(1:3);        
        
   end
    
end

%% 4 Associated	Land Mark.

%% 5 Similarity	Transform.Adapt	the ST to output the error in pose given a time.

%% Functions 

function [Pose_t,Pose_est,Pk] = PoseEstimation(L,R,V,pk0) 

global ini;
global auxCV;

persistent x_w y_w suma_theta  
persistent Pose_est1 Pose_t1 Pk1

S = 243;

%%%%%%%%%%%%%%%   variable initialitation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (ini == 1)
    ini = 0;
    x_w=0;
    y_w=0;
    suma_theta=0; 
    Pose_est1 = [x_w;y_w;suma_theta];
    Pose_t=[x_w;y_w;suma_theta]
    Pk1 = pk0;
end

%%%%%%%%%%%%%%%%%  compute the odometry    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_th = ((R-L)/(2*S));
delta_d = (R+L)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  integration of the robot pose %%%%%%%%%%%%%%%%%%%%%%%%%%

x_w = x_w + (delta_d + V(1,1))*cos(suma_theta + delta_th + V(2,2))
y_w = y_w + (delta_d + V(1,1))*sin(suma_theta + delta_th + V(2,2))
suma_theta=mod((suma_theta + delta_th + V(2,2)) ,2*pi); 

Pose_t=[x_w;y_w;suma_theta];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Build the Jacocian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3

F_x = [[1 0 -(delta_d*sin(suma_theta+delta_th))];....
      [0 1 -(delta_d*cos(suma_theta+delta_th))];....
      [0 0 1]];

F_v = [[cos(suma_theta+delta_th) -delta_d*sin(suma_theta+delta_th)];....
      [sin(suma_theta+delta_th) delta_d*cos(suma_theta+delta_th)];[0 1]];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Compute the EKF Trucated Taylor %%%%%%%%%%%%%%%%%%%%%%%%

 Pose_est1 = Pose_est1 + F_x*(Pose_t - Pose_est1) + F_v*diag(V);
 Pk1 = F_x*Pk1*F_x' + F_v*V*F_v';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  output variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Pose_t=[x_w;y_w;suma_theta];
Pose_est= Pose_est1;
Pk=Pk1;
auxCV = Pk1;

end

function [tx,ty,tita,s]=SimilarityTransform(LandMark,newLM)
%     %assert(size( LandMark , 2) == size(detected, 2));
% 
%     %Build Matrix A
%     A = [];
%     for i=1:size( LandMark , 2)
%         A = [A;[ LandMark (1,i), LandMark (2,i),1,0]];
%         A = [A;[ LandMark (2,i),-LandMark (1,i),0,1]];
%     end
% 
%     %Build Matrix B
%     B = [];
%     for i=1:size( newLM , 2)
%         B = [B; newLM (1,i); newLM (2,i)];
%     end
% 
%     %Compute tx ty i tita
%     X = inv((A'*A))*A'*B;
%     Tx_ST = X(3);
%     Ty_ST= X(4);
%     alpha_ST = atan2(X(2),X(1))*180/pi;
end
