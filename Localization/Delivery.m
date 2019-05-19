%% TODO Localization

clear;
clc;
load('Work_Space_Localization_Short_project.mat');

%% 1. Plot calculated Pose Estimation with Encoder data

for compress=1:1
    
POSE = zeros(522,3); 

x_w=0;
y_w=0;
suma_theta=pi/2; 

CoV = zeros(3,3,524);

%PoseEstimation calculus
for i=1:522        
    
    dL = data_enc(i+1,6)-data_enc(i,6);
    dR = data_enc(i+1,7)-data_enc(i,7);
    
    delta_th = ((dR-dL)/(2*width));
    delta_d = (dR+dL)/2;

    x_w = x_w + (delta_d*cos(suma_theta));
    y_w = y_w + (delta_d*sin(suma_theta));
    suma_theta=suma_theta + delta_th;
   
    POSE(i,:) = [x_w/1000 y_w/1000 suma_theta];
end

%Display
title("Trajectoria creada")
plot (POSE(:,1), POSE(:,2), 'r.','LineWidth',1.5) % Plotting the trajectory estimated



%Display
% for index=1:522 
%     
%     plot (pCalc(:,1), pCalc(:,2), 'r.','LineWidth',1.5) 
%     
%     % Plotting the covariance matrix
%     hold on;
%     plot_ellipse(CoV(1:2,1:2,index),[pCalc(:,1), pCalc(:,2)],'g'); 
%     
%     pause(0.1);
%     %clf
% end

end

%% 2 Ploting Land Marks in Robot Reference Frame 

for compress=1:1

t = 0: 2*pi/359 : 2*pi;
laserRx = zeros(522,360);
laserRy = zeros(522,360);


for i=1:522
    [X,Y] = pol2cart(t,lds_dis (i,2:361)); 
    
    laserRx(i,:)=X/1000;
    laserRy(i,:)=Y/1000;
    
    hold on
    scatter(laserRx(1,:),laserRy(1,:));
end

end
    
%% 3 Plot the Land Marks in World Reference Frame

for compress=1:1
    
    TWR = zeros(4,4,522);
    laserWx = zeros(522,360);
    laserWy = zeros(522,360);
   
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
        for j=1:size(laserRx,2)
            
            Q = [laserRx(i,j) laserRy(i,j) 0 1]';
            R=TWR(:,:,i)*Q;
            
            laserWx(i,j)=R(1,1); 
            laserWy(i,j)=R(2,1);
        end
        
   end
    
for i=1:522
    hold on
    scatter(laserWx(1,:),laserWy(1,:));    
end


end

%% 4 Associated	Land Mark.



%% 5 Similarity	Transform.Adapt	the ST to output the error in pose given a time.

%% Functions 

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
