% TODO


load('Work_Space_Localization_Short_project.mat');

%Introducing the snapshot to visualize
x = inputdlg('Enter step time to visualize','Input', [1 20]);
index = str2num(x{:});

x_w = 0;
y_w = 0;
suma_theta = 0;

%Robot is a triangle
Robot= [0 -0.2 0 1;0.4 0 0 1;0 0.2 0 1]';

%Simplify Data Encoders
newP = zeros(523,3);

%Initial position
L  = data_enc(:,[1,6]);
R  = data_enc(:,[1,7]);
th = trajec(:,3);

for index=1:522 
    t = 0: 2*pi/359 : 2*pi;
    
    % 1. Pose	estimated.    
    subplot(1,2,2);
    title ('Data on Wordl Reference Frame', 'FontWeight','bold','FontSize',16)
    axis([-3 3 -2 4])
    grid on
    hold on
    
    % plotting the 4 Land Marks
    for i=1:4 
        circle (LandMark(i,:)',0.15)
    end
    
    %Calculate new pose with poseintegration%%%%%%

    delta_th= (R(index,2)-L(index,2))/(2*width);
    delta_d = (R(index,2)+L(index,2))/2;

    x_w = x_w + (delta_d*cos(suma_theta));
    y_w = y_w + (delta_d*sin(suma_theta));
    suma_theta=mod(suma_theta+delta_th,2*pi);
    
    trajec(index,1) = x_w;    
    trajec(index,2) = y_w;
    trajec(index,3) = suma_theta;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % PreCalculated LanMarks seen by the Robot wrt wordl reference frame
    scatter(ldx(index,:), ldy(index,:)) 
    % Plotting the trajectory
    plot (trajec(:,1), trajec(:,2), 'r.','LineWidth',1.5) 
    
    % moving the robot
    Robot_tr=transl(trajec(index,1),trajec(index,2),0)*trotz(mod(trajec(index,3)+pi/2,2*pi))*Robot;    
    patch(Robot_tr(1,:), Robot_tr(2,:),'b');
    
    % Plotting the covariance matrix
    plot_ellipse(pk.signals.values(1:2,1:2,index),[trajec(index,1),trajec(index,2)],'g'); 
    
    
    % 2. Polar to cartesian
    subplot(1,2,1);
    P = polar(t, 4.5 * ones(size(t))); % to fix the limits
    set(P, 'Visible', 'off')
    polar(t, lds_dis (index,2:361), '--g') % Ploting the laser data wrt Robot frame
    title ('Laser data at Robot Reference Frame','FontWeight','bold','FontSize',16)
    
    pause(0.1);
    clf
end





