load('Work_Space_Localization_Short_project.mat');

%Introducing the snapshot to visualize
x = inputdlg('Enter step time to visualize','Input', [1 20]);
index = str2num(x{:});

%Robot is a triangle
Robot= [0 -0.2 0 1;0.4 0 0 1;0 0.2 0 1]';

for index=1:522 
    t = 0: 2*pi/359 : 2*pi;
    
    %Polar 
    subplot(1,2,1);
    P = polar(t, 4.5 * ones(size(t))); % to fix the limits
    set(P, 'Visible', 'off')
    polar(t, lds_dis (index,2:361), '--g'); % Ploting the laser data wrt Robot frame
    title ('Laser data at Robot Reference Frame','FontWeight','bold','FontSize',16)
   
    
    subplot(1,2,2);
    title ('Data on Wordl Reference Frame', 'FontWeight','bold','FontSize',16)
    axis([-3 3 -2 4])
    grid on
    hold on
    for i=1:4 % plotting the 4 Land Marks
        circle (LandMark(i,:)',0.15)
    end
    % plotting the land mark seen by the Robot wrt wordl reference frame
    scatter(ldx(index,:), ldy(index,:)) 
    % Plotting the trajectory
    plot (trajec(:,1), trajec(:,2), 'r.','LineWidth',1.5) 
    
    % moving the robot
    Robot_tr=transl(trajec(index,1),trajec(index,2),0)*trotz(mod(trajec(index,3)+pi/2,2*pi))*Robot;
    
    patch(Robot_tr(1,:), Robot_tr(2,:),'b');
    % Plotting the covariance matrix
    plot_ellipse(pk.signals.values(1:2,1:2,index),[trajec(index,1),trajec(index,2)],'g'); 
    pause(0.1);
    clf
end