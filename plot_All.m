%% Status of figure
load("VA.mat","used_VAs");
load("pos.mat","pos_centers");
hold on;
title('Điểm thu phát');
grid("on");
xlabel('x - axis');
ylabel('y- axis');

%% Plot Pos_centers Điểm thu
%plot(pos_centers(1,:), pos_centers(2,:), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
plot(pos_centers(1,11:20), pos_centers(2,11:20), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
%% Plot used_VAs Điểm phát
VAtmp = used_VAs{1,1}(:, 1);
VAtmp(1, 1) = -used_VAs{1,1}(1, 1);
plot(VAtmp(1, 1),used_VAs{1,1}(2, 1),  'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
legend({'Pos center (Điểm thu)','used VAs \{1,1\} (Điểm phát)'});


%%
clear VAtmp;