numPlotsPerFigure = 10;
subplotRows = 5;
subplotCols = 2;
% Gọi hàm cho từng ma trận
createMultiSubplotFigure(1, IR_12, 'IR_{12}', numPlotsPerFigure, subplotRows, subplotCols, tau);
createMultiSubplotFigure(2, IR_3, 'IR_3', numPlotsPerFigure, subplotRows, subplotCols, tau);
createMultiSubplotFigure(3, testr, 'r', 3, subplotRows, subplotCols, tau);
disp('Đã hoàn tất việc vẽ biểu đồ.');






function createMultiSubplotFigure(figureHandle, dataMatrix, baseTitle, numPlots, rows, cols, tau)
    figure(figureHandle); % Chọn hoặc tạo figure được chỉ định
    sgtitle(['Biểu đồ cho ', baseTitle]); % Tiêu đề chung

    for i = 1:numPlots
        subplot(rows, cols, i);
        plot( real(dataMatrix(:, i)));
    end

    for i = 1:numPlots
        subplot(rows, cols, i);
        plot( abs(dataMatrix(:, i)));
    end
end