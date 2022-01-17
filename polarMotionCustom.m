function [ys] = polarMotionCustom(mjd)
[numbers, string, raw] = xlsread("./data ERP.xlsx");
mjd_index = find(numbers==mjd);
y_x_test_pred = numbers(mjd_index, 4);
y_y_test_pred = numbers(mjd_index, 5);
ys = [y_x_test_pred y_y_test_pred];
end

