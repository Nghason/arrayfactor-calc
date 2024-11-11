% Định nghĩa các biến quyết định ban đầu
x = sdpvar(1, 1);
y = sdpvar(1, 1);

% Định nghĩa các biến mới
z = x + y;
j = x - y;

% Tạo các ràng buộc sử dụng các biến mới
constraints = [z <= 1, j >= 2];

% Định nghĩa hàm mục tiêu
objective = x^2 + y^2;

% Giải bài toán tối ưu hóa
diagnostics = optimize(constraints, objective);

% Kiểm tra và trích xuất kết quả
if diagnostics.problem == 0
    % Giải thành công
    optimal_x = value(x);
    optimal_y = value(y);
    disp(['Optimal x: ', num2str(optimal_x)]);
    disp(['Optimal y: ', num2str(optimal_y)]);
else
    disp('Problem could not be solved');
end
