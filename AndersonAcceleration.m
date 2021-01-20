% --------------------------------------------------
% Anderson 加速模块
% Author:   Xiaowei Jia
% Date:     2020-01-08
% Version:  0.1
% TODO: 
%           Damping 模块;
%           条件数监控模块;
%           不同类型残差定义模块;
% --------------------------------------------------

function AA = AndersonAcceleration()
    % --------------------------------------------------
    % ----  Initialize: 初始化向量部分
    % ----  setSpaceSize: 设置 Anderson 加速深度
    % ----  updateSolution: Anderson 加速更新近似解
    % --------------------------------------------------
    AA.Initialize = @Initialize;
    AA.setSpaceSize = @setSpaceSize;
    AA.updateSolution = @updateSolution;
end

function [] = Initialize();
    % ----  Initialize: 初始化向量部分
    % 需要保留历史迭代信息, 因此讲涉及到的变量定义为 global 类型.
    % 主程序运行前需要 clear all 清除 global 变量.

    % Anderson 加速深度以及随迭代指标变化参数
    global mMax mAA;
    mAA = 0;

    % Anderson 最小二乘问题 Q R 分解
    global Q R;

    % 存储历史向量
    global g_old f_old;
    global DG;

end

function [] = setSpaceSize(m)
    % ----  setSpaceSize: 设置 Anderson 加速深度
    global mMax;
    mMax = m;
end

function x = updateSolution(iter, x, gval)
    % ----  updateSolution: Anderson 加速更新近似解
    % 输入参数:
    %   iter:   当前迭代指标 (最外层)
    %   x:      上一次近似解
    %   gval:   不动点迭代得到的近似解

    % 每次迭代首先提取全局变量
    global mAA;
    global Q R;
    global g_old f_old;
    global DG;
    global mMax;

    if mMax == 0
        x = gval;
    else

        % 计算残差
        fval = gval - x;

        if iter > 0
            dg = gval - g_old;
            df = fval -f_old;

            if mAA < mMax
                DG = [DG, dg];
            else
                DG = [DG(:, 2:mAA), dg];
            end

            mAA = mAA + 1;

        end

        f_old = fval;
        g_old = gval;

        if mAA == 0
            x = gval;
        else

            if mAA == 1
                R(1, 1) = norm(df);
                Q = R(1, 1) \ df;
            else

                % QR delete
                if mAA > mMax
                    [Q, R] = qrdelete(Q, R, 1);
                    mAA = mAA - 1;

                    if size(R, 1) ~= size(R, 2)
                        Q = Q(:, 1:mAA - 1);
                        R = R(1:mAA - 1, :);
                    end

                end

                % QR append
                for j = 1:mAA - 1
                    R(j, mAA) = Q(:, j)' * df;
                    df = df - R(j, mAA) * Q(:, j);
                end

                % 组装 Q R 矩阵
                R(mAA, mAA) = norm(df);
                Q = [Q, R(mAA, mAA) \ df];

            end

            % 求解最小二乘问题, 得到组合参数
            gamma = R \ (Q' * fval);
            % 更新 Anderson 加速近似解.
            x = gval - DG * gamma;
        end

    end

end
