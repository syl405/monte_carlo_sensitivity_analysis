clear

%% Define Global Variables
MEASUREMENT_OFFSET_ERROR = 0;
MEASUREMENT_STD = 0.00;

%% Parse in raw data
raw_filename = 'C:/Users/shien/Documents/NVBOTS/piecewise_constant_compensation/error_data/printer_accuracy_arjun.xlsx';
[~,~,headers] = xlsread(raw_filename, 'Z run charts', 'A1:D1','basic'); 
headers = replace(headers,' ','_'); %strip whitespace
headers = replace(headers,'%','Pct'); %strip special characters
[~,~,R] = xlsread(raw_filename, 'Z run charts', 'A2:D176','basic');
R = cell2table(R,'VariableNames',headers);
R.Measured_Length = R.Nominal_Length + R.Raw_Deviation;

% correct for quantization error
D = table(changem(R.Nominal_Length,[30-20.1 60-20.1 99.9-20.1 180-20.1 240-20.1],[10 40 80 160 220]),'VariableNames',{'Nominal_Length'});
D.Measured_Length = R.Measured_Length;
D.Absolute_Error = D.Measured_Length - D.Nominal_Length;

%% Calculate Summary Stats for Measured Sample
% S = table(unique(D.Nominal_Length+20.1),'VariableNames',{'Nominal_Length'},'RowNames',cellstr(num2str(unique(D.Nominal_Length))));
% nom_length_group = findgroups(D.Nominal_Length);
% S.Mean_Measured_Height = splitapply(@mean,D.Measured_Height,nom_length_group);
S = table([20.1;30;60;99.9;180;240],'VariableNames',{'Target_Height'});
S.Mean_Measured_Height(S.Target_Height==20.1)= 20.1;
S.Std_Measured_Height(S.Target_Height==20.1)= 0;
for i = 2:numel(S.Target_Height)
    S.Mean_Measured_Height(i) = 20.1 + mean(D.Measured_Length(D.Nominal_Length+20.1==S.Target_Height(i)));
    S.Std_Measured_Height(i) = std(D.Measured_Length(D.Nominal_Length+20.1==S.Target_Height(i)));
end

% % Calculate mean absolute error
% S.Mean_Abs_Err = S.Mean_Measured_Height-S.Target_Height;
% 
% % Calculate total error in each block
% S.Error_This_Block(1) = 0;
% S.Num_Layers_This_Block(1) = 20.1/0.3;
% for i = 2:size(S,1)
%     S.Error_This_Block(i) = S.Mean_Abs_Err(i)-S.Mean_Abs_Err(i-1);
%     S.Num_Layers_This_Block(i) = (S.Target_Height(i)-S.Target_Height(i-1))/0.3;
% end
% S.Error_Per_Layer_This_Block = S.Error_This_Block./S.Num_Layers_This_Block;

%% Monte Carlo data generation

%true population data
rng(12345) % set random seed
TRUE_RAND = cell(size(S,1),2); % columns: target_height, measured_height_vect
TRUE_RAND{1,1} = 20.1;
TRUE_RAND{1,2} = NaN;
for i = 2:size(S,1)
    TRUE_RAND{i,1} = S.Target_Height(i);
    TRUE_RAND{i,2} = S.Mean_Measured_Height(i) - S.Target_Height(i) + randn(100,1)*S.Std_Measured_Height(i);
end

%imperfectly measured data
MEASURED_RAND = TRUE_RAND;
for i = 2:size(S,1)
    MEASURED_RAND{i,2} = TRUE_RAND{i,2} + MEASUREMENT_STD * randn(size(TRUE_RAND{i,2},1),1) + MEASUREMENT_OFFSET_ERROR;
end

%% Fit confounded error models for different n
num_levels = size(MEASURED_RAND,1);
poly3_models = cell(10,1);
poly2_models = cell(10,1);
lin_models = cell(10,1);
piecewise_models = cell(10,1);
for n = 1:10 % 100 different sample sizes
       
    % piecewise error model
    breaks = nan(num_levels,1);
    coefs = nan(num_levels-1,2);
    means = nan(num_levels,1);
    breaks(1) = 0;
    means(1) = 0;
    for i = 2:num_levels
        breaks(i) = MEASURED_RAND{i,1};
        means(i) = mean(MEASURED_RAND{i,2}(1:n));
        %coefs(i-1,:) = [(means(i)-means(i-1))/(breaks(i)-breaks(i-1)), means(i-1)-(means(i)-means(i-1))/(breaks(i)-breaks(i-1))*breaks(i-1)];
        coefs(i-1,:) = [(means(i)-means(i-1))/(breaks(i)-breaks(i-1)), means(i-1)];
    end
    piecewise_models{n} = mkpp(breaks,coefs);
        
    % cubic error model
    cur_data = [];
    cur_targets = [];
    for i = 2:num_levels
        cur_targets = [cur_targets; MEASURED_RAND{i,1}*ones(n,1)];
        cur_data = [cur_data; MEASURED_RAND{i,2}(1:n)];
    end
    poly3_models{n} = fitlm(cur_targets,cur_data,'poly3','Intercept',false);
    
    % quadratic error model
    cur_data = [];
    cur_targets = [];
    for i = 2:num_levels
        cur_targets = [cur_targets; MEASURED_RAND{i,1}*ones(n,1)];
        cur_data = [cur_data; MEASURED_RAND{i,2}(1:n)];
    end
    poly2_models{n} = fitlm(cur_targets,cur_data,'quadratic','Intercept',false);
    
    % linear error model
    cur_data = [];
    cur_targets = [];
    for i = 2:num_levels
        cur_targets = [cur_targets; MEASURED_RAND{i,1}*ones(n,1)];
        cur_data = [cur_data; MEASURED_RAND{i,2}(1:n)];
    end
    lin_models{n} = fitlm(cur_targets,cur_data,'linear','Intercept',false);
end

% clear temporary variables
clear cur_data cur_targets i n breaks coefs

%% Fit perfect error model (cubic, n=20, with non-confounded data)
cur_data = [];
cur_targets = [];
for i = 2:num_levels
    cur_targets = [cur_targets; TRUE_RAND{i,1}*ones(size(TRUE_RAND{i,2},1),1)];
    cur_data = [cur_data; TRUE_RAND{i,2}];
end
true_model = fitlm(cur_targets,cur_data,'poly3','Intercept',false);

%% Calculate prediction errors for confounded piecewise and poly models
probing_points = linspace(20.1,240,100)';
true_prediction = predict(true_model,probing_points);
pred_error = nan(10,100,4); %dims: n_used,probing_point(1-pw,2-poly2,3-poly=3),model_type; value: prediction_error
for n = 1:10
    pred_error(n,:,1) = ppval(piecewise_models{n},probing_points)-true_prediction;
    pred_error(n,:,2) = predict(lin_models{n},probing_points)-true_prediction;
    pred_error(n,:,3) = predict(poly2_models{n},probing_points)-true_prediction;
    pred_error(n,:,4) = predict(poly3_models{n},probing_points)-true_prediction;
end

%% Plot predictions errors
figure
subplot(2,2,1)
surf(pred_error(:,:,1))
title('Piecewise Model')
xlabel('z')
ylabel('n')
axis([0 100 0 10 -0.4 0.4])
subplot(2,2,2)
surf(pred_error(:,:,2))
xlabel('z')
ylabel('n')
axis([0 100 0 10 -0.4 0.4])
title('Linear Model')
subplot(2,2,3)
surf(pred_error(:,:,3))
xlabel('z')
ylabel('n')
axis([0 100 0 10 -0.4 0.4])
title('Quadratic Model')
subplot(2,2,4)
surf(pred_error(:,:,4))
xlabel('z')
ylabel('n')
axis([0 100 0 10 -0.4 0.4])
title('Cubic Model')

