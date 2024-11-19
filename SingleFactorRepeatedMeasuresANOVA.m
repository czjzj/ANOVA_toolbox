function [ANOVA_F, ANOVA_p, ANOVA_Table] = SingleFactorRepeatedMeasuresANOVA(DependentMatrix, OutputName)

% Single Factor Repeated Measures ANCOVA Analysis
% Input:
%   DependentMatrix - A matrix of dependent variable, where each column is a measurement and each row is a subject.
%   OutputName - Name of the output file to save the ANOVA table (optional).
% Output:
%   ANCOVA_F - F values for measure.
%   ANCOVA_p - p values for measure.
%   ANCOVA_Table - A table containing the ANOVA results.
%___________________________________________________________________________
% Written by CUI Zhen-Jiang, 
% JingShi Hotel 9418, State Key Laboratory of Cognitive Neuroscience and Learning, 
% Beijing Normal University, China, 100875, November 17, 2024
% Refer to the second edition of Shu Hua's Multi Factor Experimental Design in Psychology and Education, page 36

    % Default for OutputName
    if nargin < 2
        OutputName = 'SGRBD_ANOVA_Results.txt';
    end

    % AS table
    Sum_Across_Measures = sum(DependentMatrix, 1);
    Sum_Across_Subjects = sum(DependentMatrix, 2);
    DependentVector = DependentMatrix(:);
    numMeasures = size(DependentMatrix, 2);
    numSubjects = size(DependentMatrix, 1);

    % Calculation of various basic quantities
    Yij = sum(DependentVector);
    Y = Yij ^ 2/(numMeasures * numSubjects);
    AS = sum(DependentVector .^ 2);
    A = sum((Sum_Across_Measures .^ 2)/numSubjects);
    S = sum((Sum_Across_Subjects .^ 2)/numMeasures);

    % Decomposition and calculation of sum of squares
    % (1) Decomposition mode
    % SS_total = SS_between + SS_within 
    %          = SS_A(between) + (SS_group + SS_residual)
    % (2) Calculation of sum of squares
    
    SS_total = AS - Y;
    df_total = numMeasures * numSubjects - 1;
    
    SS_A = A - Y;
    df_A = numMeasures - 1;
    
    SS_within = SS_total - SS_A;
    df_within = numMeasures * (numSubjects - 1);

    SS_group = S - Y;
    df_group = numSubjects - 1;

    SS_residual = SS_total - SS_A - SS_group;
    df_residual = (numMeasures - 1)*(numSubjects - 1);

    MS_A = SS_A/df_A;
    MS_group = SS_group/df_group;
    MS_residual = SS_residual/df_residual;
    F_A = MS_A/MS_residual;
    ANOVA_F = F_A;
    p_A = 1-fcdf(F_A, df_A, df_residual);
    ANOVA_p = p_A;
    F_group = MS_group/MS_residual;
    p_group = 1-fcdf(F_group, df_group, df_residual);
    
    % Prepare results for output
    Source = {'Measure', 'Within', 'Group', 'Residual', 'Total'};
    SS = [SS_A, SS_within, SS_group, SS_residual, SS_total];
    df = [df_A, df_within, df_group, df_residual, df_total];
    MS = [MS_A, NaN, MS_group, MS_residual, NaN];
    F = [F_A, NaN, F_group, NaN, NaN];
    p = [p_A, NaN, p_group, NaN, NaN];

    % Replace NaN with empty string
    MS = num2cell(MS); F = num2cell(F); p = num2cell(p);
    MS(cellfun(@(x) isnan(x), MS)) = {''};
    F(cellfun(@(x) isnan(x), F)) = {''};
    p(cellfun(@(x) isnan(x), p)) = {''};
    % Create a table for output
    ANOVA_Table = table(Source', SS', df', MS', F', p', ...
                         'VariableNames', {'Source', 'SS', 'df', 'MS', 'F_value', 'p_value'});

    % Print table to command window
    fprintf('%-15s %-15s %-15s %-15s %-15s %-15s\n', ...
        'Source', 'SS', 'df', 'MSe', 'F-value', 'p-value');
    for i = 1:length(Source)
        
        % Format p-value: If p is too small, display as 0.000
        if isnumeric(p{i}) && ~isempty(p{i}) && p{i} < 1e-4
            p_str = '0.0000';
        else
            p_str = p{i};
        end
        % Print each row with formatted p-value
        fprintf('%-15s %-15.4f %-15d %-15s %-15s %-15s\n', ...
        Source{i}, SS(i), df(i), MS{i}, F{i}, p_str);
    end

    % Open file for writing
    fid = fopen(OutputName, 'w');

    % Write headers
    fprintf(fid, '%-15s %-15s %-15s %-15s %-15s %-15s\n', ...
        'Source', 'SS', 'df', 'MSe', 'F-value', 'p-value');

    % Write rows to file
    for i = 1:length(Source)
        fprintf(fid, '%-15s %-15.4f %-15d %-15s %-15s %-15s\n', ...
            Source{i}, SS(i), df(i), MS{i}, F{i}, p{i});
    end

    % Close file
    fclose(fid);
    disp(['SGRBD_ANOVA table saved to ', OutputName]);
end