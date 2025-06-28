function GA_Line_Registration()
    % LINE_REGISTRATION
    % 
    % This script performs 3D line registration in 4D spherical space using
    % both the Clifford Algebra (GA) method and SVD.
    %
    % Test conditions:
    %   - TranslationGT = [-1 -2 -3]
    %   - Rotation angles: x = -10, y = -20, z = -30 (in degrees)
    %   - For noise-free tests, use noiseStrength = 0 and L = 100000000.
    
    %% INITIAL SETUP
    clc; clear; close all;
    rng(10);
    clifford_signature(6, 0);
    
    % Simulation parameters
    [lineCounts, nPointsOptions, noiseStrength, L] = setParameters();
    
    % Preallocate results structure
    results = initializeResults(length(lineCounts));
    
    % Ground truth transformation parameters
    [TranslationGT, angles] = groundTruthParams();
    T_a = translation_rotor(TranslationGT, L);
    [RotationMatrixGT, MotorGT] = motorCalculation(angles(1), angles(2), angles(3), T_a);
    
    %% MAIN LOOP: Process each configuration
    for idx = 1:length(lineCounts)
        numLines = lineCounts(idx);
        nPoints  = nPointsOptions(idx);
        fprintf('Processing %d lines with %d points each...\n', numLines, nPoints);
        
        % Generate noise-free lines and their start/end points
        [linePoints, noiseFreeStartEndPoints] = generateLines(numLines, nPoints);
        
        % Add noise (if any)
        [noisyLinePoints, ~] = addNoiseToLines(linePoints, nPoints, noiseStrength, numLines);
        
        % Fit lines using TLS to obtain furthest points
        furthest_points = fitLinesToNoisyLines(noisyLinePoints, numLines);
        
       % Form spherical lines and normalize them
        [sphericalLines, sphericalLinesNoise] = computeSphericalLines_fromDirMom( ...
            noiseFreeStartEndPoints, ...   % Nx6 array: [start_x,start_y,start_z, end_x,end_y,end_z]
            furthest_points, ...            % Nx6 array of the TLS‐fitted start/end
            L);                             % curvature parameter

                                                                
        % Transform the noisy spherical lines using the ground truth motor
        sphericalLinesNoiseTransformed = bivector(MotorGT * sphericalLinesNoise * reverse(MotorGT));
        
        % Recover the motor from the spherical lines
        results.MotorRecovered{idx} = processSphericalLines(numLines, sphericalLines, ...
                                                            sphericalLinesNoiseTransformed);
        % Recover the 3D rotation (as a rotor) and translation from the recovered motor
        [results.Rotor3D{idx}, results.tRecovered(idx, :)] = Rotor_and_Translation_from_Motor(L, results.MotorRecovered{idx});
        results.RotationMatrixRecovered{idx} = Rotation_Matrix_Recover(results.Rotor3D{idx});
        results.tRecoveredfrom3DForm(idx, :) = translationFrom3DRotor(furthest_points, noiseFreeStartEndPoints, ...
                                                                      RotationMatrixGT, TranslationGT, ...
                                                                      results.RotationMatrixRecovered{idx});
        % Compute a rigid transformation using SVD for comparison
        transformedFurthestPoints = transformFurthestPoints(furthest_points, RotationMatrixGT, TranslationGT);
        [results.Rsvd{idx}, results.tsvd(idx, :)] = computeRigidTransform(transformedFurthestPoints, noiseFreeStartEndPoints);
    end
    
    %% CALCULATE ERRORS AND PLOT RESULTS
    [convErr, transErr, convErrSVD, transErrSVD] = calculateErrors(RotationMatrixGT, ...
        results.RotationMatrixRecovered, results.Rsvd, TranslationGT, ...
        results.tRecoveredfrom3DForm, results.tsvd, lineCounts);
    plotErrors(lineCounts, nPointsOptions, convErr, transErr, convErrSVD, transErrSVD);
end

%% ======= PARAMETER AND RESULT INITIALIZATION FUNCTIONS =======

function [lineCounts, nPointsOptions, noiseStrength, L] = setParameters()
    lineCounts     = [50, 100, 200, 400, 800, 1600, 3200];
    nPointsOptions = [5, 10, 20, 40, 80, 160, 320];
    noiseStrength  = 0;            % e.g. 0.01 for Gaussian noise, 0 for noise-free
    L              = 100000000;    % Use 100000000 for no noise; lower 
    % values (e.g. 1000 for noisy tests)
end

function results = initializeResults(nConfigs)
    results.MotorRecovered          = cell(1, nConfigs);
    results.Rotor3D                 = cell(1, nConfigs);
    results.tRecovered              = zeros(nConfigs, 3);
    results.tRecoveredfrom3DForm    = zeros(nConfigs, 3);
    results.RotationMatrixRecovered = cell(1, nConfigs);
    results.Rsvd                    = cell(1, nConfigs);
    results.tsvd                    = zeros(nConfigs, 3);
end

function [TranslationGT, angles] = groundTruthParams()
    TranslationGT = [-1 -2 -3];
    angles = [-10, -20, -30]; % Rotation angles [x, y, z] in degrees
end

%% ======= MAPPING AND TRANSFORMATION FUNCTIONS =======

function [pointsStart, pointsEnd, furthestPointsStart, furthestPointsEnd] = mapPointsToSpherical(L, startEndPoints, furthest_points)
    pointsStart       = up_1D(L, startEndPoints(:, 1:3));
    pointsEnd         = up_1D(L, startEndPoints(:, 4:6));
    furthestPointsStart = up_1D(L, furthest_points(:, 1:3));
    furthestPointsEnd   = up_1D(L, furthest_points(:, 4:6));
end

function [sphericalLines, sphericalLinesNoise] = computeSphericalLines_fromDirMom(...
        startEndPoints, furthestPointsSE, L)
  % computeSphericalLines_fromDirMom  Build 4D “1DUp” line bivectors from 3D Plücker data
  %
  % Inputs:
  %   startEndPoints   N×6 [p0_x p0_y p0_z  p1_x p1_y p1_z]  (noise‐free endpoints)
  %   furthestPointsSE N×6 [q0_x q0_y q0_z  q1_x q1_y q1_z]  (TLS‐fitted endpoints)
  %   L                curvature λ
  %
  % Outputs:
  %   sphericalLines      1×N Clifford array of normalized line bivectors (noise‐free)
  %   sphericalLinesNoise 1×N Clifford array of normalized line bivectors (noisy)

  n     = size(startEndPoints, 1);
  delta = 1 / L;  % δ = 1/λ

  % Preallocate 1×n arrays of multivectors
  sphericalLines      = repmat(clifford(0), 1, n);
  sphericalLinesNoise = repmat(clifford(0), 1, n);

  for i = 1:n
    % –– Noise‐free line ––
    p0_0 = startEndPoints(i,1:3);   % first endpoint
    p1_0 = startEndPoints(i,4:6);   % second endpoint

    % 1) Direction d0 = (p1_0 − p0_0) / ‖p1_0 − p0_0‖
    %    i.e. d0 = (p1_0 - p0_0) / norm(p1_0 - p0_0)
    d0 = p1_0 - p0_0;
    d0 = d0 / norm(d0);

    % 2) Plücker moment m0 = p0_0 × d0
    m0 = cross(p0_0, d0);

    % 3) 3D bivector M0 = m0 ∧ I3 = m0(1)e23 + m0(2)e31 + m0(3)e12
    %    but in toolbox e31 = -e13
    Mb0 = m0(1)*e23 - m0(2)*e13 + m0(3)*e12;

    % 4) Form 4D bivector Lb0 = (δ^0 term) + (δ^1 term):
    %      -d0 ∧ e4  =  -d0(1)e14 - d0(2)e24 - d0(3)e34
    %    + (2/λ) M0 = 2*delta*Mb0
    Lb0 = -d0(1)*e14 - d0(2)*e24 - d0(3)*e34 + 2*delta*Mb0;

    % 5) Normalize so (Ln0)^2 = -1, and flip sign for convention:
    Ln0 = -unit(Lb0);
    sphericalLines(i) = Ln0;


    % –– Noisy (TLS‐fitted) line ––
    p0_1 = furthestPointsSE(i,1:3);
    p1_1 = furthestPointsSE(i,4:6);

    % Repeat steps 1–5 with the noisy endpoints
    d1 = p1_1 - p0_1;      d1 = d1 / norm(d1);
    m1 = cross(p0_1, d1);
    Mb1 = m1(1)*e23 - m1(2)*e13 + m1(3)*e12;
    Lb1 = -d1(1)*e14 - d1(2)*e24 - d1(3)*e34 + 2*delta*Mb1;
    Ln1 = -unit(Lb1);
    sphericalLinesNoise(i) = Ln1;
  end
end




%function [sphericalLines, sphericalLinesNoise] = computeSphericalLines(pointsStart, pointsEnd, furthestPointsStart, furthestPointsEnd)
    %sphericalLines      = -unit(wedge(pointsStart, pointsEnd));
    %sphericalLinesNoise = -unit(wedge(furthestPointsStart, furthestPointsEnd));
%end

function transformedPoints = transformFurthestPoints(furthest_points, RGT, TGT)
    transformedPoints = zeros(size(furthest_points));
    for i = 1:size(furthest_points, 1)
        transformedPoints(i, 1:3) = (RGT * furthest_points(i, 1:3).' + TGT').';
        transformedPoints(i, 4:6) = (RGT * furthest_points(i, 4:6).' + TGT').';
    end
end

%% ======= ERROR CALCULATION AND PLOTTING =======

function [convErr, transErr, convErrSVD, transErrSVD] = calculateErrors(RGT, RRecovered, RSVD, TGT, tRecovered, tSVD, lineCounts)
    N = length(lineCounts);
    convErr = zeros(N, 1);
    transErr = zeros(N, 1);
    convErrSVD = zeros(N, 1);
    transErrSVD = zeros(N, 1);
    for i = 1:N
        convErr(i) = norm(eye(3) - RGT * RRecovered{i}', 'fro');
        convErrSVD(i) = norm(eye(3) - RGT * RSVD{i}', 'fro');
        transErr(i) = norm(TGT - tRecovered(i, :));
        transErrSVD(i) = norm(TGT - tSVD(i, :));
    end
end

function plotErrors(lineCounts, nPointsOptions, convErr, transErr, convErrSVD, transErrSVD)
    tickLabels = arrayfun(@(lc, np) sprintf('%d lines / %d pts', lc, np), ...
                          lineCounts, nPointsOptions, 'UniformOutput', false);
    figure;
    subplot(1,2,1);
    plot(1:length(lineCounts), convErr, '-o', 'LineWidth', 2); hold on;
    plot(1:length(lineCounts), convErrSVD, '-^', 'LineWidth', 2); hold off;
    ylabel('Rotation Error (rad)');
    legend('GA', 'SVD', 'Location', 'northeast');
    grid on;
    set(gca, 'FontSize', 16, 'XTick', 1:length(lineCounts), 'XTickLabel', tickLabels, 'XTickLabelRotation', 45);
    
    subplot(1,2,2);
    plot(1:length(lineCounts), transErr, '-o', 'LineWidth', 2); hold on;
    plot(1:length(lineCounts), transErrSVD, '-^', 'LineWidth', 2); hold off;
    ylabel('Translation Error (cm)');
    legend('GA', 'SVD', 'Location', 'northeast');
    grid on;
    set(gca, 'FontSize', 16, 'XTick', 1:length(lineCounts), 'XTickLabel', tickLabels, 'XTickLabelRotation', 45);
end

%% ======= HELPER FUNCTIONS =======

function [linePoints, startEndPoints] = generateLines(numLines, nPoints)
    % GENERATELINES creates numLines 3D lines, each with nPoints.
    linePoints = zeros(numLines, nPoints, 3);
    startEndPoints = zeros(numLines, 6);
    for i = 1:numLines
        startPt = rand(1, 3) * 20 - 10;
        direction = rand(1, 3) * 2 - 1;
        direction = direction / norm(direction);
        t = linspace(0, 10, nPoints);
        points = bsxfun(@plus, startPt, bsxfun(@times, t', direction));
        linePoints(i, :, :) = points;
        startEndPoints(i, :) = [startPt, points(end, :)];
    end
end

function [noisyLinePoints, noisyStartEndPoints] = addNoiseToLines(linePoints, nPoints, noiseStrength, numLines)
    % ADDNOISETOLINES adds isotropic Gaussian noise to line points.
    noisyLinePoints = linePoints;
    noisyStartEndPoints = zeros(numLines, 6);
    noise = noiseStrength * randn(numLines, nPoints, 3);
    noisyLinePoints = noisyLinePoints + noise;
    for i = 1:numLines
        noisyStartEndPoints(i, 1:3) = noisyLinePoints(i, 1, :);
        noisyStartEndPoints(i, 4:6) = noisyLinePoints(i, end, :);
    end
end

function furthest_points = fitLinesToNoisyLines(noisyLinePoints, numLines)
    % FITLINESTONOISYLINES fits each noisy line using TLS and returns furthest points.
    furthest_points = zeros(numLines, 6);
    for i = 1:numLines
        points = squeeze(noisyLinePoints(i, :, :));
        [direction, point, ~, ~] = fitLineTLS(points);
        centeredPoints = bsxfun(@minus, points, point);
        projections = centeredPoints * direction;
        minProj = min(projections);
        maxProj = max(projections);
        linePoint1 = point + minProj * direction';
        linePoint2 = point + maxProj * direction';
        furthest_points(i, :) = [linePoint1, linePoint2];
    end
end

function [direction, point, furthestPoint1, furthestPoint2] = fitLineTLS(points)
    % FITLINETLS fits a line via Total Least Squares.
    meanPoint = mean(points, 1);
    centeredPoints = points - meanPoint;
    [~, ~, V] = svd(centeredPoints, 'econ');
    direction = V(:, 1);
    point = meanPoint;
    projections = centeredPoints * direction;
    [minProj, ~] = min(projections);
    [maxProj, ~] = max(projections);
    furthestPoint1 = meanPoint + minProj * direction';
    furthestPoint2 = meanPoint + maxProj * direction';
end

function data_spherical = up_1D(L, data)
    % UP_1D maps 3D Euclidean points to 4D spherical space.
    data = double(data);
    if size(data,2) == 3
        data = data'; % Ensure 3xN format
    end
    E = [e1, e2, e3];
    data_clifford = E * data;
    first_term = 2*L ./ (L^2 + data_clifford.^2) .* data_clifford;
    second_term = (L^2 - data_clifford.^2) ./ (L^2 + data_clifford.^2) * e4;
    data_spherical = first_term + second_term;
    data_spherical = vector(data_spherical);
end

function translation_euclidean = down_1D(L, data_spherical)
    % DOWN_1D maps 4D spherical points back to 3D Euclidean space.
    w4 = scalar_product(data_spherical, e4);
    first_term_down = L ./ (1 + w4);
    second_term_down = scalar_product(data_spherical, e1)*e1 + ...
                       scalar_product(data_spherical, e2)*e2 + ...
                       scalar_product(data_spherical, e3)*e3;
    translation_euclidean = cell2mat(coefficients(vector(first_term_down .* second_term_down)));
end

function [OriginalRotationMatrix, OriginalMotor] = motorCalculation(angleX_degrees, angleY_degrees, angleZ_degrees, translationRotor)
    % MOTORCALCULATION computes the ground truth motor.
    OriginalRotationMatrix = create3DRotationMatrix(angleX_degrees, angleY_degrees, angleZ_degrees);
    Rotor = rotation_matrix_to_rotor(OriginalRotationMatrix);
    OriginalMotor = translationRotor * Rotor;
end

function R = create3DRotationMatrix(angleX_degrees, angleY_degrees, angleZ_degrees)
    % CREATE3DROTATIONMATRIX returns a 3D rotation matrix from angles (in degrees).
    angleX = deg2rad(angleX_degrees);
    angleY = deg2rad(angleY_degrees);
    angleZ = deg2rad(angleZ_degrees);
    Rx = [1, 0, 0; 0, cos(angleX), -sin(angleX); 0, sin(angleX), cos(angleX)];
    Ry = [cos(angleY), 0, sin(angleY); 0, 1, 0; -sin(angleY), 0, cos(angleY)];
    Rz = [cos(angleZ), -sin(angleZ), 0; sin(angleZ), cos(angleZ), 0; 0, 0, 1];
    R = Rz * Ry * Rx;
end

function Rotor = rotation_matrix_to_rotor(Rotation_Matrix)
    % ROTATION_MATRIX_TO_ROTOR converts a 3D rotation matrix to a Clifford rotor.
    Rotation_Matrix = double(Rotation_Matrix);
    Frame = [Rotation_Matrix(1,1)*e1, Rotation_Matrix(1,2)*e2, Rotation_Matrix(1,3)*e3;
             Rotation_Matrix(2,1)*e1, Rotation_Matrix(2,2)*e2, Rotation_Matrix(2,3)*e3;
             Rotation_Matrix(3,1)*e1, Rotation_Matrix(3,2)*e2, Rotation_Matrix(3,3)*e3];
    E = [e1, e2, e3];
    Rotor = 1 + sum(E * Frame);
    Rotor = Rotor / abs(Rotor);
end

function T_a = translation_rotor(t, L)
    % TRANSLATION_ROTOR constructs a translation rotor from a 1x3 translation vector.
    t = double(t);
    if size(t,1) ~= 1, t = t'; end
    E = [e1, e2, e3];
    t_clifford = E * t';
    T_a = (L + t_clifford*e4) ./ (sqrt(L^2 + t_clifford.^2));
end

function MotorRecovered = processSphericalLines(noOfLines, sphericalLines, transformedSphericalLines)
    % PROCESSSPHERICALLINES recovers the motor from spherical lines.
    M = reshape(cell2mat(coefficients([sphericalLines(:)])), noOfLines, []).';
    M_dash = reshape(cell2mat(coefficients([transformedSphericalLines(:)])), noOfLines, []).';
    euclidean_basis = [e1, e2, e3, e4, e5, e6];
    v_i = euclidean_basis * M;
    v_i_dash = euclidean_basis * M_dash;
    [F, G] = computeFandGMatrices(v_i, v_i_dash);
    [G_basis_biv, F_L_i_upstairs] = formReciprocalsInF(F, G);
    X = (1/8) * (sum(G_basis_biv .* F_L_i_upstairs) + 2);
    X_reverse = reverse(X);
    Y = X * X_reverse;
    Y_coeff = cell2mat(coefficients(Y));
    u = cell2mat(coefficients(scalar(Y)));  % Scalar part
    if length(Y_coeff) == 8
        v = Y_coeff(8);
    else
        exprStr = evalc('disp(Y)');
        pattern = '([-+\s]*[0-9.]+)\s*e1234';
        matches = regexp(exprStr, pattern, 'tokens');
        if ~isempty(matches)
            quadrivectorCoeffStr = strtrim(matches{1}{1});
            quadrivectorCoeffStr = regexprep(quadrivectorCoeffStr, '[^\d.-]', '');
            v = str2double(quadrivectorCoeffStr);
        end
    end
    w_positive = sqrt(u + v);
    w_negative = -sqrt(u + v);
    w = [w_positive, w_negative];
    num_solutions = 2;
    b_all = zeros(num_solutions*2, 1);
    a_all = zeros(num_solutions*2, 1);
    for i = 1:length(w)
        discriminant = w(i)^2 - 2*v;
        if discriminant >= 0
            sqrt_discriminant = sqrt(discriminant);
        else
            error('Discriminant is negative. No real solutions exist for beta.');
        end
        b1 = (w(i) + sqrt_discriminant) / 2;
        b2 = (w(i) - sqrt_discriminant) / 2;
        a1 = w(i) - b1;
        a2 = w(i) - b2;
        index = (i-1)*2 + 1;
        b_all(index) = b1;
        a_all(index) = a1;
        b_all(index+1) = b2;
        a_all(index+1) = a2;
    end
    I = e1234;
    MotorRecoveredALL = cell(length(a_all), 1);
    for i = 1:length(a_all)
        alpha = a_all(i);
        beta  = b_all(i);
        MotorRecoveredALL{i} = (1/u) * (alpha - beta * I) * X;
    end
    % for i = 1:length(MotorRecoveredALL)
        % fprintf('Motor %d: %s\n', i, char(MotorRecoveredALL{i}));
    % end
    MotorRecovered = MotorRecoveredALL{2};
end

function [F, G] = computeFandGMatrices(v_i, v_i_dash)
    % COMPUTEFANDGMATRICES computes the F and G matrices.
    F = zeros(6,6);
    G = zeros(6,6);
    e_basis = {e1, e2, e3, e4, e5, e6};
    for i = 1:6
        for j = 1:6
            Fija = bsxfun(@scalar_product, v_i(1,:), e_basis{i});
            Fijb = bsxfun(@scalar_product, v_i_dash(1,:), e_basis{j});
            Fija_coef = cell2mat(coefficients(Fija));
            Fijb_coef = cell2mat(coefficients(Fijb));
            maxSizeF = max(length(Fija_coef), length(Fijb_coef));
            vector_ija_F = [Fija_coef, zeros(1, maxSizeF - length(Fija_coef))];
            vector_ijb_F = [Fijb_coef, zeros(1, maxSizeF - length(Fijb_coef))];
            F(i,j) = sum(vector_ija_F .* vector_ijb_F);
            
            Gija = bsxfun(@scalar_product, v_i_dash(1,:), e_basis{i});
            Gijb = bsxfun(@scalar_product, v_i_dash(1,:), e_basis{j});
            Gija_coef = cell2mat(coefficients(Gija));
            Gijb_coef = cell2mat(coefficients(Gijb));
            maxSizeG = max(length(Gija_coef), length(Gijb_coef));
            vector_ija_G = [Gija_coef, zeros(1, maxSizeG - length(Gija_coef))];
            vector_ijb_G = [Gijb_coef, zeros(1, maxSizeG - length(Gijb_coef))];
            G(i,j) = sum(vector_ija_G .* vector_ijb_G);
        end
    end
end

function [G_basis_biv, F_L_i_upstairs] = formReciprocalsInF(F, G)
    % FORMRECIPROCALSINF computes the reciprocals in the F matrix.
    e_basis = [e1, e2, e3, e4, e5, e6];
    bivector_basis = [e12, e13, e14, e23, e24, e34];
    F_basis = e_basis * F;
    G_basis = e_basis * G;
    vol6D = wedge(F_basis(1), F_basis(2), F_basis(3), F_basis(4), F_basis(5), F_basis(6));
    vol6D_inv_sq = 1 / abs(vol6D)^2;
    F_i_upstairs = clifford(zeros(1,6));
    for i = 1:6
        indices = [1:(i-1), (i+1):6];
        argsForWedge = arrayfun(@(idx) F_basis(idx), indices, 'UniformOutput', false);
        wedgeResult = wedge(argsForWedge{:});
        F_i_upstairs(i) = ((-1)^(i+1)) * vol6D_inv_sq * wedgeResult * vol6D;
    end
    M = reshape(cell2mat(coefficients(F_i_upstairs)), 6, [])';
    M2 = reshape(cell2mat(coefficients(G_basis)), 6, [])';
    F_L_i_upstairs = bivector_basis * M;
    G_basis_biv = bivector_basis * M2;
end

function R = Rotation_Matrix_Recover(New_Rotor)
    % ROTATION_MATRIX_RECOVER recovers the 3D rotation matrix from a rotor.
    F = zeros(3, 3);
    for i = 1:3
        fi_F = vector(New_Rotor * eval(sprintf('e%d', i)) * reverse(New_Rotor));
        coef_fi_F = cell2mat(coefficients(fi_F));
        if length(coef_fi_F) < 3
            coef_fi_F = [coef_fi_F, zeros(1, 3 - length(coef_fi_F))];
        else
            coef_fi_F = coef_fi_F(1:3);
        end
        F(:, i) = coef_fi_F(:);
    end
    R = F;
end

function [Rotor, Translation] = Rotor_and_Translation_from_Motor(L, Motor_4D)
    % ROTOR_AND_TRANSLATION_FROM_MOTOR extracts the 3D rotor and translation from a 4D motor.
    zero_vector = [0,0,0];
    zero_vector_spherical = up_1D(L, zero_vector);
    t_in_spherical = Motor_4D * zero_vector_spherical * reverse(Motor_4D);
    t_in_spherical = vector(t_in_spherical);
    Translation = down_1D(L, t_in_spherical);
    if numel(Translation) < 3
        Translation = [Translation(:); zeros(3 - numel(Translation), 1)];
    end
    T_a = translation_rotor(Translation, L);
    Rotor = reverse(T_a) * Motor_4D;
end

function tRecovered = translationFrom3DRotor(pointsToBeTransformed, originalPoints, GroundTruthRotationMatrix, GroundTruthTranslation, RotationMatrixFrom3DRotor)
    % TRANSLATIONFROM3DROTOR computes the translation vector from the 3D rotor.
    transformedStart = (GroundTruthRotationMatrix * originalPoints(:, 1:3).' + GroundTruthTranslation').';
    transformedEnd = (GroundTruthRotationMatrix * originalPoints(:, 4:6).' + GroundTruthTranslation').';
    originalPointsALL = [originalPoints(:, 1:3); originalPoints(:, 4:6)];
    transformedPointsALL = [transformedStart; transformedEnd];
    originalCentroid = mean(originalPointsALL, 1);
    transformedCentroid = mean(transformedPointsALL, 1);
    tRecovered = transformedCentroid' - RotationMatrixFrom3DRotor * originalCentroid';
end

function [R, t] = computeRigidTransform(transformedPoints, originalPoints)
    % COMPUTERIGIDTRANSFORM computes the best-fit rigid transformation using SVD.
    originalPointsALL = [originalPoints(:, 1:3); originalPoints(:, 4:6)];
    transformedPointsALL = [transformedPoints(:, 1:3); transformedPoints(:, 4:6)];
    normOriginalPoints = bsxfun(@minus, originalPointsALL, mean(originalPointsALL, 1));
    normTransformedPoints = bsxfun(@minus, transformedPointsALL, mean(transformedPointsALL, 1));
    H = normOriginalPoints.' * normTransformedPoints;
    [U, ~, V] = svd(H);
    R = V * U';
    if det(R) < 0
        V(:, end) = -V(:, end);
        R = V * U';
    end
    t = mean(transformedPointsALL, 1)' - R * mean(originalPointsALL, 1)';
end
