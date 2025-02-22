% --- Task 1 --- 
% % %prompt user to enter input values
% 
% A = input('value for A = ');
% B = input('value for B = ');
% C = input('value for C = ');
% D = input('value for D = ');
% E = input('value for E = ');
% F = input('value for F = ');
% G = input('value for G = ');
% H = input('value for H = ');
% I = input('value for I = ');
% 
% 
% 
% % Construct matrix seperate columns with space and rows with semicolon
% matrix0 = [A D G; B E H; C F I];
% %disp (matrix0)

%construct a function to find a minor 

function [minor] = find_minor(matrix, element) 
    minor_matrix = matrix; % set varaible for minor matrix as copy of original matrix        
   
    % search matrix for element, assign row and column variables to row and column of element 
    [row, column] = find(matrix == element);  
    minor_matrix(row,:) = []; % remove row of element from matrix (replaces each element in row with empty element)
    minor_matrix(:,column) = []; % remove column of element from matrix (replaces each element in column with empty element)
    % calculate minor of new array using equation for minor
    minor = minor_matrix(1,1)* minor_matrix(2,2) - minor_matrix(1,2)* minor_matrix(2,1);
end


% find determinate using cofactor technique
% create 3 2x2 matrixes with 3 minors and determinates
function [det] = find_determinate(matrix)  % find determinate using elements in top row
    det = 0; % set initial determinate value as 0,
    for I = 1:3 % for each top row element of array 

        % use find_minor function to find the minor of top three elements
        % of matrix, multiply that number by cofactor = -1^(I+1) where I
        % cycles from 2 to 4
        % det =+ matrix(1,I) * find_minor(matrix, matrix(1,I)) * cofactor

        Cofactor = -1^(I+1); % Cofactor variable
        Element = matrix(1,I); % Element variable
        minor = find_minor(matrix, Element); % use find minor function to find minor
        det =+ Element * Cofactor * minor; % multiply and increment determiante

    end
end

determinate0 = find_determinate(matrix0);

% --- Task 2 ---
% calculation of branch current by mesh analysis and Cramer's Rule

V1 = 6; % Set voltage value as 6V

 % Set Value for Resistors in 𝛺
R1 = 330; % Set Value for R1 
R2 = 330; % Set Value for R2
R3 = 330; % Set Value for R3
R4 = 330; % Set Value for R4
R5 = 330; % Set Value for R5
R6 = 330; % Set Value for R6

% create a function which returns the 3 Current values in the 3 branch
% current mesh

% Create new matrix for Mesh matrix
A1 = R1 + R2;
B1 = R2;
C1 = 0;
D1 = -R2;
E1 = R2 + R3 + R4;
F1 = -R4;
G1 = 0;
H1 = R4;
I1 = R4 + R5 + R6;


% Construct matrix seperate columns with space and rows with semicolon
matrix1 = [A1 D1 G1; B1 E1 H1; C1 F1 I1];
disp (matrix1)

% create function to find characteristic determinate 
function[char_det] = find_characteristic_determinate(matrix)
    char_det = 0; % set initial determinate value as 0,
    for I = 1:3 % for each left column element of array 

        % use find_minor function to find the minor of top three elements
        % of matrix, multiply that number by cofactor = -1^(I+1) where I
        % cycles from 2 to 4
        % det =+ matrix(1,I) * find_minor(matrix, matrix(1,I)) * cofactor

        Cofactor = 1^(I+1); % Cofactor variable
        Element = matrix(I,1); % Element variable
        minor = find_minor(matrix, Element); % use find minor function to find minor
        char_det =+ (Element * Cofactor * minor); % multiply and increment determiante
    end
end
% Use function from before to employ Cramers rule


function [I1, I2, I3] = Branch_Current_Calculation(voltage, matrix)
    % Set starting I values so a list can be initialized with them
    I1 = 0;
    I2 = 0;
    I3 = 0;
    Current = [I1,I2,I3]; 

    %char_det = find_characteristic_determiante(matrix);
    char_det = find_determinate(matrix);

    for I = 1:3 % for each top row element of array 

        % use find_minor function to find the minor of top three elements
        % of matrix, multiply that number by cofactor = -1^(I+1) where I
        % cycles from 2 to 4
        % det =+ matrix(1,I) * find_minor(matrix, matrix(1,I)) * cofactor

        cofactor = 1^(I+1); % Cofactor variable
        element = matrix(1,I); % Element variable
        minor = find_minor(matrix, element); % use find minor function to find minor
        Current(I) = (voltage * minor * cofactor) / (char_det); % use equation for I
    end
end
I_values1 = Branch_Current_Calculation(V1, matrix0);
disp(I_values1);
% 
% I_values = Branch_Current_Calculation(V1, matrix1);
% disp (I_values);