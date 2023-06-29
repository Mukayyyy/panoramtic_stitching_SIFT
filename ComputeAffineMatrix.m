function H = ComputeAffineMatrix( Pt1, Pt2 )
%ComputeAffineMatrix 
%   Computes the transformation matrix that transforms a point from
%   coordinate frame 1 to coordinate frame 2
%Input:Pt1,Pt2
%Output:H(3 * 3 affine transformation matrix)

    N = size(Pt1,1);
    if size(Pt1, 1) ~= size(Pt2, 1),
        error('Dimensions unmatched.');
    elseif N<3
        error('At least 3 points are required.');
    end
    %转换成公共坐标
    P1 = [Pt1';ones(1,N)];
    P2 = [Pt2';ones(1,N)];

    H = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                YOUR CODE HERE:                               %
%        Use MATLAB's "A\b" syntax to solve for H_transpose as discussed       %
%                     above, then convert it to the final H                    %

    H_transpose=P1'\P2';
    H=H_transpose';

%                                END OF YOUR CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sometimes numerical issues cause least-squares to produce a bottom
    % row which is not exactly [0 0 1], which confuses some of the later
    % code. So we'll ensure the bottom row is exactly [0 0 1].
    H(3,:) = [0 0 1];
end
