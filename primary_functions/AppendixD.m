% Appendix D Function

% The "Appendix D" function calculates the effectiveness of a control
% surface based on the ratio of the moving component's chord to the chord
% of the surface it is placed on, at the same spanwise station. It
% primarily works only for wings that have a constant chord fraction across
% the dealt-with surface. It *may* work for average chord fraction,
% end-to-end.

%% INPUTS & OUTPUTS

% tau = the non-dimensional effectiveness multiplier of the control surface

% flapRatio = the constant chord fraction of the surface

function tau = AppendixD(flapRatio)
    tau = 1.340933 + (0.00003390316 - 1.340933)/(1 + (flapRatio/0.4437918)^1.331642);
end