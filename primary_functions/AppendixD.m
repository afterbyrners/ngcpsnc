function tau = AppendixD(flapRatio) % Effectivness calculation
    tau = 1.340933 + (0.00003390316 - 1.340933)/(1 + (flapRatio/0.4437918)^1.331642); % from online curve fit, accurate except at really low ratios like <0.05
end