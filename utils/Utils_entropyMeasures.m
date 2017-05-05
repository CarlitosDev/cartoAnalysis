function [entropyA, entropyB, mutualInfo] =  Utils_entropyMeasures( signalA, signalB)


%% Import JavaMI (http://www.cs.man.ac.uk/~pococka4/JavaMI.html)

mdFolder   = 'D:\Work\MATLABPromotionsModel\External Code\Mutual Information\JavaMI';
mdJarFile  = 'JavaMI.jar';

% Load the metadata java library
mdJarPath   = fullfile(mdFolder, mdJarFile);

jpa = javaclasspath('-all');

if ~any(strcmp(mdJarPath, jpa))
    oldWarn = warning('off', 'MATLAB:Java:DuplicateClass');
    javaaddpath(mdJarPath);
    warning(oldWarn);
end

import JavaMI.*


%% Remove any non numerical and calculate measures

idxNull      = isnan(signalA) | isnan(signalB);
inputSignalA = signalA(~idxNull);
inputSignalB = signalB(~idxNull);

entropyA   = JavaMI.Entropy.calculateEntropy(inputSignalA);
entropyB   = JavaMI.Entropy.calculateEntropy(inputSignalB);
mutualInfo = JavaMI.MutualInformation.calculateMutualInformation(inputSignalA, inputSignalB);

fprintf('\nH(x) signal A is %3.2f\n', entropyA);
fprintf('H(x) signal B is %3.2f\n', entropyB);
fprintf('I(x,y) %3.2f\n', mutualInfo);
