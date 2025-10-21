function[SindivR, TindivR, SmeanA_subj, TmeanA_subj] = set_stim_size(setsize, wf, nstim, stdr)

% for standard set
SmeanA_phys = (stdr^2) * pi; % standard mean area (pixel), physical
SmeanA_subj = SmeanA_phys^0.76; % conversion to the perceived size (area), this takes account the perceived differences in stimulus intensity

% % random jitter from a log-uniform distribution between -16 and 19%
% pd = makedist("Loguniform", SmeanA_subj_init * 0.84, SmeanA_subj_init * 1.19); 
% SmeanA_subj = random(pd, 1);

getlast = true;
while getlast
    smallset = (SmeanA_subj * 0.7) + (SmeanA_subj - (SmeanA_subj * 0.7)) * rand(1, (nstim/2)); % up to -30% 
    bigset = SmeanA_subj + ((SmeanA_subj * 1.42) - SmeanA_subj) * rand(1, ((nstim/2) - 1)); % up to 42% (should be always smaller than MaxA)
    SindivA = sort([smallset, bigset]);

    lastindivA = (SmeanA_subj * nstim) - sum(SindivA);

    if (lastindivA > (SmeanA_subj * 0.7)) && (lastindivA < (SmeanA_subj * 1.42)) % falls in -30~42% range
        SindivA = sort([SindivA, lastindivA]);
        getlast = false;
    end
end

SindivD_phys = 2 * sqrt(SindivA.^(1/0.76) ./ pi); % conversion to the physical size (diameter)
SindivR = SindivD_phys ./ 2; % radii

% for test set
% Weber fraction : (St - Ss)/Ss (apparent mean size)
% (−0.244, −0.131, −0.068, −0.034, 0, 0.036, 0.073, 0.150, and 0.323).

TmeanA_subj = (SmeanA_subj * wf) + SmeanA_subj; % target mean size of a test set

getlast = true;
while getlast
    smallset = (TmeanA_subj * 0.7) + (TmeanA_subj - (TmeanA_subj * 0.7)) * rand(1, (setsize/2)); % use the same stim size range as the standard set
    bigset = TmeanA_subj  + ((TmeanA_subj  * 1.42) - TmeanA_subj ) * rand(1, ((setsize/2) - 1)); 
    TindivA = sort([smallset, bigset]);

    lastindivA = (TmeanA_subj  * setsize) - sum(TindivA);
    if (lastindivA > (TmeanA_subj * 0.7)) && (lastindivA < (TmeanA_subj * 1.42))
        TindivA = sort([TindivA, lastindivA]);
        getlast = false;
    end
end

TindivD_phys = 2 * sqrt(TindivA.^(1/0.76) ./ pi); % conversion to the physical size (diameter)
TindivR = TindivD_phys ./ 2; % radii

end