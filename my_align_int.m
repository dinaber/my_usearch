function [score,lens,isrev] = my_align_int(intSeq1,intSeq2,fnc,gapopen,gapextend,isreverse)
if nargin<5
    gapextend = [] ;
end
if nargin<6
    isreverse = false ;
end
persistent scoringMatrix
if isempty(scoringMatrix)
    scoringMatrix = [...
        5 -4 -4 -4  1 -4 -4  1 -4  1 -4 -1 -1 -1 -2;...
        -4  5 -4 -4 -4  1 -4  1  1 -4 -1 -4 -1 -1 -2;...
        -4 -4  5 -4  1 -4  1 -4  1 -4 -1 -1 -4 -1 -2;...
        -4 -4 -4  5 -4  1  1 -4 -4  1 -1 -1 -1 -4 -2;...
        1 -4  1 -4 -1 -4 -2 -2 -2 -2 -3 -1 -3 -1 -1;...
        -4  1 -4  1 -4 -1 -2 -2 -2 -2 -1 -3 -1 -3 -1;...
        -4 -4  1  1 -2 -2 -1 -4 -2 -2 -1 -1 -3 -3 -1;...
        1  1 -4 -4 -2 -2 -4 -1 -2 -2 -3 -3 -1 -1 -1;...
        -4  1  1 -4 -2 -2 -2 -2 -1 -4 -1 -3 -3 -1 -1;...
        1 -4 -4  1 -2 -2 -2 -2 -4 -1 -3 -1 -1 -3 -1;...
        -4 -1 -1 -1 -3 -1 -1 -3 -1 -3 -1 -2 -2 -2 -1;...
        -1 -4 -1 -1 -1 -3 -1 -3 -3 -1 -2 -1 -2 -2 -1;...
        -1 -1 -4 -1 -3 -1 -3 -1 -3 -1 -2 -2 -1 -2 -1;...
        -1 -1 -1 -4 -1 -3 -3 -1 -1 -3 -2 -2 -2 -1 -1;...
        -2 -2 -2 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;...
        ];
end
scoringMatrixInfo_Scale = 0.277316;

scale = 1 ;
scale = scale * scoringMatrixInfo_Scale ;

%         switch fnc
%             case 'swalign'
%             case 'nwalign'
% later
score = zeros(isreverse+1,1) ;
lens  = zeros(isreverse+1,1) ;
rev_d = uint8([4 3 2 1]);

for irev = 1:isreverse+1    
    if irev==2
        intSeq1p = intSeq1(end:-1:1) ;
        intSeq1p = rev_d(intSeq1p) ;
    else
        intSeq1p = intSeq1 ;
    end
    if isempty(gapextend)
        [score(irev), pth2, pth1] = ...
            bioinfoprivate.simplegapmex(intSeq2, intSeq1p, gapopen, scoringMatrix, 2);
    else
        [score(irev), pth2, pth1] = ...
            bioinfoprivate.affinegapmex(intseq2, intseq1p, gapopen, gapextend, scoringMatrix, 2);
    end    
    i = pth1>0 & pth2>0 ;
    lens(irev) = sum(intSeq1(pth1(i)) == intSeq2(pth2(i))) ;
end
[score,imx] = max(score) ;
isrev = imx==2 ;
lens = lens(imx) ;

score = scale * score ;

end