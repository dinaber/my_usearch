function [GrpCenter,spacer_ind,spacer_rev, len] = my_usearch(Seqs,score0,fnc,gapopen,gapextend,isall,isreverese)
if nargin<5
    gapextend = [];
end

int2nc = 'ACGT' ;

lSeqs = length(Seqs) ;
iSeqs = {} ;
for i=1:lSeqs
    iSeqs{i} = my_nt2int(Seqs{i}) ;
end

iGrpCenter = {} ;
spacer_ind = zeros(lSeqs,1) ;
spacer_rev = zeros(lSeqs,1) ;
tic
t0 = toc ; t1 = t0 ;
delay = [];
for i=1:lSeqs
    if ~mod(i,100)
        fprintf('%4d /%4d, %5.1fsec\n',length(iGrpCenter),i,toc-t1) ;
        t1 = toc ;
    end
    iSeq = iSeqs{i} ;
    score = nan(size(iGrpCenter)) ;
    isrev = nan(size(iGrpCenter)) ;
    for j = 1:length(iGrpCenter)
        [score(j),len, isrev(j)] = my_align_int(iSeq,iGrpCenter{j},fnc,gapopen,gapextend,isreverese) ;        
        if ~isall && score(j)>score0, break; end
    end
    [mx,j] = nanmax(score) ;
    if isnan(mx); error(''); end
    if ~isempty(mx) && mx>score0 % insert into group
        spacer_ind(i) = j ;
        if j==676
            disp('yes')
        end
        spacer_rev(i) = isrev(j) ;
    elseif length(iSeq)<70 %%    % open a new group
        j = length(iGrpCenter) + 1 ;
        iGrpCenter{j} = iSeq ;
        GrpCenter{j} = int2nc(iSeq) ;
        spacer_ind(i) = length(iGrpCenter) ;
    else  %too long - delay
        delay(end+1) = i;
    end
end

for d=1:length(delay)
    if ~mod(i,100)
        fprintf('DELAYED: %4d /%4d, %5.1fsec\n',length(iGrpCenter),d,toc-t1) ;
        t1 = toc ;
    end
    iSeq = iSeqs{delay(d)} ;
    score = nan(size(iGrpCenter)) ;
    isrev = nan(size(iGrpCenter)) ;
    for j = 1:length(iGrpCenter)
        [score(j),len, isrev(j)] = my_align_int(iSeq,iGrpCenter{j},fnc,gapopen,gapextend,isreverese) ;        
        if ~isall && score(j)>score0, break; end
    end
    [mx,j] = nanmax(score) ;
    if isnan(mx); error(''); end
    if ~isempty(mx) && mx>score0 % insert into group
        spacer_ind(delay(d)) = j ;
        spacer_rev(delay(d)) = isrev(j) ;
    else                          % open a new group
        j = length(iGrpCenter) + 1 ;
        iGrpCenter{j} = iSeq ;
        GrpCenter{j}  = int2nc(iSeq) ;
        spacer_ind(delay(d)) = length(iGrpCenter) ;
    end
end
fprintf('Total time: %5.1f\n',toc-t0)


