function [nonzero_ind, nonzero_X]=select_spikes(X,sparsity,sparse_choice)

% select the highest spikes and nail down their position. 
global q;

nonzero_ind=[]; nonzero_X=[];
aa=zeros(4,q);
aa(1,:)=X(1:q)';   aa(2,:)=X(q+1:2*q)';  aa(3,:)=X(2*q+1:3*q)';  aa(4,:)=X(3*q+1:4*q)';

if sparse_choice==1;  
	%%% select highest spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	for nsort=1:floor(sparsity/4) 
	% only select the first sparsity/4 most significants and make sure that numbers next to them are zero
        for ai=1:4
            [bb,ind]=sort(abs(aa(ai,:)));  ind=ind(end+1-nsort);  
            if ind==1; aa(ai,2)=0; 
            elseif ind==q; aa(ai,q-1)=0;
            else aa(ai,[ind-1,ind+1])=0;  
            end        
        end
	end % end of nsort iteration
	
	for ai=1:4
        [bb,ind]=sort(abs(aa(ai,:)));  aa(ai,ind(1:end-floor(sparsity/4)))=0;
        nonzero_ind=[nonzero_ind, (ai-1)*q + (ind(end-floor(sparsity/4)+1:end))];
        nonzero_X=[nonzero_X, aa(ai,ind(end-floor(sparsity/4)+1:end))];
	end    
end


if sparse_choice==2
    for ai=1:4
        [bb,ind]=sort(abs(aa(ai,:))); 
        highind=ind(end-floor(sparsity/4)+1:end);  % the highest sparsity/4 spikes
        highind=setxor(union([highind,highind-1,highind+1],[0,q+1]),[0,q+1]); % select high spikes and 2 neighbors around it
        nonzero_ind=[nonzero_ind, (ai-1)*q + highind];
        nonzero_X=[nonzero_X, aa(ai,highind)];
    end
end