%statistical
A=reshape(Ace_star,31,1440);
%for t=1:24
    for n=1:24
        B(n,:,:)=A(:,1+(n-1)*60:60+(n-1)*60);
        subplot(4,6,n)
        histogram(B(n,:,:));
%         M(n)=mean(B(n,:,:));
%         Med(n)=medium(B(n,:,:));
    end
    
    %histogram(B(1,:,:))
%end