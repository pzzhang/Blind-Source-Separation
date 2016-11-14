function show_play(y,m,str,H)
   
% myfigure;
for i=1:m 
%    temp=fft(y(i,:));
%    y(i,:)=real(ifft(temp.*(abs(temp)>0.5)));
%    soundsc(y(i,:),16000);
   figure
   plot(y(i,:));
   text = ['the ', num2str(i), 'th signal; ', str];
   title(text);
end  
for i=1:m
    disp(['press any key to hear the sound from ', str]);   pause;    soundsc(y(i,:),16000);
end

% if nargin==4
%     myfigure;
%     for i=1:m
%         for j=1:m
%             subplot(m,m,(i-1)*m+j);       plot(H((i-1)*m+j,:))
%         end
%     end
% end