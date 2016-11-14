function draw_weights(X,draw_style,legend_str)
global q;

a11=X(1:q)'; a12=X(q+1:2*q)'; a21=X(2*q+1:3*q)'; a22=X(3*q+1:4*q)';

subplot(2,2,1); hold on; plot(a11,draw_style); 
legend(legend_str); title('a^{11}');  
subplot(2,2,2); hold on; plot(a12,draw_style); 
legend(legend_str); title('a^{12}');  
subplot(2,2,3); hold on; plot(a21,draw_style); 
legend(legend_str); title('a^{21}');  
subplot(2,2,4); hold on; plot(a22,draw_style); 
legend(legend_str); title('a^{22}');  
