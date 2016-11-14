function mixed_x=loadmix(mix_choice)

switch mix_choice
    case 1
        load data/lee1s1;  load data/lee1s2;  mixed_x=[lee1s1-mean(lee1s1);lee1s2-mean(lee1s2)];  

    case 2
        
        load data/lee2s1;  load data/lee2s2;  mixed_x=[lee2s1-mean(lee2s1);lee2s2-mean(lee2s2)];

    case 3
        
        mx1=wavread('data/Amixm22f36')'; 
        mx2=wavread('data/Bmixm22f36')';
        mixed_x=[mx1-mean(mx1);mx2-mean(mx2)];        

    case 4
        load -ascii data/mix1_0825.txt; load -ascii data/mix2_0825.txt; 
        mixed_x=[mix1_0825-mean(mix1_0825); mix2_0825-mean(mix2_0825)];
        
end

