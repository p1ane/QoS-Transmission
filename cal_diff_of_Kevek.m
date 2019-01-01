function [diff_of_Kevek] = cal_diff_of_Kevek(CCMev,Lambda,Nt,Nr)
%CAL_DIFF_OF_Keve
%Lambda带入单个用户的Lambda，即Lambda(:,:,k)
%返回的变量也只是一个用户的求导，diff_of_Keve(:,:,k)
%checked had been wrong!
temp_b = zeros(Nt,1);
temp_a = zeros(Nr,1);
diff_of_Kevek_ele = zeros(Nt,1);


for m = 1:Nt
    temp_b = zeros(Nt,1);%changed!
    
    
    for i = 1:Nr
        temp_a = zeros(Nr,1);
        
        for j = 1:Nt
            temp_a(i) = temp_a(i) + CCMev(i,j)*Lambda(j,j);
        end
        
        
        
        temp_b(m) = temp_b(m) + CCMev(i,m)/(1+temp_a(i));
    end
    
    
    diff_of_Kevek_ele(m) = temp_b(m);
end


diff_of_Kevek = diag(diff_of_Kevek_ele);
end

