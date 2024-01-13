function H_rect=find_effective_H_rect_new(delay_taps,Doppler_taps,chan_coef,M,N)
H_rect=zeros(M*N,M*N);
taps=length(delay_taps);
% for i=0:taps-1
%     l_i=delay_taps(i+1);
%     k_i=Doppler_taps(i+1);
%     T=zeros(M*N,M*N);
% %     for m=0:M-1
% %         for n=0:N-1
% %             p=m+n*M;
% %             q=mod(m-l_i,M) + M*mod(n-k_i,N);
% %             if m<l_i
% %                 term=exp(-1i*2*pi*n/N)*exp(1i*2*pi*k_i*mod(m-l_i,M)/(M*N));
% %             elseif m>=l_i
% %                 term=exp(1i*2*pi*k_i*mod(m-l_i,M)/(M*N));
% %             end
% %             T(p+1,q+1)=term;            
% %         end
% %     end    
%     for p=0:M*N-1
%         for q=0:M*N-1
%             n=floor(p/M);
%             m=p-n*M;
%             q_test=mod(m-l_i,M) + M*mod(n-k_i,N);
%             if q_test==q
%                 if m<l_i
%                     term=exp(-1i*2*pi*n/N)*exp(1i*2*pi*k_i*mod(m-l_i,M)/(M*N));
%                 elseif m>=l_i
%                     term=exp(1i*2*pi*k_i*mod(m-l_i,M)/(M*N));
%                 end
%                 T(p+1,q+1)=term;
%             end
%         end
%     end            
%     %T
%     H_rect=H_rect+h(i+1)*T;
% end
for ele1=1:1:M
        for ele2=1:1:N            
            for tap_no=1:taps                
                if ele1+delay_taps(tap_no)<=M
                    eff_ele1 = ele1 + delay_taps(tap_no);
                    add_term = exp(1i*2*(pi/M)*(ele1-1)*(Doppler_taps(tap_no)/N));
                    int_flag = 0;
                else
                    eff_ele1 = ele1 + delay_taps(tap_no)- M;
                    add_term = exp(1i*2*(pi/M)*(ele1-1-M)*(Doppler_taps(tap_no)/N));
                    int_flag = 1;
                end
                add_term1 = 1;
                if int_flag==1
                    add_term1 = exp(-1i*2*pi*((ele2-1)/N));
                end
                eff_ele2 = mod(ele2-1+Doppler_taps(tap_no),N) + 1;
                new_chan = add_term * add_term1 * chan_coef(tap_no);
                H_rect(N*(eff_ele1-1)+eff_ele2, N*(ele1-1)+ele2)=new_chan;
            end
        end
end
%H_rect=inv(H_rect);
end