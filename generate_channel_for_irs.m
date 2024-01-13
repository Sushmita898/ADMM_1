function [delaytaps_irs,dopplertaps_irs,chan_coef_irs] = generate_channel_for_irs(taps,delay_taps1,delay_taps2,Doppler_taps1,Doppler_taps2,chan_coef1,chan_coef2)
delaytaps_irs = [];
           for i = 1:taps
               delaytaps1 = [];
               for j = 1:taps
                   delaytaps2 = delay_taps1(i)+delay_taps2(j);
                   delaytaps1 = [delaytaps1 delaytaps2];
               end
               delaytaps_irs = [delaytaps_irs delaytaps1];
           end
           % generating doppler taps for the cascaded channel
            dopplertaps_irs = [];
           for i = 1:taps
               dopplertaps1 = [];
               for j = 1:taps
                   dopplertaps2 = Doppler_taps1(i)+Doppler_taps2(j);
                   dopplertaps1 = [dopplertaps1 dopplertaps2];
               end
               dopplertaps_irs = [dopplertaps_irs dopplertaps1];
           end
           % generating channel_coef for cascaded channel
           
           chan_coef_irs = [];
           for i = 1:taps
               chan_coef_1 = [];
               for j = 1:taps
                   chan_coef_2 = chan_coef1(i)*chan_coef2(j);
                   chan_coef_1 = [chan_coef_1 chan_coef_2];
               end
               chan_coef_irs = [chan_coef_irs chan_coef_1];
           end
end