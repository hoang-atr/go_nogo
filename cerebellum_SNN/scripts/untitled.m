clear
close all

a = 0.3;
b = 16.0;
MAX_LICK_COUNT = 6;
rate = 0:100;
y=MAX_LICK_COUNT - MAX_LICK_COUNT ./ (1+exp(-a.*rate+b));
plot(rate,y)

NEGATIVE_ERROR_START_RATE = 7.0;
NEGATIVE_ERROR_END_RATE = 2.0;
POSITIVE_ERROR_START_RATE = 2.0;
POSITIVE_ERROR_END_RATE = 6.0;
BASE_ERROR_RATE = 2.0;

x = []; y = [];
for error_lick = -MAX_LICK_COUNT:0.1:MAX_LICK_COUNT
    if error_lick < 0
        slope = -(NEGATIVE_ERROR_START_RATE-NEGATIVE_ERROR_END_RATE)/MAX_LICK_COUNT;
        rate = slope*error_lick + NEGATIVE_ERROR_END_RATE;
    elseif error_lick>0
        slope = (POSITIVE_ERROR_START_RATE-POSITIVE_ERROR_END_RATE)/(MAX_LICK_COUNT/2);
        rate = slope*error_lick + POSITIVE_ERROR_END_RATE;
        
%         if (rate<POSITIVE_ERROR_START_RATE) 
%             rate = POSITIVE_ERROR_START_RATE;
%         end
    else
        rate=BASE_ERROR_RATE;
    end
    x = [x; error_lick];
    y = [y; rate];
end

figure
plot(x,y)
xlim([-MAX_LICK_COUNT-0.5 MAX_LICK_COUNT+0.5])
