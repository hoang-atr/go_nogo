function rate = error2rate(error_lick)

MAX_LICK_COUNT = 15;
NEGATIVE_ERROR_START_RATE = 6;
NEGATIVE_ERROR_END_RATE = 2;
POSITIVE_ERROR_START_RATE = 4;
POSITIVE_ERROR_END_RATE = 8;
BASE_ERROR_RATE = 2;

rate = zeros(size(error_lick));

for i = 1:numel(error_lick)    
    if (error_lick(i) < 0) 
        slope = -(NEGATIVE_ERROR_START_RATE-NEGATIVE_ERROR_END_RATE)/MAX_LICK_COUNT;
        rate(i) = slope * error_lick(i) + NEGATIVE_ERROR_END_RATE;
    elseif (error_lick(i) > 0)     
        slope = (POSITIVE_ERROR_START_RATE-POSITIVE_ERROR_END_RATE)/MAX_LICK_COUNT;
        rate(i) = slope * error_lick(i) + POSITIVE_ERROR_END_RATE;
    else
        rate(i)=2;
    end
end