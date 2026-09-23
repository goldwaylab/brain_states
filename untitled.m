count = [];
count_i = 1;
state = [];
for i = 1:size(partition,1)
    if i == 6237
        start_dynamor = size(count,2);
        disp(start_dynamor);
    elseif i == 23463
        end_dynamor = size(count,2);
        disp(end_dynamor);
    end
    if i==1
        continue;
    end
    if partition(i) == partition(i-1)
        count_i = count_i+1;
    else
        count = [count, count_i];
        state = [state,partition(i-1)];
        count_i = 1;
    end
end

count_dynamor = count(start_dynamor:end_dynamor);
state_dynamor = state(start_dynamor:end_dynamor);
small_states_num = sum(count_dynamor == 1) + sum(count_dynamor == 2)*2;
disp(small_states_num/(23463 - 6237));
state_dynamor_small_1 = state_dynamor(count_dynamor == 1);
state_dynamor_small_2 = state_dynamor(count_dynamor == 2);
state_1 = sum(state_dynamor_small_1==1) + sum(state_dynamor_small_2==1);
state_2 = sum(state_dynamor_small_1==2) + sum(state_dynamor_small_2==2);
state_3 = sum(state_dynamor_small_1==3) + sum(state_dynamor_small_2==3);
state_4 = sum(state_dynamor_small_1==4) + sum(state_dynamor_small_2==4);
state_5 = sum(state_dynamor_small_1==5) + sum(state_dynamor_small_2==5);
state_6 = sum(state_dynamor_small_1==6) + sum(state_dynamor_small_2==6);
disp([state_1, state_2,state_3,state_4,state_5,state_6,])
