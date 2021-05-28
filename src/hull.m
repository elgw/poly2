function hull()
close all

% Melkman
P=[[101560.000000, 102289.000000]; [101561.000000, 102290.000000]; [101564.000000, 102290.000000]; [101565.000000, 102291.000000]; [101567.000000, 102291.000000]; [101569.000000, 102293.000000]; [101569.000000, 102294.000000]; [101570.000000, 102295.000000]; [101570.000000, 102299.000000]; [101568.000000, 102301.000000]; [101565.000000, 102301.000000]; [101564.000000, 102302.000000]; [101563.000000, 102301.000000]; [101563.000000, 102300.000000]; [101562.000000, 102299.000000]; [101559.000000, 102299.000000]; [101558.000000, 102300.000000]; [101558.000000, 102301.000000]; [101559.000000, 102302.000000]; [101559.000000, 102304.000000]; [101558.000000, 102305.000000]; [101555.000000, 102305.000000]; [101553.000000, 102303.000000]; [101553.000000, 102301.000000]; [101551.000000, 102299.000000]; [101551.000000, 102294.000000]; [101552.000000, 102293.000000]; [101552.000000, 102292.000000]; [101553.000000, 102291.000000]; [101553.000000, 102290.000000]; [101555.000000, 102290.000000]; [101556.000000, 102289.000000]]
%P = cat(1, P, P(1,:));
plotP(P);

H = melkman(P);
figure,
plotP(P);
plotP(H);

end

function H = melkman(P)
n = size(P,1);

% An array implementation of a dequeue
dq = -1*ones(2*n, 1);
% db and dt will always refer to the same vertex
b = n; % Bottom
t = b-1; % Top 

v1 = 1; v2 = 2; v3 = 3;
if rcl(P(v1, :), P(v2, :), P(v3, :)) > 0
    [dq, b, t] = push(dq, b, t, v1);
    [dq, b, t] = push(dq, b, t, v2);
    % push v1, push v2
else
    % push v2, push v1
    [dq, b, t] = push(dq, b, t, v2);
    [dq, b, t] = push(dq, b, t, v1);
end
% push v3
[dq, b, t] = push(dq, b, t, v3);
% insert v3
[dq, b, t] = insert(dq, b, t, v3);


next = 4;

while next < n
    v = next; next = next+1;
    %% 2
while ~( rcl( P(v, :),P(dq(b),:),P(dq(b+1),:)) < 0 || ...
    rcl(P(dq(t-1), :), P(dq(t), :), P(v,:)) < 0)
    if next == n+1
        H = P(dq(b:t), :);
        return;        
    end
    v = next;
    next = next+1;    
end

%% 3
while ~(rcl(P(dq(t-1), :), P(dq(t), :), P(v, :)) > 0)
    t = t-1; %pop
end

%% 4
% push v
[dq, b, t] = push(dq, b, t, v);

while ~(rcl(P(v, :), P(dq(b), :), P(dq(b+1),:)) > 0)
    % remove d_b
    b = b+1;
end
[dq, b, t] = insert(dq, b, t, v);

plotP(P(dq(b:t), :))
end

H = P(dq(b:t), :);
return;
end

function [dq, b, t] = insert(dq, b, t, e)
b = b - 1;
dq(b) = e;
end

function [dq, b, t] = push(dq, b, t, e)
t = t+1;
dq(t) = e;
end

function v = rcl(A, B, C)
% right = 1;
% colinear = 0;
% left = -1;

Q = B - A;
R = C - A;
v = (Q(1)*R(2) - Q(2)*R(1));

end

function plotP(P)
P = cat(1, P, P(1,:));
hold on
plot(P(:,1), P(:,2))
end