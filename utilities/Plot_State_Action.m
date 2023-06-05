function Plot_State_Action(State, action, RANGE_domain, RANGE_polygon)

% State = DeleteInvalidPoints(State);

for i =1:size(State,1)
    State(i,:) = AntiNormalize(State(i,:), RANGE_polygon);
end

if ~isempty(action)
    action = AntiNormalize(action, RANGE_polygon);
end

for i = 1:size(State,1)
    plot(State(i,1),State(i,2),'k.','MarkerSize',20)
    hold on;
end

if ~isempty(action)
%     plot(action(1), action(2), 'b+','MarkerSize',20);
end
axis(RANGE_domain);
end