function [value, isterminal, direction ] = odeJumpEvent( ~, y)

value = single(~D_inverterLucas(y));
isterminal = 1;
direction = -1;
end

