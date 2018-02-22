function Clamp
								%CLAMP All clamp functions.
								%
								% default
								% fancy
								%
								% See also CLAMP>DEFAULT, CLAMP>FANCY
  
  global CLAMP;

  CLAMP.default = @default;
  CLAMP.fancy = @fancy;
end

function result = default(val, lowest, biggest)
							   %DEFAULT If VAL < LOWEST, RESULT = LOWEST.
							   % If VAL > BIGGEST, RESULT = BIGGEST.
  
  result = max(min(val, biggest), lowest);
end

%% NOTE: Btw, this fails if val is too low or too high.
function result = fancy(val, lowest, biggest)
  %FANCY If VAL < LOWEST, RESULT = 2 * LOWEST - VAL.
  % If VAL > BIGGEST, RESULT = 2 * BIGGEST - VAL.
  
  below = val < lowest;
  above = val > biggest;
  in_bounds =  ~below & ~above;

  result = val .* in_bounds + (2 * lowest - val) .* below + (2 * biggest - val) .* above;
end
