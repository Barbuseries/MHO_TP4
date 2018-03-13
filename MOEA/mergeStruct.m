function result = mergeStruct(a, b)
  result = a;

  if (~isempty(b))
	f = fieldnames(b);
	for i = 1:length(f)
      result.(f{i}) = b.(f{i});
	end
  end
end
