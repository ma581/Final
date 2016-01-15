function y = imshift(x,shift);
y = x;
y(:,:,:,:) = minmin(y); % why more dimensions than the image has

% first shift x:
shx = shift(1);
if shx>0,
	y(1+shx:end,:,:,:) = x(1:end-shx,:,:,:);
elseif shx<0,
	shx = abs(shx);
	y(1:end-shx,:,:,:) = x(1+shx:end,:,:,:);
end;

% then shift y:
shy = shift(2);
if shx~=0,
	x = y;
end;
if shy>0,
	y(:,1+shy:end,:,:) = x(:,1:end-shy,:,:);
elseif shy<0,
	shy = abs(shy);
	y(:,1:end-shy,:,:) = x(:,1+shy:end,:,:);
end;

if max(abs(shift))==0,
	y = x;
end;

end

