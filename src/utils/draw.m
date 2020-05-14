function draw(varargin)
clf;
hold on
H = {};
for i = 1:nargin
    plot(1:varargin{i}.epoch, varargin{i}.loss(1:varargin{i}.epoch));
    H{i} = varargin{i}.name;
end
legend(H);
end